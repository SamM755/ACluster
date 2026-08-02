#include "cluster_common.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <unordered_set>
#include <vector>

namespace ClusterLogic {

namespace {

struct VectorHash {
    std::size_t operator()(const std::vector<long long>& values) const {
        std::size_t seed = 0;
        for (long long value : values) {
            seed ^= std::hash<long long>{}(value) + 0x9e3779b9U + (seed << 6U) + (seed >> 2U);
        }
        return seed;
    }
};

std::vector<std::vector<long long>> to_rows(const SoAData<long long>& data) {
    std::vector<std::vector<long long>> rows(
        data.num_rows, std::vector<long long>(data.dim, 0)
    );
    for (int column = 0; column < data.dim; ++column) {
        for (int row = 0; row < data.num_rows; ++row) rows[row][column] = data.columns[column][row];
    }
    return rows;
}

ClusteringResult sorted_result(
    const std::vector<std::vector<long long>>& centroids,
    const std::vector<int>& assignment
) {
    std::vector<long long> original_sizes(centroids.size(), 0);
    for (int cluster : assignment) {
        if (cluster >= 0) ++original_sizes[cluster];
    }
    std::vector<int> order(centroids.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int left, int right) {
        return original_sizes[left] < original_sizes[right];
    });
    std::vector<int> remap(centroids.size());
    std::vector<std::vector<long long>> sorted_centroids;
    std::vector<long long> sorted_sizes;
    sorted_centroids.reserve(centroids.size());
    sorted_sizes.reserve(centroids.size());
    for (std::size_t index = 0; index < order.size(); ++index) {
        remap[order[index]] = static_cast<int>(index);
        sorted_centroids.push_back(centroids[order[index]]);
        sorted_sizes.push_back(original_sizes[order[index]]);
    }
    std::vector<int> sorted_assignment = assignment;
    for (int& cluster : sorted_assignment) {
        if (cluster >= 0) cluster = remap[cluster];
    }
    return {std::move(sorted_centroids), std::move(sorted_assignment), std::move(sorted_sizes)};
}

}

ClusteringResult kMeansEuclidean(
    const SoAData<long long>& data,
    int k,
    int max_iter,
    double tolerance,
    int dim
) {
    if (data.num_rows == 0) return {{}, {}, {}};
    auto rows = to_rows(data);
    auto distinct = rows;
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    k = std::min(k, static_cast<int>(distinct.size()));
    if (k <= 0) return {{}, {}, {}};

    std::mt19937 generator(20260406U);
    std::vector<int> indices(data.num_rows);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), generator);
    std::unordered_set<std::vector<long long>, VectorHash> seen;
    std::vector<std::vector<long long>> centroids;
    for (int index : indices) {
        if (seen.insert(rows[index]).second) centroids.push_back(rows[index]);
        if (static_cast<int>(centroids.size()) == k) break;
    }
    k = static_cast<int>(centroids.size());
    std::vector<int> assignment(data.num_rows, -1);

    for (int iteration = 0; iteration < max_iter; ++iteration) {
        for (int row = 0; row < data.num_rows; ++row) {
            long double best_distance = std::numeric_limits<long double>::infinity();
            int best_cluster = 0;
            for (int cluster = 0; cluster < k; ++cluster) {
                long double distance = 0.0L;
                for (int column = 0; column < dim; ++column) {
                    const long double difference = static_cast<long double>(rows[row][column])
                        - static_cast<long double>(centroids[cluster][column]);
                    distance += difference * difference;
                }
                if (distance < best_distance) {
                    best_distance = distance;
                    best_cluster = cluster;
                }
            }
            assignment[row] = best_cluster;
        }

        std::vector<std::vector<long double>> sums(k, std::vector<long double>(dim, 0.0L));
        std::vector<int> counts(k, 0);
        for (int row = 0; row < data.num_rows; ++row) {
            const int cluster = assignment[row];
            ++counts[cluster];
            for (int column = 0; column < dim; ++column) sums[cluster][column] += rows[row][column];
        }

        long double shift = 0.0L;
        auto updated = centroids;
        for (int cluster = 0; cluster < k; ++cluster) {
            if (counts[cluster] == 0) continue;
            for (int column = 0; column < dim; ++column) {
                const long long value = static_cast<long long>(
                    std::llround(sums[cluster][column] / counts[cluster])
                );
                shift += std::abs(static_cast<long double>(value - centroids[cluster][column]));
                updated[cluster][column] = value;
            }
        }
        centroids.swap(updated);
        if (shift <= tolerance) break;
    }
    return sorted_result(centroids, assignment);
}

}
