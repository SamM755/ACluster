#include "vortex.h"

#include "../utils/type_converter.h"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace {

using RankColPair = std::pair<uint32_t, uint16_t>;

std::vector<std::unordered_map<uint64_t, uint32_t>> build_rank_maps(const std::vector<std::vector<double>>& rows) {
    const size_t num_rows = rows.size();
    const size_t dim = rows[0].size();

    std::vector<std::unordered_map<uint64_t, uint32_t>> value_counts(dim);
    for (size_t j = 0; j < dim; ++j) {
        value_counts[j].reserve(num_rows * 2);
    }

    for (const auto& row : rows) {
        if (row.size() != dim) {
            throw std::invalid_argument("Vortex requires all rows to have the same dimension.");
        }
        for (size_t j = 0; j < dim; ++j) {
            uint64_t bits = double_to_long(row[j]);
            value_counts[j][bits] += 1;
        }
    }

    std::vector<std::unordered_map<uint64_t, uint32_t>> rank_maps(dim);
    for (size_t j = 0; j < dim; ++j) {
        std::vector<std::pair<uint64_t, uint32_t>> freq;
        freq.reserve(value_counts[j].size());
        for (const auto& kv : value_counts[j]) {
            freq.push_back(kv);
        }

        std::sort(freq.begin(), freq.end(), [](const auto& a, const auto& b) {
            if (a.second != b.second) return a.second > b.second;
            return a.first < b.first;
        });

        auto& ranks = rank_maps[j];
        ranks.reserve(freq.size() * 2);
        for (uint32_t rank = 0; rank < static_cast<uint32_t>(freq.size()); ++rank) {
            ranks[freq[rank].first] = rank + 1;
        }
    }

    return rank_maps;
}

std::vector<std::vector<RankColPair>> build_vortex_keys(
    const std::vector<std::vector<double>>& rows,
    const std::vector<std::unordered_map<uint64_t, uint32_t>>& rank_maps
) {
    const size_t num_rows = rows.size();
    const size_t dim = rows[0].size();
    if (dim > static_cast<size_t>(std::numeric_limits<uint16_t>::max())) {
        throw std::invalid_argument("Vortex dimension exceeds uint16_t range.");
    }

    std::vector<std::vector<RankColPair>> keys(num_rows, std::vector<RankColPair>(dim));
    for (size_t i = 0; i < num_rows; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            uint64_t bits = double_to_long(rows[i][j]);
            auto it = rank_maps[j].find(bits);
            if (it == rank_maps[j].end()) {
                throw std::runtime_error("Vortex rank lookup failed.");
            }
            keys[i][j] = {it->second, static_cast<uint16_t>(j)};
        }

        std::sort(keys[i].begin(), keys[i].end(), [](const RankColPair& a, const RankColPair& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });
    }

    return keys;
}

bool vortex_less(const std::vector<RankColPair>& lhs, const std::vector<RankColPair>& rhs) {
    const size_t n = lhs.size();
    for (size_t i = 0; i < n; ++i) {
        if (lhs[i] == rhs[i]) continue;
        if ((i & 1U) == 0U) {
            return lhs[i] < rhs[i];
        }
        return rhs[i] < lhs[i];
    }
    return false;
}

} 

std::vector<uint32_t> Vortex::compute_permutation(const std::vector<std::vector<double>>& rows) const {
    if (rows.empty()) return {};
    if (rows.size() > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
        throw std::invalid_argument("Vortex row count exceeds uint32_t range.");
    }
    if (rows[0].empty()) {
        std::vector<uint32_t> identity(rows.size());
        std::iota(identity.begin(), identity.end(), 0U);
        return identity;
    }

    const auto rank_maps = build_rank_maps(rows);
    const auto keys = build_vortex_keys(rows, rank_maps);

    std::vector<uint32_t> permutation(rows.size());
    std::iota(permutation.begin(), permutation.end(), 0U);

    std::stable_sort(permutation.begin(), permutation.end(), [&](uint32_t a, uint32_t b) {
        if (vortex_less(keys[a], keys[b])) return true;
        if (vortex_less(keys[b], keys[a])) return false;
        return a < b;
    });

    return permutation;
}

std::vector<std::vector<double>> Vortex::reorder_rows(const std::vector<std::vector<double>>& rows) const {
    if (rows.empty()) return {};
    auto permutation = compute_permutation(rows);
    std::vector<std::vector<double>> reordered;
    reordered.reserve(rows.size());
    for (uint32_t idx : permutation) {
        reordered.push_back(rows[idx]);
    }
    return reordered;
}
