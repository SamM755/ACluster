#pragma once

#include "algorithms/acluster.h"
#include "algorithms/kcluster.h"
#include "standalone_support.h"

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

struct ClusterCaseResult {
    double compressed_bytes;
    double compression_ratio;
    double encode_ns_per_row;
    double decode_ns_per_row;
    double references;
};

template<typename Compressor>
ClusterCaseResult benchmark_cluster_case(
    const std::string& path,
    int dim,
    int k,
    int page_size,
    int pack_size,
    int runs
) {
    const auto original = read_csv_rows(path);
    if (original.empty()) throw std::runtime_error("Empty dataset: " + path);
    double size_sum = 0.0;
    double encode_sum = 0.0;
    double decode_sum = 0.0;
    double reference_sum = 0.0;

    for (int run = 0; run < runs; ++run) {
        Compressor compressor;
        int mutable_dim = dim;
        std::vector<uint8_t> encoded;
        long long references = 0;
        const auto encode_start = std::chrono::high_resolution_clock::now();
        if constexpr (std::is_same_v<Compressor, KCluster> || std::is_same_v<Compressor, ACluster>) {
            auto result = compressor.encode_multidim_with_diag(
                path, mutable_dim, k, page_size, pack_size, 10
            );
            encoded = std::move(result.first);
            references = result.second.total_reference_points;
        } else {
            encoded = compressor.encode_multidim(path, mutable_dim, k, page_size, pack_size, 10);
        }
        const auto encode_end = std::chrono::high_resolution_clock::now();
        const auto decode_start = std::chrono::high_resolution_clock::now();
        auto decoded = compressor.decode_multidim(encoded);
        const auto decode_end = std::chrono::high_resolution_clock::now();
        if (run == 0) verify_rows(original, decoded);
        size_sum += static_cast<double>(encoded.size());
        encode_sum += std::chrono::duration_cast<std::chrono::nanoseconds>(encode_end - encode_start).count();
        decode_sum += std::chrono::duration_cast<std::chrono::nanoseconds>(decode_end - decode_start).count();
        reference_sum += static_cast<double>(references);
    }

    const double run_count = static_cast<double>(runs);
    const double row_count = static_cast<double>(original.size());
    const double average_size = size_sum / run_count;
    return {
        average_size,
        row_count * dim * sizeof(double) / average_size,
        encode_sum / run_count / row_count,
        decode_sum / run_count / row_count,
        reference_sum / run_count,
    };
}
