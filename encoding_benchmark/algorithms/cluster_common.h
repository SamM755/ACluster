// in algorithms/cluster_common.h
#ifndef CLUSTER_COMMON_H
#define CLUSTER_COMMON_H

#include <vector>
#include <string>
#include <cstdint>
#include <memory>
#include <iostream>
#include <iomanip> 
#include <chrono>

template<typename T>
struct SoAData {
    std::vector<std::vector<T>> columns;
    int dim = 0;
    int num_rows = 0;
};

struct ClusteringResult {
    std::vector<std::vector<long long>> medoids;
    std::vector<int> cluster_assignment;
    std::vector<long long> cluster_sizes;
};

struct PreprocessorResult {
    std::unique_ptr<SoAData<long long>> data;
    std::vector<int> max_decimal_places;
    std::vector<long long> min_values;
};


struct StageTimings {
    long long preprocessing_ns = 0;
    long long clustering_ns = 0;
    long long residual_calc_ns = 0;
    long long bitstream_gen_ns = 0;
};


struct ComponentSizes {
    size_t header_bits = 0;
    size_t medoids_bits = 0;
    size_t frequencies_bits = 0;
    size_t residuals_bits = 0;
};

struct KClusterCosts {
    long long final_residual_cost_bits = 0; 
};

struct AClusterCosts {
    long long final_total_cost_bits = 0; 
};

struct PageDiagnostics {
    int page_index;
    StageTimings timings;
    ComponentSizes sizes;
    
    KClusterCosts kcluster_costs;
    AClusterCosts acluster_costs;

    std::string algorithm_name; 
};

struct GlobalDiagnostics {
    StageTimings total_timings;
    ComponentSizes total_sizes;
    

    void print_summary(size_t num_pages) {
        long long total_time = total_timings.preprocessing_ns + total_timings.clustering_ns +
                               total_timings.residual_calc_ns + total_timings.bitstream_gen_ns;
        if (total_time == 0) total_time = 1;

        std::cout << "\n--- Performance Diagnostics Summary (" << num_pages << " pages) ---\n";
        std::cout << "1. Stage Timing Breakdown:\n";
        std::cout << "   - Preprocessing:  " << std::setw(8) << total_timings.preprocessing_ns / 1e6 << " ms ("
                  << (total_timings.preprocessing_ns * 100.0 / total_time) << "%)\n";
        std::cout << "   - Clustering:     " << std::setw(8) << total_timings.clustering_ns / 1e6 << " ms ("
                  << (total_timings.clustering_ns * 100.0 / total_time) << "%)\n";
        std::cout << "   - Residual Calc:  " << std::setw(8) << total_timings.residual_calc_ns / 1e6 << " ms ("
                  << (total_timings.residual_calc_ns * 100.0 / total_time) << "%)\n";
        std::cout << "   - Bitstream Gen:  " << std::setw(8) << total_timings.bitstream_gen_ns / 1e6 << " ms ("
                  << (total_timings.bitstream_gen_ns * 100.0 / total_time) << "%)\n";

        size_t total_bits = total_sizes.header_bits + total_sizes.medoids_bits +
                            total_sizes.frequencies_bits + total_sizes.residuals_bits;
        if (total_bits == 0) total_bits = 1;
        
        std::cout << "2. Encoding Component Size Breakdown:\n";
        std::cout << "   - Header:         " << std::setw(8) << total_sizes.header_bits / 8 << " bytes ("
                  << (total_sizes.header_bits * 100.0 / total_bits) << "%)\n";
        std::cout << "   - Medoids:        " << std::setw(8) << total_sizes.medoids_bits / 8 << " bytes ("
                  << (total_sizes.medoids_bits * 100.0 / total_bits) << "%)\n";
        std::cout << "   - Frequencies:    " << std::setw(8) << total_sizes.frequencies_bits / 8 << " bytes ("
                  << (total_sizes.frequencies_bits * 100.0 / total_bits) << "%)\n";
        std::cout << "   - Residuals:      " << std::setw(8) << total_sizes.residuals_bits / 8 << " bytes ("
                  << (total_sizes.residuals_bits * 100.0 / total_bits) << "%)\n";
        std::cout << "--------------------------------------------------------\n";
    }
};

#endif // CLUSTER_COMMON_H