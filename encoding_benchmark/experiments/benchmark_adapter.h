#pragma once

#include "experiment_api.h"

#include <string>
#include <vector>

void set_benchmark_runs(int runs);
void begin_benchmark_output(const std::string& output_dir, const std::string& scenario, bool memory);
void finish_benchmark_output(const std::string& output_dir, bool memory);
void run_baseline_dataset(
    const DatasetConfig& dataset,
    const std::string& input_path,
    const std::string& output_dir,
    bool gzip,
    bool include_vortex
);
void run_kmeans_dataset(
    const DatasetConfig& dataset,
    const std::string& input_path,
    const std::string& output_dir
);
void run_reference_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
);
void run_dimensionality_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
);
void run_subspace_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims,
    int group_size
);
