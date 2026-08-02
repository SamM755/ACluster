#pragma once

#include <string>
#include <vector>

struct DatasetConfig {
    std::string name;
    std::string file_name;
    int dim;
    int rows;
};

struct ExperimentContext {
    std::string data_root;
    std::string output_root;
    std::vector<DatasetConfig> datasets;
    int runs;
};

std::vector<DatasetConfig> public_datasets();
std::string dataset_path(const ExperimentContext& context, const DatasetConfig& dataset);
std::string experiment_output(const ExperimentContext& context, const std::string& name);

void run_rebuttal_1d(const ExperimentContext& context);
void run_baseline_comparison(const ExperimentContext& context);
void run_gzip_comparison(const ExperimentContext& context);
void run_order_test(const ExperimentContext& context);
void run_encoding_components(const ExperimentContext& context);
void run_memory_usage(const ExperimentContext& context);
void run_outlier_experiment(const ExperimentContext& context);
void run_kmeans_comparison(const ExperimentContext& context);
void run_parameter_sensitivity(const ExperimentContext& context);
void run_reference_comparison(const ExperimentContext& context);
void run_dimensionality_experiment(const ExperimentContext& context);
void run_subspace_grouping(const ExperimentContext& context);
