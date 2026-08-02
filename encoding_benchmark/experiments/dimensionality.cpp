#include "benchmark_adapter.h"
#include "experiment_api.h"

#include <filesystem>

void run_dimensionality_experiment(const ExperimentContext& context) {
    std::vector<std::string> names;
    std::vector<int> dims;
    for (const DatasetConfig& dataset : context.datasets) {
        if (dataset.dim <= 1) continue;
        names.push_back(dataset.name);
        dims.push_back(dataset.dim);
    }
    if (names.empty()) return;
    const std::string output = experiment_output(context, "dimensionality");
    std::filesystem::create_directories(output);
    std::string root = context.data_root;
    if (!root.empty() && root.back() != '/' && root.back() != '\\') root.push_back('/');
    run_dimensionality_core(root, output, names, dims);
}
