#include "experiments/benchmark_adapter.h"
#include "experiments/experiment_api.h"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char* argv[]) {
    try {
        ExperimentContext context{"../data/final", "output", public_datasets(), 100};
        std::string experiment = "all";
        std::string dataset_filter;

        for (int i = 1; i < argc; ++i) {
            const std::string option = argv[i];
            if (option == "--experiment" && i + 1 < argc) experiment = argv[++i];
            else if (option == "--runs" && i + 1 < argc) context.runs = std::stoi(argv[++i]);
            else if (option == "--dataset" && i + 1 < argc) dataset_filter = argv[++i];
            else if (option == "--data-root" && i + 1 < argc) context.data_root = argv[++i];
            else if (option == "--output-root" && i + 1 < argc) context.output_root = argv[++i];
            else throw std::runtime_error("Unknown or incomplete option: " + option);
        }

        if (context.runs <= 0) throw std::runtime_error("--runs must be positive");
        if (!dataset_filter.empty()) {
            context.datasets.erase(
                std::remove_if(
                    context.datasets.begin(), context.datasets.end(),
                    [&](const DatasetConfig& dataset) { return dataset.name != dataset_filter; }
                ),
                context.datasets.end()
            );
            if (context.datasets.empty()) throw std::runtime_error("Unknown dataset: " + dataset_filter);
        }

        set_benchmark_runs(context.runs);
        auto selected = [&](const std::string& name) {
            return experiment == "all" || experiment == name;
        };

        if (selected("rebuttal-1d")) run_rebuttal_1d(context);
        if (selected("baseline")) run_baseline_comparison(context);
        if (selected("gzip")) run_gzip_comparison(context);
        if (selected("order")) run_order_test(context);
        if (selected("components")) run_encoding_components(context);
        if (selected("memory")) run_memory_usage(context);
        if (selected("outlier")) run_outlier_experiment(context);
        if (selected("kmeans")) run_kmeans_comparison(context);
        if (selected("sensitivity")) run_parameter_sensitivity(context);
        if (selected("reference")) run_reference_comparison(context);
        if (selected("dimensionality")) run_dimensionality_experiment(context);
        if (selected("subspace")) run_subspace_grouping(context);

        const std::vector<std::string> valid = {
            "all", "rebuttal-1d", "baseline", "gzip", "order", "components",
            "memory", "outlier", "kmeans", "sensitivity", "reference",
            "dimensionality", "subspace"
        };
        if (std::find(valid.begin(), valid.end(), experiment) == valid.end()) {
            throw std::runtime_error("Unknown experiment: " + experiment);
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
