#include "benchmark_adapter.h"
#include "experiment_api.h"

void run_kmeans_comparison(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "kmeans");
    begin_benchmark_output(output, "kmeans", false);
    for (const DatasetConfig& dataset : context.datasets) {
        run_kmeans_dataset(dataset, dataset_path(context, dataset), output);
    }
    finish_benchmark_output(output, false);
}
