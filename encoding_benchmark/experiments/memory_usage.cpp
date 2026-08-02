#include "benchmark_adapter.h"
#include "experiment_api.h"

void run_memory_usage(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "memory");
    begin_benchmark_output(output, "memory", true);
    for (const DatasetConfig& dataset : context.datasets) {
        run_baseline_dataset(dataset, dataset_path(context, dataset), output, false, true);
    }
    finish_benchmark_output(output, true);
}
