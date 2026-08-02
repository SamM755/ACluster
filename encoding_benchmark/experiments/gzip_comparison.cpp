#include "benchmark_adapter.h"
#include "experiment_api.h"

void run_gzip_comparison(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "gzip");
    begin_benchmark_output(output, "gzip", false);
    for (const DatasetConfig& dataset : context.datasets) {
        run_baseline_dataset(dataset, dataset_path(context, dataset), output, true, false);
    }
    finish_benchmark_output(output, false);
}
