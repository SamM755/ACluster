#include "cluster_case.h"
#include "experiment_api.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace {

void write_case(
    std::ofstream& csv,
    const DatasetConfig& dataset,
    const std::string& experiment,
    const std::string& method,
    int value,
    int runs,
    const ClusterCaseResult& result
) {
    csv << dataset.name << ',' << experiment << ',' << method << ',' << value << ',' << runs << ','
        << std::fixed << std::setprecision(2) << result.compressed_bytes << ','
        << std::setprecision(4) << result.compression_ratio << ','
        << std::setprecision(2) << result.encode_ns_per_row << ','
        << result.decode_ns_per_row << ',' << result.references << '\n';
}

}

void run_parameter_sensitivity(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "sensitivity");
    std::filesystem::create_directories(output);
    std::ofstream csv(output + "/results.csv");
    csv << "Dataset,Experiment,Method,Value,Runs,CompressedBytes,CompressionRatio,"
           "EncodeNsPerRow,DecodeNsPerRow,References\n";

    for (const DatasetConfig& dataset : context.datasets) {
        const std::string path = dataset_path(context, dataset);
        for (int k : {10, 100, 1000}) {
            write_case(
                csv, dataset, "k", "KCluster", k, context.runs,
                benchmark_cluster_case<KCluster>(path, dataset.dim, k, 10000, 10, context.runs)
            );
        }
        for (int page_size : {1000, 10000, 20000}) {
            write_case(
                csv, dataset, "page", "KCluster", page_size, context.runs,
                benchmark_cluster_case<KCluster>(path, dataset.dim, 100, page_size, 10, context.runs)
            );
            write_case(
                csv, dataset, "page", "ACluster", page_size, context.runs,
                benchmark_cluster_case<ACluster>(path, dataset.dim, 100, page_size, 10, context.runs)
            );
        }
        for (int pack_size : {1, 10, 100}) {
            write_case(
                csv, dataset, "pack", "KCluster", pack_size, context.runs,
                benchmark_cluster_case<KCluster>(path, dataset.dim, 100, 10000, pack_size, context.runs)
            );
            write_case(
                csv, dataset, "pack", "ACluster", pack_size, context.runs,
                benchmark_cluster_case<ACluster>(path, dataset.dim, 100, 10000, pack_size, context.runs)
            );
        }
    }
}
