#include "cluster_case.h"
#include "experiment_api.h"

#include <filesystem>
#include <fstream>
#include <iomanip>

void run_reference_comparison(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "reference");
    std::filesystem::create_directories(output);
    std::ofstream csv(output + "/results.csv");
    csv << "Dataset,Method,K,Runs,CompressedBytes,CompressionRatio,EncodeNsPerRow,DecodeNsPerRow,References\n";

    for (const DatasetConfig& dataset : context.datasets) {
        const std::string path = dataset_path(context, dataset);
        const ClusterCaseResult adaptive = benchmark_cluster_case<ACluster>(
            path, dataset.dim, 100, 10000, 10, context.runs
        );
        csv << dataset.name << ",ACluster,adaptive," << context.runs << ','
            << std::fixed << std::setprecision(2) << adaptive.compressed_bytes << ','
            << std::setprecision(4) << adaptive.compression_ratio << ','
            << std::setprecision(2) << adaptive.encode_ns_per_row << ','
            << adaptive.decode_ns_per_row << ',' << adaptive.references << '\n';
        for (int k : {100, 200, 500}) {
            const ClusterCaseResult fixed = benchmark_cluster_case<KCluster>(
                path, dataset.dim, k, 10000, 10, context.runs
            );
            csv << dataset.name << ",KCluster," << k << ',' << context.runs << ','
                << std::fixed << std::setprecision(2) << fixed.compressed_bytes << ','
                << std::setprecision(4) << fixed.compression_ratio << ','
                << std::setprecision(2) << fixed.encode_ns_per_row << ','
                << fixed.decode_ns_per_row << ',' << fixed.references << '\n';
        }
    }
}
