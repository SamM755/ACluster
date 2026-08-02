#include "experiment_api.h"
#include "cluster_case.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>

void run_rebuttal_1d(const ExperimentContext& context) {
    const std::map<std::string, std::size_t> expected = {
        {"SSD", 3397}, {"WS", 51929}, {"PI", 4205}, {"BC", 91873}
    };
    const std::string output = experiment_output(context, "rebuttal_1d");
    std::filesystem::create_directories(output);
    std::ofstream csv(output + "/results.csv");
    csv << "Dataset,Runs,CompressedBytes,CompressionRatio,EncodeNsPerRow,DecodeNsPerRow,References\n";

    for (const DatasetConfig& dataset : context.datasets) {
        const auto match = expected.find(dataset.name);
        if (match == expected.end()) continue;
        const ClusterCaseResult result = benchmark_cluster_case<ACluster>(
            dataset_path(context, dataset), 1, 100, 10000, 10, context.runs
        );
        if (static_cast<std::size_t>(result.compressed_bytes + 0.5) != match->second) {
            throw std::runtime_error("Rebuttal byte regression failed for " + dataset.name);
        }
        csv << dataset.name << ',' << context.runs << ','
            << std::fixed << std::setprecision(2) << result.compressed_bytes << ','
            << std::setprecision(4) << result.compression_ratio << ','
            << std::setprecision(2) << result.encode_ns_per_row << ','
            << result.decode_ns_per_row << ',' << result.references << '\n';
    }
}
