#include "experiment_api.h"
#include "standalone_support.h"

#include "algorithms/acluster.h"
#include "algorithms/kcluster.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace {

template<typename Compressor>
void write_components(
    std::ofstream& csv,
    const std::string& method,
    const DatasetConfig& dataset,
    const std::string& path,
    int runs
) {
    const auto original = read_csv_rows(path);
    double preprocessing = 0.0;
    double clustering = 0.0;
    double residual = 0.0;
    double bitstream = 0.0;
    double header = 0.0;
    double medoids = 0.0;
    double frequencies = 0.0;
    double residuals = 0.0;

    for (int run = 0; run < runs; ++run) {
        Compressor compressor;
        int mutable_dim = dataset.dim;
        auto encoded = compressor.encode_multidim_with_diag(
            path, mutable_dim, 100, 10000, 10, 10
        );
        if (run == 0) verify_rows(original, compressor.decode_multidim(encoded.first));
        preprocessing += encoded.second.total_timings.preprocessing_ns;
        clustering += encoded.second.total_timings.clustering_ns;
        residual += encoded.second.total_timings.residual_calc_ns;
        bitstream += encoded.second.total_timings.bitstream_gen_ns;
        header += encoded.second.total_sizes.header_bits;
        medoids += encoded.second.total_sizes.medoids_bits;
        frequencies += encoded.second.total_sizes.frequencies_bits;
        residuals += encoded.second.total_sizes.residuals_bits;
    }

    const double scale = static_cast<double>(runs) * original.size();
    csv << dataset.name << ',' << method << ',' << runs << ','
        << std::fixed << std::setprecision(4)
        << preprocessing / scale << ',' << clustering / scale << ','
        << residual / scale << ',' << bitstream / scale << ','
        << header / scale << ',' << medoids / scale << ','
        << frequencies / scale << ',' << residuals / scale << '\n';
}

}

void run_encoding_components(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "components");
    std::filesystem::create_directories(output);
    std::ofstream csv(output + "/results.csv");
    csv << "Dataset,Method,Runs,PreprocessingNsPerRow,ClusteringNsPerRow,ResidualNsPerRow,"
           "BitstreamNsPerRow,HeaderBitsPerRow,ReferenceBitsPerRow,FrequencyBitsPerRow,ResidualBitsPerRow\n";
    for (const DatasetConfig& dataset : context.datasets) {
        const std::string path = dataset_path(context, dataset);
        write_components<KCluster>(csv, "KCluster", dataset, path, context.runs);
        write_components<ACluster>(csv, "ACluster", dataset, path, context.runs);
    }
}
