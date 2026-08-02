#include "experiment_api.h"
#include "standalone_support.h"

#include "algorithms/acluster.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>

void run_outlier_experiment(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "outlier");
    const std::string cache = output + "/cache";
    std::filesystem::create_directories(cache);
    std::ofstream csv(output + "/results.csv");
    csv << "Dataset,OutlierRatio,Runs,CompressedBytes,CompressionRatio,EncodeNsPerRow,"
           "DecodeNsPerRow,References\n";

    for (const DatasetConfig& dataset : context.datasets) {
        auto base_rows = read_csv_rows_as_strings(dataset_path(context, dataset));
        if (base_rows.size() > 10000) base_rows.resize(10000);
        if (base_rows.empty()) continue;
        std::vector<double> minimum(dataset.dim, std::numeric_limits<double>::infinity());
        std::vector<double> maximum(dataset.dim, -std::numeric_limits<double>::infinity());
        for (const auto& row : base_rows) {
            for (int column = 0; column < dataset.dim; ++column) {
                const double value = column < static_cast<int>(row.size()) ? parse_number(row[column]) : 0.0;
                minimum[column] = std::min(minimum[column], value);
                maximum[column] = std::max(maximum[column], value);
            }
        }

        for (double ratio : {0.0, 0.001, 0.005, 0.01, 0.05}) {
            double size_sum = 0.0;
            double encode_sum = 0.0;
            double decode_sum = 0.0;
            double reference_sum = 0.0;
            for (int run = 0; run < context.runs; ++run) {
                auto rows = base_rows;
                const int outlier_count = static_cast<int>(std::round(rows.size() * ratio));
                std::vector<int> indices(rows.size());
                std::iota(indices.begin(), indices.end(), 0);
                const std::size_t seed = 20260406U + std::hash<std::string>{}(dataset.name)
                    + static_cast<std::size_t>(ratio * 1000000.0) + static_cast<std::size_t>(run);
                std::mt19937 generator(static_cast<unsigned>(seed));
                std::shuffle(indices.begin(), indices.end(), generator);
                std::uniform_int_distribution<int> column_distribution(0, dataset.dim - 1);
                for (int index = 0; index < outlier_count; ++index) {
                    const int row_index = indices[index];
                    const int column = column_distribution(generator);
                    if (column >= static_cast<int>(rows[row_index].size())) continue;
                    double range = maximum[column] - minimum[column];
                    if (!(range > 0.0)) range = std::max(1.0, std::abs(maximum[column]));
                    const double value = index % 2 == 0
                        ? maximum[column] + range * 2.0
                        : minimum[column] - range * 2.0;
                    std::ostringstream formatted;
                    formatted << std::setprecision(17) << value;
                    rows[row_index][column] = formatted.str();
                }

                const std::string input = cache + "/" + dataset.name + "_"
                    + std::to_string(static_cast<int>(ratio * 1000000.0)) + "_"
                    + std::to_string(run) + ".csv";
                write_string_rows(rows, input);
                ACluster compressor;
                int mutable_dim = dataset.dim;
                const auto encode_start = std::chrono::high_resolution_clock::now();
                auto encoded = compressor.encode_multidim_with_diag(
                    input, mutable_dim, 100, 10000, 10, 10
                );
                const auto encode_end = std::chrono::high_resolution_clock::now();
                const auto decode_start = std::chrono::high_resolution_clock::now();
                auto decoded = compressor.decode_multidim(encoded.first);
                const auto decode_end = std::chrono::high_resolution_clock::now();
                if (run == 0) verify_rows(read_csv_rows(input), decoded);
                size_sum += encoded.first.size();
                encode_sum += std::chrono::duration_cast<std::chrono::nanoseconds>(encode_end - encode_start).count();
                decode_sum += std::chrono::duration_cast<std::chrono::nanoseconds>(decode_end - decode_start).count();
                reference_sum += encoded.second.total_reference_points;
                std::filesystem::remove(input);
            }

            const double runs = static_cast<double>(context.runs);
            const double row_count = static_cast<double>(base_rows.size());
            const double average_size = size_sum / runs;
            csv << dataset.name << ',' << std::fixed << std::setprecision(3) << ratio << ','
                << context.runs << ',' << std::setprecision(2) << average_size << ','
                << std::setprecision(4) << (row_count * dataset.dim * sizeof(double) / average_size) << ','
                << std::setprecision(2) << (encode_sum / runs / row_count) << ','
                << (decode_sum / runs / row_count) << ',' << (reference_sum / runs) << '\n';
        }
    }
    std::filesystem::remove_all(cache);
}
