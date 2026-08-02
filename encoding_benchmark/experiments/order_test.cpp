#include "experiment_api.h"
#include "standalone_support.h"

#include "algorithms/acluster.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>

void run_order_test(const ExperimentContext& context) {
    const std::string output = experiment_output(context, "order");
    const std::string cache = output + "/cache";
    std::filesystem::create_directories(cache);
    std::ofstream csv(output + "/results.csv");
    csv << "Scenario,Dataset,Runs,CompressedBytes,CompressionRatio,EncodeNsPerRow,DecodeNsPerRow,References\n";

    for (const DatasetConfig& dataset : context.datasets) {
        const std::string source = dataset_path(context, dataset);
        const auto source_rows = read_csv_rows_as_strings(source);
        for (const std::string scenario : {"original", "random_shuffle", "sort_by_first_col"}) {
            double size_sum = 0.0;
            double encode_sum = 0.0;
            double decode_sum = 0.0;
            double reference_sum = 0.0;
            for (int run = 0; run < context.runs; ++run) {
                std::string input = source;
                if (scenario != "original") {
                    auto rows = source_rows;
                    if (scenario == "random_shuffle") {
                        std::mt19937 generator(20260406U + static_cast<unsigned>(run));
                        std::shuffle(rows.begin(), rows.end(), generator);
                    } else {
                        std::stable_sort(rows.begin(), rows.end(), [](const auto& left, const auto& right) {
                            const double left_value = left.empty() ? 0.0 : parse_number(left.front());
                            const double right_value = right.empty() ? 0.0 : parse_number(right.front());
                            return left_value < right_value;
                        });
                    }
                    input = cache + "/" + dataset.name + "_" + scenario + "_" + std::to_string(run) + ".csv";
                    write_string_rows(rows, input);
                }

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
                if (input != source) std::filesystem::remove(input);
            }
            const double runs = static_cast<double>(context.runs);
            const double rows = static_cast<double>(source_rows.size());
            const double average_size = size_sum / runs;
            csv << scenario << ',' << dataset.name << ',' << context.runs << ','
                << std::fixed << std::setprecision(2) << average_size << ','
                << std::setprecision(4) << (rows * dataset.dim * sizeof(double) / average_size) << ','
                << std::setprecision(2) << (encode_sum / runs / rows) << ','
                << (decode_sum / runs / rows) << ',' << (reference_sum / runs) << '\n';
        }
    }
    std::filesystem::remove_all(cache);
}
