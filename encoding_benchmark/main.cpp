#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <numeric>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <fstream>      // For file I/O
#include <filesystem>   // For creating directories (C++17)
#include <algorithm>
#include <map>
#include <set>

#include "utils/csv_reader.h"
#include "algorithms/chimp.h"
#include "algorithms/elf.h"
#include "algorithms/gorilla.h"
#include "algorithms/rle.h"
#include "algorithms/huffman.h"
#include "algorithms/reger.h"
#include "algorithms/acluster.h"
#include "algorithms/kcluster.h"
#include "utils/type_converter.h"
#include "utils/debug.h" // Keep debug for now

void write_binary_file(const std::string& filename, const std::vector<uint8_t>& data) {
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile) {
        throw std::runtime_error("Could not open file for writing: " + filename);
    }
    outfile.write(reinterpret_cast<const char*>(data.data()), data.size());
}

std::vector<uint8_t> read_binary_file(const std::string& filename) {
    std::ifstream infile(filename, std::ios::binary | std::ios::ate);
    if (!infile) {
        throw std::runtime_error("Could not open file for reading: " + filename);
    }
    std::streamsize size = infile.tellg();
    infile.seekg(0, std::ios::beg);
    std::vector<uint8_t> buffer(size);
    infile.read(reinterpret_cast<char*>(buffer.data()), size);
    return buffer;
}


// For order-preserving algorithms (Chimp, Elf, Gorilla, RLE, Huffman)
void verify_data_ordered(const std::vector<double>& original, const std::vector<double>& decoded, const std::string& col_name) {
    if (original.size() != decoded.size()) {
        std::cerr << "Verification failed for " << col_name << ": original size (" << original.size() 
                  << ") != decoded size (" << decoded.size() << ")" << std::endl;
        throw std::runtime_error("Size mismatch during verification.");
    }
    constexpr double epsilon = 1e-6;
    for (size_t i = 0; i < original.size(); ++i) {
        bool error = std::isnan(original[i]) != std::isnan(decoded[i]) || 
                     (!std::isnan(original[i]) && std::abs(original[i] - decoded[i]) > epsilon);
        if (error) {
            std::cerr << "Verification failed for " << col_name << " at index " << i << ": original " 
                      << std::fixed << std::setprecision(8) << original[i] 
                      << " != decoded " << decoded[i] << std::endl;
            throw std::runtime_error("Value mismatch during verification.");
        }
    }
}

// For order-oblivious algorithms (Reger, KCluster, ACluster)
void verify_data_order_oblivious(std::vector<std::vector<double>>& original_rows, std::vector<std::vector<double>>& decoded_rows, const std::string& dataset_name) {
    if (original_rows.size() != decoded_rows.size()) {
        std::cerr << "Verification failed for " << dataset_name << ": original row count (" << original_rows.size() 
                  << ") != decoded row count (" << decoded_rows.size() << ")" << std::endl;
        throw std::runtime_error("Row count mismatch during verification.");
    }
    if (original_rows.empty()) return;

    // Sorting is a stable way to compare for equality, works for duplicates.
    std::sort(original_rows.begin(), original_rows.end());
    std::sort(decoded_rows.begin(), decoded_rows.end());
    
    constexpr double epsilon = 1e-9;
    for (size_t i = 0; i < original_rows.size(); ++i) {
        
        if (original_rows[i].size() != decoded_rows[i].size()) {
            std::cerr << "Verification FAILED for " << dataset_name << " after sorting." << std::endl;
            std::cerr << "Dimension mismatch at sorted row index " << i << ": original dim=" << original_rows[i].size()
                      << ", decoded dim=" << decoded_rows[i].size() << std::endl;
            throw std::runtime_error("Data mismatch: Row dimension mismatch.");
        }

        
        for (size_t j = 0; j < original_rows[i].size(); ++j) {
            
            bool is_nan_mismatch = std::isnan(original_rows[i][j]) != std::isnan(decoded_rows[i][j]);
            bool is_value_mismatch = !std::isnan(original_rows[i][j]) && std::abs(original_rows[i][j] - decoded_rows[i][j]) > epsilon;

            if (is_nan_mismatch || is_value_mismatch) {
                 std::cerr << "Verification FAILED for " << dataset_name << " after sorting." << std::endl;
                 std::cerr << "Mismatch found at sorted row index " << i << ", column " << j << "." << std::endl;
                 
                 std::cerr << std::fixed << std::setprecision(15);
                 std::cerr << "Original row [" << i << "]: {";
                 for(size_t k = 0; k < original_rows[i].size(); ++k) {
                     std::cerr << original_rows[i][k] << (k == original_rows[i].size() - 1 ? "" : ", ");
                 }
                 std::cerr << "}" << std::endl;

                 std::cerr << "Decoded row  [" << i << "]: {";
                 for(size_t k = 0; k < decoded_rows[i].size(); ++k) {
                     std::cerr << decoded_rows[i][k] << (k == decoded_rows[i].size() - 1 ? "" : ", ");
                 }
                 std::cerr << "}" << std::endl;
                 
                 std::cerr << "--> Mismatching value: Original=" << original_rows[i][j] << ", Decoded=" << decoded_rows[i][j] << std::endl;
                 
                 
                 if (i > 0) {
                     std::cerr << "Previous original row [" << i - 1 << "]: {";
                     for(size_t k = 0; k < original_rows[i-1].size(); ++k) std::cerr << original_rows[i-1][k] << (k == original_rows[i-1].size() - 1 ? "" : ", ");
                     std::cerr << "}" << std::endl;
                     std::cerr << "Previous decoded row  [" << i - 1 << "]: {";
                     for(size_t k = 0; k < decoded_rows[i-1].size(); ++k) std::cerr << decoded_rows[i-1][k] << (k == decoded_rows[i-1].size() - 1 ? "" : ", ");
                     std::cerr << "}" << std::endl;
                 }
                 if (i + 1 < original_rows.size()) {
                     std::cerr << "Next original row [" << i + 1 << "]: {";
                     for(size_t k = 0; k < original_rows[i+1].size(); ++k) std::cerr << original_rows[i+1][k] << (k == original_rows[i+1].size() - 1 ? "" : ", ");
                     std::cerr << "}" << std::endl;
                     std::cerr << "Next decoded row  [" << i + 1 << "]: {";
                     for(size_t k = 0; k < decoded_rows[i+1].size(); ++k) std::cerr << decoded_rows[i+1][k] << (k == decoded_rows[i+1].size() - 1 ? "" : ", ");
                     std::cerr << "}" << std::endl;
                 }


                 throw std::runtime_error("Data mismatch for order-oblivious algorithm.");
            }
        }
    }
}

// template<typename Compressor>
// void run_columnar_encoder_test(
//     const std::string& dataset_name, 
//     int expected_dim, 
//     int expected_num_rows, 
//     const std::string& csv_file_path,
//     const std::string& output_dir
// ) {
//     Compressor compressor;
//     std::string algorithm_name = compressor.get_name();

//     std::cout << "====================================================\n";
//     std::cout << "Dataset: " << dataset_name << ".csv, Method: " << algorithm_name << "\n";
//     std::cout << "----------------------------------------------------\n";

//     long long total_compressed_size_bytes = 0;
//     long long total_encoding_time_ns = 0;
//     long long total_decoding_time_ns = 0;
    

//     auto start_full_encode = std::chrono::high_resolution_clock::now();
//     std::vector<std::vector<double>> all_columns = read_csv_columns(csv_file_path);


//     std::string algo_output_dir = output_dir + "/" + algorithm_name;
//     if (!std::filesystem::exists(algo_output_dir)) {
//         std::filesystem::create_directory(algo_output_dir);
//     }

//     for (int i = 0; i < all_columns.size(); ++i) {
//         std::vector<uint8_t> compressed_data = compressor.encode(all_columns[i]);
//         total_compressed_size_bytes += compressed_data.size();
//         std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
//         write_binary_file(compressed_filename, compressed_data);
//     }
//     auto end_full_encode = std::chrono::high_resolution_clock::now();
//     total_encoding_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count();


//     auto start_full_decode = std::chrono::high_resolution_clock::now();
//     std::vector<std::vector<double>> decoded_columns;
//     for (int i = 0; i < all_columns.size(); ++i) {
//         std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
//         std::vector<uint8_t> compressed_data = read_binary_file(compressed_filename);
//         decoded_columns.push_back(compressor.decode(compressed_data));
//     }
//     auto end_full_decode = std::chrono::high_resolution_clock::now();
//     total_decoding_time_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count();
    
//      // FIX 1: Pass all required arguments to verify_data_ordered
//     for(size_t i = 0; i < all_columns.size(); ++i) {
//         verify_data_ordered(all_columns[i], decoded_columns[i], algorithm_name + " col_" + std::to_string(i));
//     }



//     double per_point_encoding_ns = static_cast<double>(total_encoding_time_ns) / expected_num_rows;
//     double per_point_decoding_ns = static_cast<double>(total_decoding_time_ns) / expected_num_rows;
//     long long total_points = (long long)expected_num_rows * expected_dim;
//     double bytes_per_point = static_cast<double>(total_compressed_size_bytes) / total_points;

//     std::cout << "Verification successful!" << std::endl;
//     std::cout << std::fixed << std::setprecision(2);
//     std::cout << "Total Byte Cost:       " << total_compressed_size_bytes << " bytes" << std::endl;
//     std::cout << "Bytes per point:       " << bytes_per_point << std::endl;
//     std::cout << "Encoding time (IO incl): " << per_point_encoding_ns << " ns per point" << std::endl;
//     std::cout << "Decoding time (IO incl): " << per_point_decoding_ns << " ns per point" << std::endl;
//     std::cout << "====================================================\n" << std::endl;
// }



std::vector<uint8_t> columns_to_bytes(const std::vector<std::vector<double>>& columns) {
    if (columns.empty()) return {};
    size_t total_doubles = columns.size() * columns[0].size();
    std::vector<uint8_t> byte_stream;
    byte_stream.reserve(total_doubles * sizeof(double));

    for (const auto& col : columns) {
        for (double val : col) {
            uint64_t long_val = double_to_long(val); // Assuming you have this helper from before
            for (size_t i = 0; i < sizeof(double); ++i) {
                byte_stream.push_back((long_val >> (i * 8)) & 0xFF);
            }
        }
    }
    return byte_stream;
}


std::vector<std::vector<double>> bytes_to_columns(const std::vector<uint8_t>& byte_stream, int num_cols) {
    if (byte_stream.empty() || num_cols == 0) return {};
    size_t total_bytes = byte_stream.size();
    if (total_bytes % sizeof(double) != 0) throw std::runtime_error("Invalid byte stream size for double conversion.");
    
    size_t total_doubles = total_bytes / sizeof(double);
    if (total_doubles % num_cols != 0) throw std::runtime_error("Cannot divide doubles evenly into columns.");
    
    size_t num_rows = total_doubles / num_cols;
    std::vector<std::vector<double>> columns(num_cols, std::vector<double>(num_rows));
    
    size_t byte_idx = 0;
    for (int c = 0; c < num_cols; ++c) {
        for (int r = 0; r < num_rows; ++r) {
            uint64_t long_val = 0;
            for (size_t i = 0; i < sizeof(double); ++i) {
                long_val |= (static_cast<uint64_t>(byte_stream[byte_idx++]) << (i * 8));
            }
            columns[c][r] = long_to_double(long_val); // Assuming you have this helper
        }
    }
    return columns;
}

template<typename Compressor>
void run_columnar_encoder_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir
) {
    Compressor compressor;
    std::string algorithm_name = compressor.get_name();

    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ".csv, Method: " << algorithm_name << "\n";
    std::cout << "----------------------------------------------------\n";

    const int TOTAL_RUNS = 15;
    const int AVG_RUNS = 10;

    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns;

    for (int run = 0; run < TOTAL_RUNS; ++run) {
        long long total_compressed_size_bytes = 0;
        
        // --- Encoding ---
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> all_columns = read_csv_columns(csv_file_path);
        if (all_columns.empty() || all_columns[0].size() != expected_num_rows) {
             std::cerr << "Warning: Data loading issue for " << dataset_name << ". Skipping test." << std::endl;
            return;
        }
        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        std::filesystem::create_directories(algo_output_dir);

        for (size_t i = 0; i < all_columns.size(); ++i) {
            std::vector<uint8_t> compressed_data = compressor.encode(all_columns[i]);
            total_compressed_size_bytes += compressed_data.size();
            std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
            write_binary_file(compressed_filename, compressed_data);
        }
        auto end_full_encode = std::chrono::high_resolution_clock::now();
        
        // --- Decoding ---
        auto start_full_decode = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<uint8_t>> all_compressed_data;
        all_compressed_data.reserve(all_columns.size());
        for (size_t i = 0; i < all_columns.size(); ++i) {
            std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
            all_compressed_data.push_back(read_binary_file(compressed_filename));
        }
        std::vector<std::vector<double>> decoded_columns;
        decoded_columns.reserve(all_columns.size());
        for (const auto& compressed_col_data : all_compressed_data) {
            decoded_columns.push_back(compressor.decode(compressed_col_data));
        }
        auto end_full_decode = std::chrono::high_resolution_clock::now();

        // --- Verification (run only once) ---
        if (run == 0) {
            for(size_t i = 0; i < all_columns.size(); ++i) {
                verify_data_ordered(all_columns[i], decoded_columns[i], algorithm_name + " col_" + std::to_string(i));
            }
            std::cout << "Verification successful!" << std::endl;
        }

        all_compressed_sizes.push_back(total_compressed_size_bytes);
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count());
    }

    // --- Calculate and Print Averages ---
    double avg_compressed_size = static_cast<double>(std::accumulate(all_compressed_sizes.end() - AVG_RUNS, all_compressed_sizes.end(), 0LL)) / AVG_RUNS;
    double avg_encoding_time_ns = static_cast<double>(std::accumulate(all_encoding_times_ns.end() - AVG_RUNS, all_encoding_times_ns.end(), 0LL)) / AVG_RUNS;
    double avg_decoding_time_ns = static_cast<double>(std::accumulate(all_decoding_times_ns.end() - AVG_RUNS, all_decoding_times_ns.end(), 0LL)) / AVG_RUNS;
    
    long long total_points = (long long)expected_num_rows * expected_dim;
    double avg_bytes_per_point = avg_compressed_size / total_points;
    double avg_per_point_encoding_ns = avg_encoding_time_ns / expected_num_rows;
    double avg_per_point_decoding_ns = avg_decoding_time_ns / expected_num_rows;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes (avg of last " << AVG_RUNS << " runs)" << std::endl;
    std::cout << "Bytes per point:       " << avg_bytes_per_point << std::endl;
    std::cout << "Encoding time (IO incl): " << avg_per_point_encoding_ns << " ns per point" << std::endl;
    std::cout << "Decoding time (IO incl): " << avg_per_point_decoding_ns << " ns per point" << std::endl;
    std::cout << "====================================================\n" << std::endl;
}


void run_multidim_reger_test(
    const std::string& dataset_name, int expected_dim, int expected_num_rows, 
    const std::string& csv_file_path, const std::string& output_dir)
{
    Reger compressor;
    std::string algorithm_name = compressor.get_name();
    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ", Method: " << algorithm_name << " (Multi-dim)\n";
    std::cout << "----------------------------------------------------\n";

    const int TOTAL_RUNS = 15;
    const int AVG_RUNS = 10;
    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns;

    for (int run = 0; run < TOTAL_RUNS; ++run) {
        // --- ENCODING ---
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<std::string>> str_rows = read_csv_rows_as_strings(csv_file_path);
        if (str_rows.empty()) { std::cout << "No data, skipping.\n"; return; }
        int num_rows = str_rows.size(); int num_cols = str_rows[0].size();
        std::vector<double> scaling_factors(num_cols);
        for (int j = 0; j < num_cols; ++j) {
            int max_d = 0;
            for (int i = 0; i < num_rows; ++i) {
                auto p = str_rows[i][j].find('.');
                if (p != std::string::npos) max_d = std::max(max_d, (int)(str_rows[i][j].length() - p - 1));
            }
            scaling_factors[j] = std::pow(10.0, max_d);
        }
        std::vector<std::vector<int64_t>> scaled_rows(num_rows, std::vector<int64_t>(num_cols));
        for (int i = 0; i < num_rows; ++i) {
            for (int j = 0; j < num_cols; ++j) {
                scaled_rows[i][j] = static_cast<int64_t>(std::round(std::stod(str_rows[i][j]) * scaling_factors[j]));
            }
        }
        std::vector<uint8_t> compressed_data = compressor.encode_multidim(scaled_rows, scaling_factors);
        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        std::filesystem::create_directories(algo_output_dir);
        std::string filename = algo_output_dir + "/" + dataset_name + ".bin";
        write_binary_file(filename, compressed_data);
        auto end_full_encode = std::chrono::high_resolution_clock::now();

        // --- DECODING ---
        auto start_full_decode = std::chrono::high_resolution_clock::now();
        auto compressed_read = read_binary_file(filename);
        auto [decoded_scaled_rows, decoded_factors] = compressor.decode_multidim(compressed_read);
        std::vector<std::vector<double>> decoded_rows;
        if (!decoded_scaled_rows.empty()) {
            decoded_rows.resize(decoded_scaled_rows.size(), std::vector<double>(decoded_scaled_rows[0].size()));
            for (size_t i = 0; i < decoded_scaled_rows.size(); ++i) {
                for (size_t j = 0; j < decoded_scaled_rows[i].size(); ++j) {
                    decoded_rows[i][j] = static_cast<double>(decoded_scaled_rows[i][j]) / decoded_factors[j];
                }
            }
        }
        auto end_full_decode = std::chrono::high_resolution_clock::now();

        // --- VERIFICATION (run only once) ---
        if (run == 0) {
            std::vector<std::vector<double>> original_rows = read_csv_rows(csv_file_path);
            verify_data_order_oblivious(original_rows, decoded_rows, dataset_name);
            std::cout << "Verification successful!\n";
        }

        all_compressed_sizes.push_back(compressed_data.size());
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count());
    }

    // --- Calculate and Print Averages ---
    double avg_compressed_size = static_cast<double>(std::accumulate(all_compressed_sizes.end() - AVG_RUNS, all_compressed_sizes.end(), 0LL)) / AVG_RUNS;
    double avg_encoding_time_ns = static_cast<double>(std::accumulate(all_encoding_times_ns.end() - AVG_RUNS, all_encoding_times_ns.end(), 0LL)) / AVG_RUNS;
    double avg_decoding_time_ns = static_cast<double>(std::accumulate(all_decoding_times_ns.end() - AVG_RUNS, all_decoding_times_ns.end(), 0LL)) / AVG_RUNS;
    
    long long total_points = (long long)expected_num_rows * expected_dim;
    double avg_bytes_per_point = avg_compressed_size / total_points;
    double avg_per_row_encoding_ns = avg_encoding_time_ns / expected_num_rows;
    double avg_per_row_decoding_ns = avg_decoding_time_ns / expected_num_rows;
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes (avg of last " << AVG_RUNS << " runs)\n";
    std::cout << "Bytes per point:       " << avg_bytes_per_point << "\n";
    std::cout << "Encoding time (IO incl): " << avg_per_row_encoding_ns << " ns per row\n";
    std::cout << "Decoding time (IO incl): " << avg_per_row_decoding_ns << " ns per row\n";
    std::cout << "====================================================\n\n";
}

template<typename ClusterCompressor>
void run_cluster_encoder_test(
    const std::string& dataset_name, int dim, int num_rows, 
    const std::string& csv_file_path, const std::string& output_dir)
{
    ClusterCompressor compressor;
    std::string algorithm_name = compressor.get_name();
    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ", Method: " << algorithm_name << " (Multi-dim)\n";
    std::cout << "----------------------------------------------------\n";

    const int TOTAL_RUNS = 15;
    const int AVG_RUNS = 10;
    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns;
    int k = 100, page_size = 10000, pack_size = 10, block_size = 10;
    GlobalDiagnostics final_diagnostics;

    for(int run = 0; run < TOTAL_RUNS; ++run) {
        
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        int mutable_dim = dim;
        
        // original version
        auto compressed_data = compressor.encode_multidim(csv_file_path, mutable_dim, k, page_size, pack_size, block_size);
        
        // diagnostic version 
        // auto [compressed_data, current_run_diagnostics] = compressor.encode_multidim_with_diag(
            // csv_file_path, mutable_dim, k, page_size, pack_size, block_size
        // );
        // final_diagnostics = current_run_diagnostics;

        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        std::filesystem::create_directories(algo_output_dir);
        std::string filename = algo_output_dir + "/" + dataset_name + ".bin";
        write_binary_file(filename, compressed_data);
        auto end_full_encode = std::chrono::high_resolution_clock::now();

        auto start_full_decode = std::chrono::high_resolution_clock::now();
        auto compressed_read = read_binary_file(filename);
        auto decoded_rows = compressor.decode_multidim(compressed_read);
        auto end_full_decode = std::chrono::high_resolution_clock::now();

        if (run == 0) {
            std::vector<std::vector<double>> original_rows = read_csv_rows(csv_file_path);
            verify_data_order_oblivious(original_rows, decoded_rows, dataset_name);
            std::cout << "Verification successful!\n";
        }

        all_compressed_sizes.push_back(compressed_data.size());
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count());
    }

    double avg_compressed_size = static_cast<double>(std::accumulate(all_compressed_sizes.end() - AVG_RUNS, all_compressed_sizes.end(), 0LL)) / AVG_RUNS;
    double avg_encoding_time_ns = static_cast<double>(std::accumulate(all_encoding_times_ns.end() - AVG_RUNS, all_encoding_times_ns.end(), 0LL)) / AVG_RUNS;
    double avg_decoding_time_ns = static_cast<double>(std::accumulate(all_decoding_times_ns.end() - AVG_RUNS, all_decoding_times_ns.end(), 0LL)) / AVG_RUNS;
    
    long long total_points = (long long)num_rows * dim;
    double avg_bytes_per_point = avg_compressed_size / total_points;
    double avg_per_row_encoding_ns = avg_encoding_time_ns / num_rows;
    double avg_per_row_decoding_ns = avg_decoding_time_ns / num_rows;

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes (avg of last " << AVG_RUNS << " runs)\n";
    std::cout << "Bytes per point:       " << avg_bytes_per_point << "\n";
    std::cout << "Encoding time (IO incl): " << avg_per_row_encoding_ns << " ns per row\n";
    std::cout << "Decoding time (IO incl): " << avg_per_row_decoding_ns << " ns per row\n";
    std::cout << "====================================================\n\n";

    size_t num_pages = (num_rows + page_size - 1) / page_size;
    // final_diagnostics.print_summary(num_pages);
    
    std::cout << "====================================================\n\n";
}


template<typename Compressor>
void run_byte_compression_test(const std::string& dataset, int dim, int num_rows, const std::string& path, const std::string& out_dir) {
    Compressor c; std::string name = c.get_name();
    std::cout << "====================================================\n" << "Dataset: " << dataset << ", Method: " << name << "\n" << "----------------------------------------------------\n";
    auto start_enc = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<double>> all_cols = read_csv_columns(path);
    auto byte_stream = columns_to_bytes(all_cols);
    auto compressed = c.encode(byte_stream);
    std::string algo_dir = out_dir + "/" + name; std::filesystem::create_directories(algo_dir);
    write_binary_file(algo_dir + "/" + dataset + ".bin", compressed);
    auto end_enc = std::chrono::high_resolution_clock::now();
    auto start_dec = std::chrono::high_resolution_clock::now();
    auto compressed_read = read_binary_file(algo_dir + "/" + dataset + ".bin");
    auto decoded_bytes = c.decode(compressed_read);
    auto decoded_cols = bytes_to_columns(decoded_bytes, dim);
    auto end_dec = std::chrono::high_resolution_clock::now();

    // FIX 2: Also fix the verify_data_ordered call here
    for(size_t i = 0; i < all_cols.size(); ++i) verify_data_ordered(all_cols[i], decoded_cols[i], name + " col_" + std::to_string(i));

    long long enc_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count();
    long long dec_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count();
    long long total_points = (long long)num_rows * dim;
    std::cout << "Verification successful!\n" << std::fixed << std::setprecision(2)
              << "Total Byte Cost:       " << compressed.size() << " bytes\n"
              << "Bytes per point:       " << (double)compressed.size() / total_points << "\n"
              << "Encoding time (IO incl): " << (double)enc_ns / num_rows << " ns per row\n"
              << "Decoding time (IO incl): " << (double)dec_ns / num_rows << " ns per row\n"
              << "====================================================\n\n";
}



template<typename ClusterCompressor>
void run_cluster_encoder_experiment_instance(
    const std::string& dataset_name, int dim, int num_rows, 
    const std::string& csv_file_path, const std::string& output_dir,
    int k, int page_size) // Parameters for the experiment
{
    ClusterCompressor compressor;
    std::string algorithm_name = compressor.get_name();
    
    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name 
              << ", Method: " << algorithm_name
              << ", k: " << k 
              << ", page_size: " << page_size << "\n";
    std::cout << "----------------------------------------------------\n";

    // These are fixed for simplicity in these experiments
    int pack_size = 10;
    int block_size = 10;
    
    auto start_full_encode = std::chrono::high_resolution_clock::now();
    int mutable_dim = dim;
    auto compressed_data = compressor.encode_multidim(csv_file_path, mutable_dim, k, page_size, pack_size, block_size);
    
    // Create a unique directory for this specific parameter combination to avoid file overwrites
    std::string param_string = algorithm_name + "_k" + std::to_string(k) + "_p" + std::to_string(page_size);
    std::string algo_output_dir = output_dir + "/" + dataset_name + "/" + param_string;
    std::filesystem::create_directories(algo_output_dir);
    std::string filename = algo_output_dir + "/" + dataset_name + ".bin";
    write_binary_file(filename, compressed_data);
    auto end_full_encode = std::chrono::high_resolution_clock::now();

    auto start_full_decode = std::chrono::high_resolution_clock::now();
    auto compressed_read = read_binary_file(filename);
    auto decoded_rows = compressor.decode_multidim(compressed_read);
    auto end_full_decode = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> original_rows = read_csv_rows(csv_file_path);
    verify_data_order_oblivious(original_rows, decoded_rows, dataset_name);

    long long enc_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count();
    long long dec_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count();
    long long total_points = (long long)num_rows * mutable_dim;
    
    std::cout << "Verification successful!\n" << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << compressed_data.size() << " bytes\n";
    std::cout << "Bytes per point:       " << (double)compressed_data.size() / total_points << "\n";
    std::cout << "Encoding time (IO incl): " << (double)enc_ns / num_rows << " ns per row\n";
    std::cout << "Decoding time (IO incl): " << (double)dec_ns / num_rows << " ns per row\n";
    std::cout << "====================================================\n\n";
}


void run_k_sensitivity_experiment(const std::string& file_root, const std::string& output_dir) {
    std::cout << "\n\n\n####################################################\n";
    std::cout << "###   STARTING K-VALUE SENSITIVITY EXPERIMENT  ###\n";
    std::cout << "####################################################\n\n\n";

    std::map<std::string, std::pair<int, int>> datasets_to_test = {
        {"Crop",    {46, 24000}},
        {"gas",     {64, 13910}}
    };

    std::vector<int> k_values = {10, 50, 100, 500, 1000};
    int fixed_page_size = 10000; // Page size is fixed for this experiment

    for (const auto& ds_pair : datasets_to_test) {
        const std::string& dataset_name = ds_pair.first;
        int dim = ds_pair.second.first;
        int num_rows = ds_pair.second.second;
        std::string file_path = file_root + dataset_name + ".csv";

        for (int k : k_values) {
            run_cluster_encoder_experiment_instance<KCluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                k, fixed_page_size
            );
        }
    }
    std::cout << "\n\n### K-VALUE SENSITIVITY EXPERIMENT COMPLETE ###\n\n";
}

void run_page_size_sensitivity_experiment(const std::string& file_root, const std::string& output_dir) {
    std::cout << "\n\n\n#######################################################\n";
    std::cout << "###   STARTING PAGE SIZE SENSITIVITY EXPERIMENT   ###\n";
    std::cout << "#######################################################\n\n\n";
    
    std::map<std::string, std::pair<int, int>> datasets_to_test = {
        {"Crop",    {46, 24000}},
        {"gas",     {64, 13910}}
    };

    std::vector<int> page_size_values = {1000, 2000, 5000, 10000, 20000};
    int fixed_k_for_kcluster = 100; // K is fixed for KCluster in this experiment

    for (const auto& ds_pair : datasets_to_test) {
        const std::string& dataset_name = ds_pair.first;
        int dim = ds_pair.second.first;
        int num_rows = ds_pair.second.second;
        std::string file_path = file_root + dataset_name + ".csv";

        for (int page_size : page_size_values) {
            // Run for KCluster
            run_cluster_encoder_experiment_instance<KCluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                fixed_k_for_kcluster, page_size
            );

            // Run for ACluster (k value is ignored, pass 0 as placeholder)
            run_cluster_encoder_experiment_instance<ACluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                0, page_size
            );
        }
    }
    std::cout << "\n\n### PAGE SIZE SENSITIVITY EXPERIMENT COMPLETE ###\n\n";
}


int main() {
    try {
        std::string fileRoot = "../data/";
        std::string outputDir = "output";
        std::filesystem::create_directories(outputDir);

        std::vector<std::string> datasets = {"SSD-bench", "profile-income", "Blockchain-tr", "Crop", "gas"};
        std::vector<int> data_size = {8926, 14825, 99999, 24000, 13910};
        
        std::vector<int> dims = {1, 1, 1, 46, 64};
        run_k_sensitivity_experiment(fileRoot, outputDir);
        run_page_size_sensitivity_experiment(fileRoot, outputDir);

        for (size_t i = 0; i < datasets.size(); ++i) {
            std::string file_path = fileRoot + datasets[i] + ".csv";
            
            run_columnar_encoder_test<Chimp>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            run_columnar_encoder_test<Elf>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            run_columnar_encoder_test<Gorilla>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            
            run_columnar_encoder_test<Rle>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            
            run_columnar_encoder_test<Huffman>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            
            run_multidim_reger_test(datasets[i], dims[i], data_size[i], file_path, outputDir);
            
            run_cluster_encoder_test<KCluster>(datasets[i], dims[i], data_size[i], file_path, outputDir);
            run_cluster_encoder_test<ACluster>(datasets[i], dims[i], data_size[i], file_path, outputDir);
        }
        
    } catch (const std::exception& e) {
        std::cerr << "A critical error occurred in main: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}