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
#include <functional>
#include <thread>
#include <atomic>
#include <random>
#include <limits>
#include <exception>
#include <type_traits>
#include <windows.h>
#include <psapi.h>

#include "utils/csv_reader.h"
#include "algorithms/chimp.h"
#include "algorithms/elf.h"
#include "algorithms/gorilla.h"
#include "algorithms/rle.h"
#include "algorithms/huffman.h"
#include "algorithms/reger.h"
#include "algorithms/acluster.h"
#include "algorithms/kcluster.h"
#include "algorithms/kcluster_kmeans.h"
#include "algorithms/acluster_sub.h"
#include "algorithms/kcluster_sub.h"
#include "utils/type_converter.h"
#include "algorithms/alp.h"
#include "algorithms/secondary_compressor.h"
#include "algorithms/vortex.h"
#include "utils/debug.h" 

// ==========================================
// Global CSV Output Handler & Stats Utils
// ==========================================
std::ofstream resultsFile;
std::ofstream memoryResultsFile;
std::map<std::string, std::vector<double>> methodMemoryBytes;
const int PAGE_POINTS = 10000;
std::string currentScenario = "original";

struct ScenarioResultRecord {
    std::string scenario;
    std::string dataset;
    std::string algorithm;
    double compression_ratio;
    double enc_avg_ns;
};

std::vector<ScenarioResultRecord> scenarioResults;

struct LoggedResultRecord {
    std::string scenario;
    std::string dataset;
    std::string algorithm;
    int dim;
    int rows;
    double compression_ratio;
    double enc_avg_ns;
    double dec_avg_ns;
};

std::vector<LoggedResultRecord> loggedResults;

void init_results_csv(const std::string& filepath) {
    resultsFile.open(filepath);
    if (!resultsFile.is_open()) {
        std::cerr << "Error: Could not open results file for writing: " << filepath << std::endl;
        exit(1);
    }
    // Write Header
    resultsFile << "Scenario,Dataset,Algorithm,Dimensions,Rows,OriginalSize_Bytes,CompressedSize_Bytes,CompressionRatio,BytesPerPoint,"
                << "EncTime_Avg_ns,EncTime_Min_ns,EncTime_Max_ns,"
                << "DecTime_Avg_ns,DecTime_Min_ns,DecTime_Max_ns\n";
}

void init_memory_csv(const std::string& filepath) {
    memoryResultsFile.open(filepath);
    if (!memoryResultsFile.is_open()) {
        std::cerr << "Error: Could not open memory results file for writing: " << filepath << std::endl;
        exit(1);
    }
    memoryResultsFile << "Scenario,Dataset,Algorithm,PageRows,PeakEncodePrivateBytes\n";
}

void log_memory_to_csv(const std::string& dataset, const std::string& algo, int page_rows, double peak_bytes) {
    if (!memoryResultsFile.is_open()) return;
    memoryResultsFile << currentScenario << "," << dataset << "," << algo << "," << page_rows << "," << std::fixed << std::setprecision(0) << peak_bytes << "\n";
    memoryResultsFile.flush();
    methodMemoryBytes[currentScenario + "|" + algo].push_back(peak_bytes);
}

void write_memory_method_average_csv(const std::string& filepath) {
    std::ofstream out(filepath);
    if (!out.is_open()) return;
    out << "Algorithm,AvgPeakEncodePrivateBytes\n";
    for (const auto& kv : methodMemoryBytes) {
        if (kv.second.empty()) continue;
        double avg = std::accumulate(kv.second.begin(), kv.second.end(), 0.0) / kv.second.size();
        out << kv.first << "," << std::fixed << std::setprecision(0) << avg << "\n";
    }
}

void print_memory_method_average() {
    std::cout << "=============== Memory Average (Per Method) ===============\n";
    for (const auto& kv : methodMemoryBytes) {
        if (kv.second.empty()) continue;
        double avg = std::accumulate(kv.second.begin(), kv.second.end(), 0.0) / kv.second.size();
        std::cout << kv.first << ": " << std::fixed << std::setprecision(0) << avg << " bytes\n";
    }
    std::cout << "============================================================\n";
}

void log_result_to_csv(
    const std::string& dataset, const std::string& algo, int dim, int rows, 
    long long orig_size, double comp_size, 
    double enc_avg, double enc_min, double enc_max,
    double dec_avg, double dec_min, double dec_max
) {
    if (!resultsFile.is_open()) return;
    
    double ratio = (comp_size > 0) ? (double)orig_size / comp_size : 0;
    double bpp = (rows * dim > 0) ? comp_size / (double)(rows * dim) : 0;

    resultsFile << currentScenario << "," << dataset << "," << algo << "," << dim << "," << rows << ","
                << orig_size << "," << std::fixed << std::setprecision(0) << comp_size << ","
                << std::setprecision(4) << ratio << "," << bpp << ","
                << std::setprecision(2) << enc_avg << "," << enc_min << "," << enc_max << ","
                << dec_avg << "," << dec_min << "," << dec_max << "\n";
    resultsFile.flush(); 

    scenarioResults.push_back({currentScenario, dataset, algo, ratio, enc_avg});
    loggedResults.push_back({currentScenario, dataset, algo, dim, rows, ratio, enc_avg, dec_avg});
}

uint64_t get_process_private_bytes() {
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc), sizeof(pmc))) {
        return static_cast<uint64_t>(pmc.PrivateUsage);
    }
    return 0;
}

uint64_t measure_peak_private_delta_bytes(const std::function<void()>& fn) {
    uint64_t baseline = get_process_private_bytes();
    std::atomic<bool> running(true);
    std::atomic<uint64_t> peak(baseline);
    std::thread sampler([&]() {
        while (running.load(std::memory_order_relaxed)) {
            uint64_t current = get_process_private_bytes();
            uint64_t old = peak.load(std::memory_order_relaxed);
            while (current > old && !peak.compare_exchange_weak(old, current, std::memory_order_relaxed)) {}
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    });
    std::exception_ptr eptr;
    try {
        fn();
    } catch (...) {
        eptr = std::current_exception();
    }
    running.store(false, std::memory_order_relaxed);
    if (sampler.joinable()) sampler.join();
    uint64_t final_now = get_process_private_bytes();
    uint64_t old = peak.load(std::memory_order_relaxed);
    while (final_now > old && !peak.compare_exchange_weak(old, final_now, std::memory_order_relaxed)) {}
    uint64_t peak_value = peak.load(std::memory_order_relaxed);
    if (eptr) std::rethrow_exception(eptr);
    return (peak_value > baseline) ? (peak_value - baseline) : 0;
}

struct TimeStats {
    double avg;
    double min;
    double max;
};

// Calculate stats per ROW/POINT based on divisor
TimeStats calculate_stats(const std::vector<long long>& times_ns, int divisor) {
    if (times_ns.empty() || divisor == 0) return {0.0, 0.0, 0.0};
    
    auto [min_it, max_it] = std::minmax_element(times_ns.begin(), times_ns.end());
    double sum = std::accumulate(times_ns.begin(), times_ns.end(), 0.0);
    
    return {
        (sum / times_ns.size()) / divisor, 
        (double)(*min_it) / divisor,       
        (double)(*max_it) / divisor        
    };
}

// ==========================================
// File I/O Helpers
// ==========================================

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

void cleanup_method_output(const std::string& output_dir, const std::string& algorithm_name) {
    std::filesystem::remove_all(output_dir + "/" + algorithm_name);
}

std::vector<std::vector<double>> slice_columns_to_page(std::vector<std::vector<double>> columns, int& page_rows) {
    if (columns.empty()) {
        page_rows = 0;
        return columns;
    }
    page_rows = std::min(PAGE_POINTS, static_cast<int>(columns[0].size()));
    for (auto& col : columns) col.resize(page_rows);
    return columns;
}

std::vector<std::vector<double>> slice_rows_to_page(std::vector<std::vector<double>> rows, int& page_rows) {
    page_rows = std::min(PAGE_POINTS, static_cast<int>(rows.size()));
    rows.resize(page_rows);
    return rows;
}

std::vector<std::vector<std::string>> slice_string_rows_to_page(std::vector<std::vector<std::string>> rows, int& page_rows) {
    page_rows = std::min(PAGE_POINTS, static_cast<int>(rows.size()));
    rows.resize(page_rows);
    return rows;
}

std::string write_page_csv_from_string_rows(const std::vector<std::vector<std::string>>& rows, const std::string& out_path) {
    std::ofstream out(out_path);
    if (!out.is_open()) throw std::runtime_error("Could not write page csv: " + out_path);
    for (const auto& row : rows) {
        for (size_t i = 0; i < row.size(); ++i) {
            if (i) out << ",";
            out << row[i];
        }
        out << "\n";
    }
    return out_path;
}

enum class RowOrderMode {
    Original,
    RandomShuffle,
    SortByFirstColumn
};

std::string scenario_name(RowOrderMode mode) {
    if (mode == RowOrderMode::Original) return "original";
    if (mode == RowOrderMode::RandomShuffle) return "random_shuffle";
    return "sort_by_first_col";
}

double to_double_or_zero(const std::string& s) {
    try {
        return std::stod(s);
    } catch (...) {
        return 0.0;
    }
}

std::string prepare_dataset_for_mode(
    const std::string& source_csv_path,
    const std::string& dataset_name,
    RowOrderMode mode,
    const std::string& cache_root
) {
    if (mode == RowOrderMode::Original) return source_csv_path;
    std::vector<std::vector<std::string>> rows = read_csv_rows_as_strings(source_csv_path);
    if (rows.empty()) return source_csv_path;

    if (mode == RowOrderMode::RandomShuffle) {
        std::mt19937 rng(20260406);
        std::shuffle(rows.begin(), rows.end(), rng);
    } else {
        std::stable_sort(rows.begin(), rows.end(), [](const auto& a, const auto& b) {
            double av = a.empty() ? 0.0 : to_double_or_zero(a[0]);
            double bv = b.empty() ? 0.0 : to_double_or_zero(b[0]);
            return av < bv;
        });
    }

    std::filesystem::create_directories(cache_root);
    std::string out_path = cache_root + "/" + dataset_name + "__" + scenario_name(mode) + ".csv";
    return write_page_csv_from_string_rows(rows, out_path);
}

void write_scenario_summaries(const std::string& output_dir) {
    std::ofstream out_acluster(output_dir + "/summary_acluster_by_dataset_scenario.csv");
    out_acluster << "Scenario,Dataset,CompressionRatio,EncTime_Avg_ns\n";
    for (const auto& r : scenarioResults) {
        if (r.algorithm == "ACluster") {
            out_acluster << r.scenario << "," << r.dataset << "," << std::fixed << std::setprecision(4)
                         << r.compression_ratio << "," << std::setprecision(2) << r.enc_avg_ns << "\n";
        }
    }

    std::map<std::string, std::vector<std::pair<double, double>>> agg;
    for (const auto& r : scenarioResults) {
        agg[r.scenario + "|" + r.algorithm].push_back({r.compression_ratio, r.enc_avg_ns});
    }

    std::ofstream out_method(output_dir + "/summary_method_avg_by_scenario.csv");
    out_method << "Scenario,Algorithm,AvgCompressionRatio,AvgEncTime_Avg_ns\n";
    for (const auto& kv : agg) {
        const auto& vec = kv.second;
        if (vec.empty()) continue;
        double ratio_sum = 0.0, enc_sum = 0.0;
        for (const auto& p : vec) {
            ratio_sum += p.first;
            enc_sum += p.second;
        }
        size_t split = kv.first.find('|');
        std::string scenario = kv.first.substr(0, split);
        std::string algo = kv.first.substr(split + 1);
        out_method << scenario << "," << algo << "," << std::fixed << std::setprecision(4)
                   << (ratio_sum / vec.size()) << "," << std::setprecision(2) << (enc_sum / vec.size()) << "\n";
    }
}

// ==========================================
// Verification Helpers
// ==========================================

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

void verify_data_order_oblivious(std::vector<std::vector<double>>& original_rows, std::vector<std::vector<double>>& decoded_rows, const std::string& dataset_name) {
    if (original_rows.size() != decoded_rows.size()) {
        std::cerr << "Verification failed for " << dataset_name << ": original row count (" << original_rows.size() 
                  << ") != decoded row count (" << decoded_rows.size() << ")" << std::endl;
        throw std::runtime_error("Row count mismatch during verification.");
    }
    if (original_rows.empty()) return;

    // Sorting is a stable way to compare for equality, works for duplicates.
    // Create copies to sort so we don't mess up subsequent runs if we reused data (though we reload usually)
    auto sorted_orig = original_rows;
    auto sorted_dec = decoded_rows;
    
    std::sort(sorted_orig.begin(), sorted_orig.end());
    std::sort(sorted_dec.begin(), sorted_dec.end());
    
    constexpr double epsilon = 1e-6;
    for (size_t i = 0; i < sorted_orig.size(); ++i) {
        if (sorted_orig[i].size() != sorted_dec[i].size()) {
            std::cerr << "Dimension mismatch at sorted row index " << i << std::endl;
            throw std::runtime_error("Data mismatch: Row dimension mismatch.");
        }

        for (size_t j = 0; j < sorted_orig[i].size(); ++j) {
            bool is_nan_mismatch = std::isnan(sorted_orig[i][j]) != std::isnan(sorted_dec[i][j]);
            bool is_value_mismatch = !std::isnan(sorted_orig[i][j]) && std::abs(sorted_orig[i][j] - sorted_dec[i][j]) > epsilon;

            if (is_nan_mismatch || is_value_mismatch) {
                 std::cerr << "Verification FAILED for " << dataset_name << " after sorting." << std::endl;
                 std::cerr << "Mismatch found at sorted row index " << i << ", column " << j << "." << std::endl;
                 std::cerr << "Original=" << sorted_orig[i][j] << ", Decoded=" << sorted_dec[i][j] << std::endl;
                 throw std::runtime_error("Data mismatch for order-oblivious algorithm.");
            }
        }
    }
}

std::vector<uint8_t> columns_to_bytes(const std::vector<std::vector<double>>& columns) {
    if (columns.empty()) return {};
    size_t total_doubles = columns.size() * columns[0].size();
    std::vector<uint8_t> byte_stream;
    byte_stream.reserve(total_doubles * sizeof(double));

    for (const auto& col : columns) {
        for (double val : col) {
            uint64_t long_val = double_to_long(val); 
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
            columns[c][r] = long_to_double(long_val); 
        }
    }
    return columns;
}

// ==========================================
// Test Runner Functions
// ==========================================

const int NUM_RUNS = 1; 

std::vector<uint8_t> rows_to_bytes_row_major(const std::vector<std::vector<double>>& rows) {
    if (rows.empty()) return {};
    size_t dim = rows[0].size();
    std::vector<uint8_t> bytes;
    bytes.reserve(rows.size() * dim * sizeof(double));
    for (const auto& row : rows) {
        if (row.size() != dim) {
            throw std::runtime_error("Inconsistent row width while serializing rows.");
        }
        for (double v : row) {
            uint64_t bits = double_to_long(v);
            for (size_t b = 0; b < sizeof(double); ++b) {
                bytes.push_back(static_cast<uint8_t>((bits >> (b * 8)) & 0xFF));
            }
        }
    }
    return bytes;
}

std::vector<std::vector<double>> bytes_to_rows_row_major(const std::vector<uint8_t>& bytes, int rows, int dim) {
    if (rows <= 0 || dim <= 0) return {};
    size_t expected = static_cast<size_t>(rows) * static_cast<size_t>(dim) * sizeof(double);
    if (bytes.size() != expected) {
        throw std::runtime_error("Invalid payload size for row-major restoration.");
    }
    std::vector<std::vector<double>> restored(rows, std::vector<double>(dim));
    size_t idx = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < dim; ++j) {
            uint64_t bits = 0;
            for (size_t b = 0; b < sizeof(double); ++b) {
                bits |= (static_cast<uint64_t>(bytes[idx++]) << (b * 8));
            }
            restored[i][j] = long_to_double(bits);
        }
    }
    return restored;
}

void append_u32(std::vector<uint8_t>& out, uint32_t v) {
    out.push_back(static_cast<uint8_t>(v & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
}

uint32_t read_u32(const std::vector<uint8_t>& in, size_t& pos) {
    if (pos + 4 > in.size()) throw std::runtime_error("Corrupted Vortex blob: unexpected EOF.");
    uint32_t v = static_cast<uint32_t>(in[pos]) |
                 (static_cast<uint32_t>(in[pos + 1]) << 8) |
                 (static_cast<uint32_t>(in[pos + 2]) << 16) |
                 (static_cast<uint32_t>(in[pos + 3]) << 24);
    pos += 4;
    return v;
}

uint64_t read_u64(const std::vector<uint8_t>& in, size_t& pos) {
    uint64_t lo = static_cast<uint64_t>(read_u32(in, pos));
    uint64_t hi = static_cast<uint64_t>(read_u32(in, pos));
    return lo | (hi << 32);
}

bool parse_sub_cost_metrics(
    const std::vector<uint8_t>& payload,
    double& header_bits_per_row,
    double& base_bits_per_row,
    double& residual_bits_per_row,
    double& id_extra_bits_per_row
) {
    header_bits_per_row = 0.0;
    base_bits_per_row = 0.0;
    residual_bits_per_row = 0.0;
    id_extra_bits_per_row = 0.0;
    if (payload.size() < 6) return false;
    bool is_sub = (payload[0] == 'K' || payload[0] == 'A') &&
                  payload[1] == 'S' && payload[2] == 'U' && payload[3] == 'B';
    if (!is_sub) return false;
    size_t pos = 5;
    uint8_t mode = payload[4];
    if (mode != 1) return false;
    if (pos >= payload.size()) return false;
    uint8_t strategy = payload[pos++];
    uint64_t header_bits = 0, base_bits = 0, residual_bits = 0, id_bits = 0;
    uint32_t row_count = 0;
    if (strategy == 0) {
        if (pos + 3 * 4 + 3 * 8 > payload.size()) return false;
        read_u32(payload, pos);
        row_count = read_u32(payload, pos);
        uint32_t group_count = read_u32(payload, pos);
        header_bits = read_u64(payload, pos);
        base_bits = read_u64(payload, pos);
        residual_bits = read_u64(payload, pos);
        if (row_count == 0 || group_count == 0) return false;
        if (pos + static_cast<size_t>(group_count) * 4 > payload.size()) return false;
        for (uint32_t i = 0; i < group_count; ++i) read_u32(payload, pos);
        if (pos + static_cast<size_t>(group_count) * 4 > payload.size()) return false;
        for (uint32_t i = 0; i < group_count; ++i) read_u32(payload, pos);
        if (pos + static_cast<size_t>(group_count - 1) * 4 > payload.size()) return false;
        for (uint32_t i = 0; i + 1 < group_count; ++i) id_bits += static_cast<uint64_t>(read_u32(payload, pos)) * 8ULL;
    } else if (strategy == 1) {
        if (pos + 6 * 4 + 4 * 8 > payload.size()) return false;
        read_u32(payload, pos);
        row_count = read_u32(payload, pos);
        read_u32(payload, pos);
        read_u32(payload, pos);
        read_u32(payload, pos);
        read_u32(payload, pos);
        header_bits = read_u64(payload, pos);
        base_bits = read_u64(payload, pos);
        residual_bits = read_u64(payload, pos);
        id_bits = read_u64(payload, pos);
        if (row_count == 0) return false;
    } else {
        return false;
    }
    header_bits_per_row = static_cast<double>(header_bits) / static_cast<double>(row_count);
    base_bits_per_row = static_cast<double>(base_bits) / static_cast<double>(row_count);
    residual_bits_per_row = static_cast<double>(residual_bits) / static_cast<double>(row_count);
    id_extra_bits_per_row = static_cast<double>(id_bits) / static_cast<double>(row_count);
    return true;
}

void verify_rows_ordered(const std::vector<std::vector<double>>& original, const std::vector<std::vector<double>>& decoded, const std::string& dataset) {
    if (original.size() != decoded.size()) {
        throw std::runtime_error("Verification failed for " + dataset + ": row count mismatch.");
    }
    constexpr double eps = 1e-6;
    for (size_t i = 0; i < original.size(); ++i) {
        if (original[i].size() != decoded[i].size()) {
            throw std::runtime_error("Verification failed for " + dataset + ": dim mismatch.");
        }
        for (size_t j = 0; j < original[i].size(); ++j) {
            bool nan_mismatch = std::isnan(original[i][j]) != std::isnan(decoded[i][j]);
            bool val_mismatch = !std::isnan(original[i][j]) && std::abs(original[i][j] - decoded[i][j]) > eps;
            if (nan_mismatch || val_mismatch) {
                throw std::runtime_error("Verification failed for " + dataset + ": value mismatch.");
            }
        }
    }
}

void run_vortex_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir
) {
    Vortex vortex;
    Rle rle;
    const std::string algorithm_name = vortex.get_name() + "+Rle";

    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ", Method: " << algorithm_name << "\n";
    std::cout << "----------------------------------------------------\n";

    std::vector<std::vector<double>> original_rows = read_csv_rows(csv_file_path);
    if (original_rows.empty()) {
        std::cerr << "Warning: Data loading issue for " << dataset_name << ". Skipping test." << std::endl;
        return;
    }
    int effective_num_rows = 0;
    original_rows = slice_rows_to_page(std::move(original_rows), effective_num_rows);
    if (static_cast<int>(original_rows[0].size()) != expected_dim) {
        std::cerr << "Warning: Dimension mismatch for " << dataset_name << ". Skipping test." << std::endl;
        return;
    }

    long long original_size_bytes = static_cast<long long>(effective_num_rows) * expected_dim * 8LL;
    std::vector<long long> all_compressed_sizes;
    std::vector<long long> all_encoding_times_ns;
    std::vector<long long> all_decoding_times_ns;
    std::vector<long long> all_peak_mem_bytes;

    for (int run = 0; run < NUM_RUNS; ++run) {
        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        if (run == 0) std::filesystem::create_directories(algo_output_dir);
        std::string prefix = algo_output_dir + "/" + dataset_name;

        std::vector<uint32_t> permutation;
        std::vector<std::vector<double>> reordered_rows;
        std::vector<std::vector<double>> reordered_columns(expected_dim, std::vector<double>(effective_num_rows));
        long long total_compressed_size_bytes = 0;
        auto start_enc = std::chrono::high_resolution_clock::now();
        uint64_t peak_mem_delta = measure_peak_private_delta_bytes([&]() {
            permutation = vortex.compute_permutation(original_rows);
            reordered_rows.reserve(original_rows.size());
            for (uint32_t idx : permutation) reordered_rows.push_back(original_rows[idx]);
            for (int r = 0; r < effective_num_rows; ++r) {
                for (int c = 0; c < expected_dim; ++c) {
                    reordered_columns[c][r] = reordered_rows[r][c];
                }
            }
            for (int c = 0; c < expected_dim; ++c) {
                std::vector<uint8_t> encoded = rle.encode(reordered_columns[c]);
                total_compressed_size_bytes += static_cast<long long>(encoded.size());
                std::string compressed_filename = prefix + "_col_" + std::to_string(c) + ".bin";
                write_binary_file(compressed_filename, encoded);
            }
            std::vector<double> permutation_as_double(effective_num_rows);
            for (int r = 0; r < effective_num_rows; ++r) {
                permutation_as_double[r] = static_cast<double>(permutation[r]);
            }
            std::vector<uint8_t> encoded_perm = rle.encode(permutation_as_double);
            total_compressed_size_bytes += static_cast<long long>(encoded_perm.size());
            std::string perm_filename = prefix + "_perm.bin";
            write_binary_file(perm_filename, encoded_perm);

            std::vector<uint8_t> meta(8);
            uint32_t rows_u32 = static_cast<uint32_t>(effective_num_rows);
            uint32_t dim_u32 = static_cast<uint32_t>(expected_dim);
            meta[0] = static_cast<uint8_t>(rows_u32 & 0xFF);
            meta[1] = static_cast<uint8_t>((rows_u32 >> 8) & 0xFF);
            meta[2] = static_cast<uint8_t>((rows_u32 >> 16) & 0xFF);
            meta[3] = static_cast<uint8_t>((rows_u32 >> 24) & 0xFF);
            meta[4] = static_cast<uint8_t>(dim_u32 & 0xFF);
            meta[5] = static_cast<uint8_t>((dim_u32 >> 8) & 0xFF);
            meta[6] = static_cast<uint8_t>((dim_u32 >> 16) & 0xFF);
            meta[7] = static_cast<uint8_t>((dim_u32 >> 24) & 0xFF);
            total_compressed_size_bytes += static_cast<long long>(meta.size());
            std::string meta_filename = prefix + "_meta.bin";
            write_binary_file(meta_filename, meta);
        });
        std::string perm_filename = prefix + "_perm.bin";
        std::string meta_filename = prefix + "_meta.bin";
        auto end_enc = std::chrono::high_resolution_clock::now();

        auto start_dec = std::chrono::high_resolution_clock::now();
        std::vector<uint8_t> meta_in = read_binary_file(meta_filename);
        if (meta_in.size() != 8) {
            throw std::runtime_error("Invalid Vortex meta file.");
        }
        uint32_t rows_read = static_cast<uint32_t>(meta_in[0]) |
                             (static_cast<uint32_t>(meta_in[1]) << 8) |
                             (static_cast<uint32_t>(meta_in[2]) << 16) |
                             (static_cast<uint32_t>(meta_in[3]) << 24);
        uint32_t dim_read = static_cast<uint32_t>(meta_in[4]) |
                            (static_cast<uint32_t>(meta_in[5]) << 8) |
                            (static_cast<uint32_t>(meta_in[6]) << 16) |
                            (static_cast<uint32_t>(meta_in[7]) << 24);
        if (rows_read != static_cast<uint32_t>(effective_num_rows) || dim_read != static_cast<uint32_t>(expected_dim)) {
            throw std::runtime_error("Vortex meta mismatch during decode.");
        }

        std::vector<std::vector<double>> decoded_columns;
        decoded_columns.reserve(expected_dim);
        for (int c = 0; c < expected_dim; ++c) {
            std::string compressed_filename = prefix + "_col_" + std::to_string(c) + ".bin";
            auto file_data = read_binary_file(compressed_filename);
            decoded_columns.push_back(rle.decode(file_data));
        }

        auto perm_data = read_binary_file(perm_filename);
        std::vector<double> decoded_perm_as_double = rle.decode(perm_data);
        if (decoded_perm_as_double.size() != static_cast<size_t>(effective_num_rows)) {
            throw std::runtime_error("Decoded permutation size mismatch.");
        }
        std::vector<uint32_t> decoded_perm(effective_num_rows);
        for (int r = 0; r < effective_num_rows; ++r) {
            double pv = decoded_perm_as_double[r];
            uint32_t pi = static_cast<uint32_t>(std::llround(pv));
            if (pi >= static_cast<uint32_t>(effective_num_rows) || std::abs(pv - static_cast<double>(pi)) > 1e-9) {
                throw std::runtime_error("Decoded permutation value is invalid.");
            }
            decoded_perm[r] = pi;
        }

        std::vector<std::vector<double>> decoded_reordered_rows(effective_num_rows, std::vector<double>(expected_dim));
        for (int c = 0; c < expected_dim; ++c) {
            if (decoded_columns[c].size() != static_cast<size_t>(effective_num_rows)) {
                throw std::runtime_error("Decoded column size mismatch.");
            }
            for (int r = 0; r < effective_num_rows; ++r) {
                decoded_reordered_rows[r][c] = decoded_columns[c][r];
            }
        }

        std::vector<std::vector<double>> restored_rows(effective_num_rows, std::vector<double>(expected_dim));
        for (int r = 0; r < effective_num_rows; ++r) {
            restored_rows[decoded_perm[r]] = decoded_reordered_rows[r];
        }
        auto end_dec = std::chrono::high_resolution_clock::now();

        if (run == 0) {
            verify_rows_ordered(original_rows, restored_rows, dataset_name);
            std::cout << "Verification successful!" << std::endl;
        }

        all_compressed_sizes.push_back(total_compressed_size_bytes);
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count());
        all_peak_mem_bytes.push_back(static_cast<long long>(peak_mem_delta));
    }

    double avg_compressed_size = std::accumulate(all_compressed_sizes.begin(), all_compressed_sizes.end(), 0.0) / NUM_RUNS;
    TimeStats enc_stats = calculate_stats(all_encoding_times_ns, effective_num_rows);
    TimeStats dec_stats = calculate_stats(all_decoding_times_ns, effective_num_rows);
    double compression_ratio = avg_compressed_size > 0 ? static_cast<double>(original_size_bytes) / avg_compressed_size : 0.0;
    double avg_peak_mem = std::accumulate(all_peak_mem_bytes.begin(), all_peak_mem_bytes.end(), 0.0) / NUM_RUNS;

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "Compression Ratio:     " << compression_ratio << "\n";
    std::cout << std::setprecision(2);
    std::cout << "Encoding (ns/row):     " << enc_stats.avg << "\n";
    std::cout << "Decoding (ns/row):     " << dec_stats.avg << "\n";
    std::cout << "Peak Encode Memory:    " << avg_peak_mem << " bytes\n";
    std::cout << "====================================================\n" << std::endl;

    log_result_to_csv(dataset_name, algorithm_name, expected_dim, effective_num_rows, original_size_bytes, avg_compressed_size,
                      enc_stats.avg, enc_stats.min, enc_stats.max,
                      dec_stats.avg, dec_stats.min, dec_stats.max);
    log_memory_to_csv(dataset_name, algorithm_name, effective_num_rows, avg_peak_mem);
    cleanup_method_output(output_dir, algorithm_name);
}

template<typename Compressor>
void run_columnar_encoder_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir,
    SecondaryCompressor* sec_comp = nullptr
) {
    // Handle secondary compressor defaults
    NoneCompressor default_none;
    if (!sec_comp) sec_comp = &default_none;

    Compressor compressor;
    std::string algorithm_name = compressor.get_name() + sec_comp->get_name(); // Combine names

    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ".csv, Method: " << algorithm_name << "\n";
    std::cout << "----------------------------------------------------\n";

    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns, all_peak_mem_bytes;
    
    // Load data once
    std::vector<std::vector<double>> all_columns = read_csv_columns(csv_file_path);
    if (all_columns.empty() || all_columns[0].empty()) {
         std::cerr << "Warning: Data loading issue for " << dataset_name << ". Skipping test." << std::endl;
        return;
    }
    int effective_num_rows = 0;
    all_columns = slice_columns_to_page(std::move(all_columns), effective_num_rows);
    long long original_size_bytes = (long long)effective_num_rows * expected_dim * sizeof(double);

    for (int run = 0; run < NUM_RUNS; ++run) {
        long long total_compressed_size_bytes = 0;
        
        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        if(run == 0) std::filesystem::create_directories(algo_output_dir);
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        uint64_t peak_mem_delta = measure_peak_private_delta_bytes([&]() {
            for (size_t i = 0; i < all_columns.size(); ++i) {
                std::vector<uint8_t> primary = compressor.encode(all_columns[i]);
                std::vector<uint8_t> final_data = sec_comp->compress(primary);
                total_compressed_size_bytes += final_data.size();
                std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
                write_binary_file(compressed_filename, final_data);
            }
        });
        auto end_full_encode = std::chrono::high_resolution_clock::now();
        
        // --- Decoding ---
        auto start_full_decode = std::chrono::high_resolution_clock::now();
        std::vector<std::vector<double>> decoded_columns;
        decoded_columns.reserve(all_columns.size());

        for (size_t i = 0; i < all_columns.size(); ++i) {
            std::string compressed_filename = algo_output_dir + "/" + dataset_name + "_col_" + std::to_string(i) + ".bin";
            auto file_data = read_binary_file(compressed_filename);

            // 1. Secondary Decompression
            auto primary_data = sec_comp->decompress(file_data);
            // 2. Primary Decompression
            decoded_columns.push_back(compressor.decode(primary_data));
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
        all_peak_mem_bytes.push_back(static_cast<long long>(peak_mem_delta));
    }

    // --- Stats ---
    double avg_compressed_size = std::accumulate(all_compressed_sizes.begin(), all_compressed_sizes.end(), 0.0) / NUM_RUNS;
    TimeStats enc_stats = calculate_stats(all_encoding_times_ns, effective_num_rows);
    TimeStats dec_stats = calculate_stats(all_decoding_times_ns, effective_num_rows);
    double avg_peak_mem = std::accumulate(all_peak_mem_bytes.begin(), all_peak_mem_bytes.end(), 0.0) / NUM_RUNS;
    
    double avg_bytes_per_point = avg_compressed_size / ((double)effective_num_rows * expected_dim);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes\n";
    std::cout << "Bytes per point:       " << avg_bytes_per_point << "\n";
    std::cout << "Enc (ns/pt): Avg=" << enc_stats.avg << ", Min=" << enc_stats.min << ", Max=" << enc_stats.max << "\n";
    std::cout << "Dec (ns/pt): Avg=" << dec_stats.avg << ", Min=" << dec_stats.min << ", Max=" << dec_stats.max << "\n";
    std::cout << "Peak Encode Memory:    " << avg_peak_mem << " bytes\n";
    std::cout << "====================================================\n" << std::endl;

    // Log to CSV
    log_result_to_csv(dataset_name, algorithm_name, expected_dim, effective_num_rows, original_size_bytes, avg_compressed_size,
                      enc_stats.avg, enc_stats.min, enc_stats.max,
                      dec_stats.avg, dec_stats.min, dec_stats.max);
    log_memory_to_csv(dataset_name, algorithm_name, effective_num_rows, avg_peak_mem);
    cleanup_method_output(output_dir, algorithm_name);
}


void run_multidim_reger_test(
    const std::string& dataset_name, int expected_dim, int expected_num_rows, 
    const std::string& csv_file_path, const std::string& output_dir,
    SecondaryCompressor* sec_comp = nullptr
) {
    NoneCompressor default_none;
    if (!sec_comp) sec_comp = &default_none;

    Reger compressor;
    std::string algorithm_name = compressor.get_name() + sec_comp->get_name();

    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ", Method: " << algorithm_name << " (Multi-dim)\n";
    std::cout << "----------------------------------------------------\n";

    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns, all_peak_mem_bytes;
    
    // Read data for processing
    std::vector<std::vector<std::string>> str_rows = read_csv_rows_as_strings(csv_file_path);
    if (str_rows.empty()) { std::cout << "No data, skipping.\n"; return; }
    int effective_num_rows = 0;
    str_rows = slice_string_rows_to_page(std::move(str_rows), effective_num_rows);
    long long original_size_bytes = (long long)effective_num_rows * expected_dim * sizeof(double);

    for (int run = 0; run < NUM_RUNS; ++run) {
        // --- ENCODING ---
        int num_rows = str_rows.size(); int num_cols = str_rows[0].size();
        std::vector<uint8_t> compressed_data;
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        uint64_t peak_mem_delta = measure_peak_private_delta_bytes([&]() {
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
            std::vector<uint8_t> primary = compressor.encode_multidim(scaled_rows, scaling_factors);
            compressed_data = sec_comp->compress(primary);
        });

        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        if(run == 0) std::filesystem::create_directories(algo_output_dir);
        std::string filename = algo_output_dir + "/" + dataset_name + ".bin";
        write_binary_file(filename, compressed_data);
        auto end_full_encode = std::chrono::high_resolution_clock::now();

        // --- DECODING ---
        auto start_full_decode = std::chrono::high_resolution_clock::now();
        auto file_data = read_binary_file(filename);
        
        // 1. Secondary Decompress
        auto primary_restored = sec_comp->decompress(file_data);
        // 2. Primary Decompress
        auto [decoded_scaled_rows, decoded_factors] = compressor.decode_multidim(primary_restored);

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
            int verify_rows = 0;
            original_rows = slice_rows_to_page(std::move(original_rows), verify_rows);
            verify_data_order_oblivious(original_rows, decoded_rows, dataset_name);
            std::cout << "Verification successful!\n";
        }

        all_compressed_sizes.push_back(compressed_data.size());
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count());
        all_peak_mem_bytes.push_back(static_cast<long long>(peak_mem_delta));
    }

    // --- Stats ---
    double avg_compressed_size = std::accumulate(all_compressed_sizes.begin(), all_compressed_sizes.end(), 0.0) / NUM_RUNS;
    TimeStats enc_stats = calculate_stats(all_encoding_times_ns, effective_num_rows);
    TimeStats dec_stats = calculate_stats(all_decoding_times_ns, effective_num_rows);
    double avg_peak_mem = std::accumulate(all_peak_mem_bytes.begin(), all_peak_mem_bytes.end(), 0.0) / NUM_RUNS;
    
    double avg_bytes_per_point = avg_compressed_size / ((double)effective_num_rows * expected_dim);
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes\n";
    std::cout << "Bytes per point:       " << avg_bytes_per_point << "\n";
    std::cout << "Enc (ns/row): Avg=" << enc_stats.avg << ", Min=" << enc_stats.min << ", Max=" << enc_stats.max << "\n";
    std::cout << "Dec (ns/row): Avg=" << dec_stats.avg << ", Min=" << dec_stats.min << ", Max=" << dec_stats.max << "\n";
    std::cout << "Peak Encode Memory:    " << avg_peak_mem << " bytes\n";
    std::cout << "====================================================\n\n";

    log_result_to_csv(dataset_name, algorithm_name, expected_dim, effective_num_rows, original_size_bytes, avg_compressed_size,
                      enc_stats.avg, enc_stats.min, enc_stats.max,
                      dec_stats.avg, dec_stats.min, dec_stats.max);
    log_memory_to_csv(dataset_name, algorithm_name, effective_num_rows, avg_peak_mem);
    cleanup_method_output(output_dir, algorithm_name);
}

template<typename ClusterCompressor>
void run_cluster_encoder_test(
    const std::string& dataset_name, int dim, int num_rows, 
    const std::string& csv_file_path, const std::string& output_dir,
    SecondaryCompressor* sec_comp = nullptr,
    int dimension_group_size = 0
) {
    NoneCompressor default_none;
    if (!sec_comp) sec_comp = &default_none;

    auto apply_group_size = [&](auto& comp) {
        using T = std::decay_t<decltype(comp)>;
        if constexpr (std::is_same_v<T, KClusterSub> || std::is_same_v<T, AClusterSub>) {
            if (dimension_group_size > 0) comp.set_dimension_group_size(dimension_group_size);
        }
    };
    ClusterCompressor compressor;
    apply_group_size(compressor);
    std::string algorithm_name = compressor.get_name() + sec_comp->get_name(); // Combine names

    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name << ", Method: " << algorithm_name << " (Multi-dim)\n";
    std::cout << "----------------------------------------------------\n";

    std::vector<long long> all_compressed_sizes, all_encoding_times_ns, all_decoding_times_ns, all_peak_mem_bytes;
    std::vector<std::vector<std::string>> page_rows_str = read_csv_rows_as_strings(csv_file_path);
    if (page_rows_str.empty()) { std::cout << "No data, skipping.\n"; return; }
    int effective_num_rows = 0;
    page_rows_str = slice_string_rows_to_page(std::move(page_rows_str), effective_num_rows);
    std::string cache_dir = output_dir + "/_page_cache";
    std::filesystem::create_directories(cache_dir);
    std::string page_csv_path = cache_dir + "/" + dataset_name + "_cluster_page.csv";
    write_page_csv_from_string_rows(page_rows_str, page_csv_path);
    long long original_size_bytes = (long long)effective_num_rows * dim * sizeof(double);
    int k = 100, page_size = 10000, pack_size = 10, block_size = 10;

    for(int run = 0; run < NUM_RUNS; ++run) {
        
        auto start_full_encode = std::chrono::high_resolution_clock::now();
        int mutable_dim = dim;
        std::vector<uint8_t> final_data;
        std::vector<uint8_t> primary_data;
        GlobalDiagnostics captured_diag;
        bool has_diag = false;
        uint64_t peak_mem_delta = measure_peak_private_delta_bytes([&]() {
            using C = std::decay_t<ClusterCompressor>;
            if constexpr (std::is_same_v<C, KCluster>) {
                auto with_diag = compressor.encode_multidim_with_diag(page_csv_path, mutable_dim, k, page_size, pack_size, block_size);
                primary_data = std::move(with_diag.first);
                captured_diag = std::move(with_diag.second);
                has_diag = true;
            } else if constexpr (std::is_same_v<C, ACluster>) {
                auto with_diag = compressor.encode_multidim_with_diag(page_csv_path, mutable_dim, k, page_size, pack_size, block_size);
                primary_data = std::move(with_diag.first);
                captured_diag = std::move(with_diag.second);
                has_diag = true;
            } else {
                primary_data = compressor.encode_multidim(page_csv_path, mutable_dim, k, page_size, pack_size, block_size);
            }
            final_data = sec_comp->compress(primary_data);
        });

        std::string algo_output_dir = output_dir + "/" + algorithm_name;
        if (run == 0) std::filesystem::create_directories(algo_output_dir);
        std::string filename = algo_output_dir + "/" + dataset_name + ".bin";
        write_binary_file(filename, final_data);
        if (run == 0) {
            using C = std::decay_t<ClusterCompressor>;
            if constexpr (std::is_same_v<C, KClusterSub> || std::is_same_v<C, AClusterSub>) {
                int sub_mode = -1;
                int sub_strategy = -1;
                if (final_data.size() >= 5) {
                    if ((final_data[0] == 'K' || final_data[0] == 'A') &&
                        final_data[1] == 'S' && final_data[2] == 'U' && final_data[3] == 'B') {
                        sub_mode = static_cast<int>(final_data[4]);
                        if (final_data.size() >= 6) sub_strategy = static_cast<int>(final_data[5]);
                    }
                }
                std::ofstream mode_csv(output_dir + "/sub_mode_usage.csv", std::ios::app);
                mode_csv << dataset_name << "," << algorithm_name << "," << sub_mode << "\n";
                std::ofstream strategy_csv(output_dir + "/sub_strategy_usage.csv", std::ios::app);
                strategy_csv << dataset_name << "," << algorithm_name << "," << sub_strategy << "\n";
                std::cout << "Sub mode used: " << sub_mode << " (always 1), strategy: " << sub_strategy << " (0=split,1=base-carry)\n";

                double header_bits_per_row = 0.0, base_bits_per_row = 0.0, residual_bits_per_row = 0.0, id_extra_bits_per_row = 0.0;
                if (parse_sub_cost_metrics(primary_data, header_bits_per_row, base_bits_per_row, residual_bits_per_row, id_extra_bits_per_row)) {
                    std::ofstream sub_cost_csv(output_dir + "/method_cost_components_by_dataset.csv", std::ios::app);
                    sub_cost_csv << dataset_name << "," << algorithm_name << ","
                                 << std::fixed << std::setprecision(4) << header_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << base_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << residual_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << id_extra_bits_per_row << "\n";
                }
            } else if constexpr (std::is_same_v<C, KCluster> || std::is_same_v<C, ACluster>) {
                if (has_diag && effective_num_rows > 0) {
                    double row_n = static_cast<double>(effective_num_rows);
                    double header_bits_per_row = static_cast<double>(captured_diag.total_sizes.header_bits + captured_diag.total_sizes.frequencies_bits) / row_n;
                    double base_bits_per_row = static_cast<double>(captured_diag.total_sizes.medoids_bits) / row_n;
                    double residual_bits_per_row = static_cast<double>(captured_diag.total_sizes.residuals_bits) / row_n;
                    std::ofstream sub_cost_csv(output_dir + "/method_cost_components_by_dataset.csv", std::ios::app);
                    sub_cost_csv << dataset_name << "," << algorithm_name << ","
                                 << std::fixed << std::setprecision(4) << header_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << base_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << residual_bits_per_row << ","
                                 << std::fixed << std::setprecision(4) << 0.0 << "\n";
                }
            }
        }
        auto end_full_encode = std::chrono::high_resolution_clock::now();

        auto start_full_decode = std::chrono::high_resolution_clock::now();
        auto read_data = read_binary_file(filename);

        // 1. Secondary Decode
        auto primary_data_restored = sec_comp->decompress(read_data);
        // 2. Primary Decode
        auto decoded_rows = compressor.decode_multidim(primary_data_restored);
        
        auto end_full_decode = std::chrono::high_resolution_clock::now();

        if (run == 0) {
            std::vector<std::vector<double>> original_rows = read_csv_rows(page_csv_path);
            verify_data_order_oblivious(original_rows, decoded_rows, dataset_name);
            std::cout << "Verification successful!\n";
        }

        all_compressed_sizes.push_back(final_data.size());
        all_encoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_encode - start_full_encode).count());
        all_decoding_times_ns.push_back(std::chrono::duration_cast<std::chrono::nanoseconds>(end_full_decode - start_full_decode).count());
        all_peak_mem_bytes.push_back(static_cast<long long>(peak_mem_delta));
    }

    // --- Stats ---
    double avg_compressed_size = std::accumulate(all_compressed_sizes.begin(), all_compressed_sizes.end(), 0.0) / NUM_RUNS;
    TimeStats enc_stats = calculate_stats(all_encoding_times_ns, effective_num_rows);
    TimeStats dec_stats = calculate_stats(all_decoding_times_ns, effective_num_rows);
    double avg_peak_mem = std::accumulate(all_peak_mem_bytes.begin(), all_peak_mem_bytes.end(), 0.0) / NUM_RUNS;
    
    double avg_bytes_per_point = avg_compressed_size / ((double)effective_num_rows * dim);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Byte Cost:       " << avg_compressed_size << " bytes\n";
    std::cout << "Bytes per point:       " << avg_bytes_per_point << "\n";
    std::cout << "Enc (ns/row): Avg=" << enc_stats.avg << ", Min=" << enc_stats.min << ", Max=" << enc_stats.max << "\n";
    std::cout << "Dec (ns/row): Avg=" << dec_stats.avg << ", Min=" << dec_stats.min << ", Max=" << dec_stats.max << "\n";
    std::cout << "Peak Encode Memory:    " << avg_peak_mem << " bytes\n";
    std::cout << "====================================================\n\n";

    log_result_to_csv(dataset_name, algorithm_name, dim, effective_num_rows, original_size_bytes, avg_compressed_size,
                      enc_stats.avg, enc_stats.min, enc_stats.max,
                      dec_stats.avg, dec_stats.min, dec_stats.max);
    log_memory_to_csv(dataset_name, algorithm_name, effective_num_rows, avg_peak_mem);
    cleanup_method_output(output_dir, algorithm_name);
    std::filesystem::remove(page_csv_path);
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
    int k, int page_size, int pack_size) // Parameters for the experiment
{
    ClusterCompressor compressor;
    std::string algorithm_name = compressor.get_name();
    
    std::cout << "====================================================\n";
    std::cout << "Dataset: " << dataset_name 
              << ", Method: " << algorithm_name
              << ", k: " << k 
              << ", page_size: " << page_size 
              << ", pack_size: " << pack_size << "\n";
    std::cout << "----------------------------------------------------\n";

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
                k, fixed_page_size, 10
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
                fixed_k_for_kcluster, page_size, 10
            );

            // Run for ACluster (k value is ignored, pass 0 as placeholder)
            run_cluster_encoder_experiment_instance<ACluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                0, page_size, 10
            );
        }
    }
    std::cout << "\n\n### PAGE SIZE SENSITIVITY EXPERIMENT COMPLETE ###\n\n";
}

void run_pack_size_sensitivity_experiment(const std::string& file_root, const std::string& output_dir) {
    std::cout << "\n\n\n#######################################################\n";
    std::cout << "###   STARTING PACK SIZE SENSITIVITY EXPERIMENT   ###\n";
    std::cout << "#######################################################\n\n\n";
    
    std::map<std::string, std::pair<int, int>> datasets_to_test = {
        {"jinfeng",      {27, 76512}},
        {"Crop",    {46, 24000}},
        {"gas",     {64, 13910}}
    };

    std::vector<int> pack_size_values = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    
    int fixed_k = 100;       
    int fixed_page_size = 10000; 

    for (const auto& ds_pair : datasets_to_test) {
        const std::string& dataset_name = ds_pair.first;
        int dim = ds_pair.second.first;
        int num_rows = ds_pair.second.second;
        std::string file_path = file_root + dataset_name + ".csv";

        for (int pack_size : pack_size_values) {
            // Run for KCluster
            run_cluster_encoder_experiment_instance<KCluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                fixed_k, fixed_page_size, pack_size
            );

            // Run for ACluster (k is ignored but we pass fixed_k)
            run_cluster_encoder_experiment_instance<ACluster>(
                dataset_name, dim, num_rows, file_path, output_dir,
                fixed_k, fixed_page_size, pack_size
            );
        }
    }
    std::cout << "\n\n### PACK SIZE SENSITIVITY EXPERIMENT COMPLETE ###\n\n";
}

struct ClusterReferenceResult {
    std::string algorithm;
    int k;
    long long references;
    double compression_ratio;
    double enc_ns_per_row;
    double dec_ns_per_row;
    long long compressed_bytes;
    int rows;
};

void log_progress(const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tmv{};
    localtime_s(&tmv, &t);
    std::cout << "[" << std::put_time(&tmv, "%H:%M:%S") << "] " << msg << std::endl;
}

template<typename ClusterCompressor>
ClusterReferenceResult run_cluster_reference_once(
    const std::string& algorithm_name,
    int k,
    const std::string& dataset_name,
    int dim,
    const std::string& csv_file_path
) {
    ClusterCompressor compressor;
    int mutable_dim = dim;
    int page_size = 10000, pack_size = 10, block_size = 10;

    log_progress("START " + algorithm_name + " dataset=" + dataset_name + " k=" + std::to_string(k));
    auto start_enc = std::chrono::high_resolution_clock::now();
    auto encoded_with_diag = compressor.encode_multidim_with_diag(csv_file_path, mutable_dim, k, page_size, pack_size, block_size);
    auto end_enc = std::chrono::high_resolution_clock::now();

    const auto& compressed_data = encoded_with_diag.first;
    const auto& diag = encoded_with_diag.second;

    auto start_dec = std::chrono::high_resolution_clock::now();
    auto decoded_rows = compressor.decode_multidim(compressed_data);
    auto end_dec = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<double>> original_rows = read_csv_rows(csv_file_path);
    verify_data_order_oblivious(original_rows, decoded_rows, dataset_name + "_" + algorithm_name);

    int rows = static_cast<int>(original_rows.size());
    long long original_size_bytes = static_cast<long long>(rows) * dim * sizeof(double);
    double ratio = compressed_data.empty() ? 0.0 : static_cast<double>(original_size_bytes) / static_cast<double>(compressed_data.size());
    double enc_ns_per_row = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count()) / rows : 0.0;
    double dec_ns_per_row = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count()) / rows : 0.0;
    log_progress("DONE  " + algorithm_name + " dataset=" + dataset_name + " k=" + std::to_string(k) +
                 " ratio=" + std::to_string(ratio) + " refs=" + std::to_string(diag.total_reference_points));

    return {
        algorithm_name,
        k,
        diag.total_reference_points,
        ratio,
        enc_ns_per_row,
        dec_ns_per_row,
        static_cast<long long>(compressed_data.size()),
        rows
    };
}

void run_cluster_reference_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    std::ofstream detail(output_dir + "/cluster_reference_results.csv");
    detail << "Dataset,Algorithm,K,Rows,References,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row,CompressedBytes\n";

    std::ofstream matched(output_dir + "/cluster_reference_match_summary.csv");
    matched << "Dataset,ACluster_Ratio,ACluster_Enc_ns,ACluster_Dec_ns,ACluster_Refs,BestK,BestK_Ratio,BestK_Enc_ns,BestK_Dec_ns,BestK_Refs,RatioGap,RefDelta(K-A)\n";

    std::vector<int> k_candidates = {100, 200, 500};
    double sum_a_ratio = 0.0, sum_a_enc = 0.0, sum_a_dec = 0.0, sum_a_refs = 0.0;
    double sum_k_ratio = 0.0, sum_k_enc = 0.0, sum_k_dec = 0.0, sum_k_refs = 0.0;
    double sum_gap = 0.0;
    int matched_count = 0;

    for (size_t i = 0; i < datasets.size(); ++i) {
        std::string dataset_name = datasets[i];
        std::string file_path = file_root + dataset_name + ".csv";
        int dim = dims[i];

        log_progress("DATASET BEGIN " + dataset_name);
        std::cout << "====================================================\n";
        std::cout << "Reference Experiment Dataset: " << dataset_name << "\n";
        std::cout << "----------------------------------------------------\n";

        auto a_res = run_cluster_reference_once<ACluster>("ACluster", 100, dataset_name, dim, file_path);
        detail << dataset_name << "," << a_res.algorithm << "," << a_res.k << "," << a_res.rows << "," << a_res.references << ","
               << std::fixed << std::setprecision(4) << a_res.compression_ratio << ","
               << std::setprecision(2) << a_res.enc_ns_per_row << "," << a_res.dec_ns_per_row << ","
               << a_res.compressed_bytes << "\n";

        double best_gap = std::numeric_limits<double>::max();
        ClusterReferenceResult best_k_res{};
        for (int k : k_candidates) {
            auto k_res = run_cluster_reference_once<KCluster>("KCluster", k, dataset_name, dim, file_path);
            detail << dataset_name << "," << k_res.algorithm << "," << k_res.k << "," << k_res.rows << "," << k_res.references << ","
                   << std::fixed << std::setprecision(4) << k_res.compression_ratio << ","
                   << std::setprecision(2) << k_res.enc_ns_per_row << "," << k_res.dec_ns_per_row << ","
                   << k_res.compressed_bytes << "\n";

            double gap = std::abs(k_res.compression_ratio - a_res.compression_ratio);
            if (gap < best_gap) {
                best_gap = gap;
                best_k_res = k_res;
            }
        }

        matched << dataset_name << ","
                << std::fixed << std::setprecision(4) << a_res.compression_ratio << ","
                << std::setprecision(2) << a_res.enc_ns_per_row << "," << a_res.dec_ns_per_row << ","
                << a_res.references << ","
                << best_k_res.k << ","
                << std::fixed << std::setprecision(4) << best_k_res.compression_ratio << ","
                << std::setprecision(2) << best_k_res.enc_ns_per_row << "," << best_k_res.dec_ns_per_row << ","
                << best_k_res.references << ","
                << std::fixed << std::setprecision(6) << best_gap << ","
                << (best_k_res.references - a_res.references) << "\n";

        sum_a_ratio += a_res.compression_ratio;
        sum_a_enc += a_res.enc_ns_per_row;
        sum_a_dec += a_res.dec_ns_per_row;
        sum_a_refs += static_cast<double>(a_res.references);
        sum_k_ratio += best_k_res.compression_ratio;
        sum_k_enc += best_k_res.enc_ns_per_row;
        sum_k_dec += best_k_res.dec_ns_per_row;
        sum_k_refs += static_cast<double>(best_k_res.references);
        sum_gap += best_gap;
        matched_count++;

        std::cout << "ACluster ratio=" << std::fixed << std::setprecision(4) << a_res.compression_ratio
                  << ", refs=" << a_res.references
                  << " | Best K=" << best_k_res.k
                  << ", ratio=" << best_k_res.compression_ratio
                  << ", refs=" << best_k_res.references << "\n";
        std::cout << "====================================================\n\n";
        log_progress("DATASET END   " + dataset_name);
    }

    std::ofstream avg(output_dir + "/cluster_reference_match_avg.csv");
    avg << "Metric,ACluster,KCluster(BestMatch),RatioGap\n";
    if (matched_count > 0) {
        avg << "CompressionRatio," << std::fixed << std::setprecision(4) << (sum_a_ratio / matched_count) << ","
            << (sum_k_ratio / matched_count) << "," << (sum_gap / matched_count) << "\n";
        avg << "EncTime_ns_per_row," << std::fixed << std::setprecision(2) << (sum_a_enc / matched_count) << ","
            << (sum_k_enc / matched_count) << ",\n";
        avg << "DecTime_ns_per_row," << std::fixed << std::setprecision(2) << (sum_a_dec / matched_count) << ","
            << (sum_k_dec / matched_count) << ",\n";
        avg << "ReferencePoints," << std::fixed << std::setprecision(2) << (sum_a_refs / matched_count) << ","
            << (sum_k_refs / matched_count) << ",\n";
    }
    avg.flush();
}

void run_acluster_three_mode_reference_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    std::vector<RowOrderMode> modes = {
        RowOrderMode::Original,
        RowOrderMode::RandomShuffle,
        RowOrderMode::SortByFirstColumn
    };
    std::string cache_dir = output_dir + "/_scenario_cache_acluster";
    std::filesystem::create_directories(cache_dir);

    std::ofstream detail(output_dir + "/acluster_three_modes_with_refs.csv");
    detail << "Scenario,Dataset,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row,ReferencePoints,Rows,CompressedBytes\n";

    std::map<std::string, std::vector<std::tuple<double, double, double, double>>> scenario_agg;
    ACluster compressor;
    int k_placeholder = 100, page_size = 10000, pack_size = 10, block_size = 10;

    for (auto mode : modes) {
        std::string sc = scenario_name(mode);
        for (size_t i = 0; i < datasets.size(); ++i) {
            const std::string& dataset = datasets[i];
            int dim = dims[i];
            std::string src = file_root + dataset + ".csv";
            std::string path = prepare_dataset_for_mode(src, dataset, mode, cache_dir);

            log_progress("ACluster 3Mode START scenario=" + sc + " dataset=" + dataset);
            int mutable_dim = dim;
            auto start_enc = std::chrono::high_resolution_clock::now();
            auto encoded_with_diag = compressor.encode_multidim_with_diag(path, mutable_dim, k_placeholder, page_size, pack_size, block_size);
            auto end_enc = std::chrono::high_resolution_clock::now();

            const auto& compressed = encoded_with_diag.first;
            const auto& diag = encoded_with_diag.second;

            auto start_dec = std::chrono::high_resolution_clock::now();
            auto decoded = compressor.decode_multidim(compressed);
            auto end_dec = std::chrono::high_resolution_clock::now();

            auto original_rows = read_csv_rows(path);
            verify_data_order_oblivious(original_rows, decoded, dataset + "_" + sc + "_ACluster");

            int rows = static_cast<int>(original_rows.size());
            long long original_size = static_cast<long long>(rows) * dim * sizeof(double);
            double comp_ratio = compressed.empty() ? 0.0 : static_cast<double>(original_size) / compressed.size();
            double enc_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count()) / rows : 0.0;
            double dec_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count()) / rows : 0.0;
            double refs = static_cast<double>(diag.total_reference_points);

            detail << sc << "," << dataset << ","
                   << std::fixed << std::setprecision(4) << comp_ratio << ","
                   << std::setprecision(2) << enc_ns << "," << dec_ns << ","
                   << static_cast<long long>(refs) << "," << rows << "," << compressed.size() << "\n";
            detail.flush();
            scenario_agg[sc].push_back({comp_ratio, enc_ns, dec_ns, refs});
            log_progress("ACluster 3Mode DONE  scenario=" + sc + " dataset=" + dataset +
                         " ratio=" + std::to_string(comp_ratio) + " refs=" + std::to_string(static_cast<long long>(refs)));

            if (mode != RowOrderMode::Original) std::filesystem::remove(path);
        }
    }

    std::ofstream summary(output_dir + "/acluster_three_modes_avg.csv");
    summary << "Scenario,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row,AvgReferencePoints\n";
    for (auto mode : modes) {
        std::string sc = scenario_name(mode);
        const auto& vals = scenario_agg[sc];
        double s_ratio = 0.0, s_enc = 0.0, s_dec = 0.0, s_ref = 0.0;
        for (const auto& t : vals) {
            s_ratio += std::get<0>(t);
            s_enc += std::get<1>(t);
            s_dec += std::get<2>(t);
            s_ref += std::get<3>(t);
        }
        double n = vals.empty() ? 1.0 : static_cast<double>(vals.size());
        summary << sc << ","
                << std::fixed << std::setprecision(4) << (s_ratio / n) << ","
                << std::setprecision(2) << (s_enc / n) << "," << (s_dec / n) << ","
                << std::setprecision(2) << (s_ref / n) << "\n";
    }
    summary.flush();
    std::filesystem::remove_all(cache_dir);
}

std::string make_outlier_dataset(
    const std::string& source_csv_path,
    const std::string& dataset_name,
    double outlier_ratio,
    const std::string& cache_dir,
    int& outlier_count
) {
    auto rows = read_csv_rows_as_strings(source_csv_path);
    if (rows.empty()) {
        outlier_count = 0;
        return source_csv_path;
    }
    int n = static_cast<int>(rows.size());
    outlier_count = static_cast<int>(std::round(n * outlier_ratio));
    if (outlier_count <= 0) outlier_count = 1;
    if (outlier_count > n) outlier_count = n;

    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 rng(static_cast<uint32_t>(20260406 + std::hash<std::string>{}(dataset_name) + static_cast<size_t>(outlier_ratio * 1000000)));
    std::shuffle(indices.begin(), indices.end(), rng);

    for (int i = 0; i < outlier_count; ++i) {
        int row_idx = indices[i];
        if (rows[row_idx].empty()) continue;
        double base = 0.0;
        try {
            base = std::stod(rows[row_idx][0]);
        } catch (...) {
            base = 0.0;
        }
        double bumped = base + 100.0;
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << bumped;
        rows[row_idx][0] = oss.str();
    }

    std::filesystem::create_directories(cache_dir);
    std::ostringstream name;
    name << dataset_name << "__outlier_" << std::fixed << std::setprecision(3) << (outlier_ratio * 100.0) << ".csv";
    std::string out_path = cache_dir + "/" + name.str();
    return write_page_csv_from_string_rows(rows, out_path);
}

void run_acluster_outlier_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    std::vector<double> outlier_ratios = {0.0, 0.001, 0.005, 0.01, 0.05};
    std::map<double, std::vector<std::tuple<double, double, double, double>>> agg;
    std::string cache_dir = output_dir + "/_outlier_cache";
    std::filesystem::create_directories(cache_dir);

    std::ofstream detail(output_dir + "/acluster_outlier_detail.csv");
    detail << "Dataset,OutlierRatio,OutlierCount,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row,References,Rows,CompressedBytes\n";
    detail.flush();

    ACluster compressor;
    int page_size = 10000, pack_size = 10, block_size = 10;
    int k_placeholder = 100;

    for (size_t i = 0; i < datasets.size(); ++i) {
        const std::string& dataset = datasets[i];
        int dim = dims[i];
        std::string source_path = file_root + dataset + ".csv";

        for (double ratio : outlier_ratios) {
            std::string outlier_path;
            try {
                int outlier_count = 0;
                outlier_path = make_outlier_dataset(source_path, dataset, ratio, cache_dir, outlier_count);

                log_progress("OUTLIER START dataset=" + dataset + " ratio=" + std::to_string(ratio));
                int mutable_dim = dim;
                auto start_enc = std::chrono::high_resolution_clock::now();
                auto encoded_with_diag = compressor.encode_multidim_with_diag(outlier_path, mutable_dim, k_placeholder, page_size, pack_size, block_size);
                auto end_enc = std::chrono::high_resolution_clock::now();

                const auto& compressed = encoded_with_diag.first;
                const auto& diag = encoded_with_diag.second;

                auto start_dec = std::chrono::high_resolution_clock::now();
                auto decoded = compressor.decode_multidim(compressed);
                auto end_dec = std::chrono::high_resolution_clock::now();

                auto original_rows = read_csv_rows(outlier_path);

                int rows = static_cast<int>(original_rows.size());
                long long original_size = static_cast<long long>(rows) * dim * sizeof(double);
                double comp_ratio = compressed.empty() ? 0.0 : static_cast<double>(original_size) / compressed.size();
                double enc_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count()) / rows : 0.0;
                double dec_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count()) / rows : 0.0;

                detail << dataset << "," << std::fixed << std::setprecision(3) << ratio << "," << outlier_count << ","
                       << std::fixed << std::setprecision(4) << comp_ratio << ","
                       << std::setprecision(2) << enc_ns << "," << dec_ns << ","
                       << diag.total_reference_points << "," << rows << "," << compressed.size() << "\n";
                detail.flush();
                agg[ratio].push_back({comp_ratio, enc_ns, dec_ns, static_cast<double>(diag.total_reference_points)});
                log_progress("OUTLIER DONE  dataset=" + dataset + " ratio=" + std::to_string(ratio) +
                             " comp=" + std::to_string(comp_ratio) + " enc=" + std::to_string(enc_ns));
            } catch (const std::exception& e) {
                log_progress("OUTLIER FAIL  dataset=" + dataset + " ratio=" + std::to_string(ratio) + " err=" + e.what());
            }
            if (!outlier_path.empty() && outlier_path != source_path) {
                std::filesystem::remove(outlier_path);
            }
        }
    }

    std::ofstream summary(output_dir + "/acluster_outlier_summary_avg.csv");
    summary << "OutlierRatio,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row,AvgReferencePoints\n";
    for (double ratio : outlier_ratios) {
        const auto& v = agg[ratio];
        double s_comp = 0.0, s_enc = 0.0, s_dec = 0.0, s_ref = 0.0;
        for (const auto& t : v) {
            s_comp += std::get<0>(t);
            s_enc += std::get<1>(t);
            s_dec += std::get<2>(t);
            s_ref += std::get<3>(t);
        }
        double n = v.empty() ? 1.0 : static_cast<double>(v.size());
        summary << std::fixed << std::setprecision(3) << ratio << ","
                << std::setprecision(4) << (s_comp / n) << ","
                << std::setprecision(2) << (s_enc / n) << "," << (s_dec / n) << ","
                << std::setprecision(2) << (s_ref / n) << "\n";
    }
    summary.flush();
    std::filesystem::remove_all(cache_dir);
}

int write_first_k_columns_csv(
    const std::string& source_csv_path,
    const std::string& out_csv_path,
    int k
) {
    auto rows = read_csv_rows_as_strings(source_csv_path);
    std::ofstream out(out_csv_path);
    if (!out.is_open()) throw std::runtime_error("Cannot write reduced csv: " + out_csv_path);
    int row_count = 0;
    for (const auto& row : rows) {
        if (row.empty()) continue;
        int take = std::min(k, static_cast<int>(row.size()));
        for (int i = 0; i < take; ++i) {
            if (i) out << ",";
            out << row[i];
        }
        for (int i = take; i < k; ++i) {
            if (i) out << ",";
            out << "0";
        }
        out << "\n";
        row_count++;
    }
    return row_count;
}

void run_dimensionality_impact_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    init_results_csv(output_dir + "/dimensionality_raw_results.csv");
    currentScenario = "dimensionality";
    size_t start_idx = loggedResults.size();
    std::string cache_dir = output_dir + "/_dim_cache";
    std::filesystem::create_directories(cache_dir);
    NoneCompressor none_comp;

    auto safe_run = [&](const std::string& method_name, const std::function<void()>& fn) {
        try {
            fn();
        } catch (const std::exception& e) {
            std::cerr << "Method failed [" << method_name << "]: " << e.what() << std::endl;
        } catch (...) {
            std::cerr << "Method failed [" << method_name << "]: unknown error" << std::endl;
        }
    };

    for (size_t i = 0; i < datasets.size(); ++i) {
        if (dims[i] <= 1) continue;
        std::string src = file_root + datasets[i] + ".csv";
        int max_k = std::min(5, dims[i]);
        for (int k = 1; k <= max_k; ++k) {
            std::string reduced_path = cache_dir + "/" + datasets[i] + "_first" + std::to_string(k) + ".csv";
            int rows = write_first_k_columns_csv(src, reduced_path, k);
            if (rows <= 0) continue;
            log_progress("DIM START dataset=" + datasets[i] + " k=" + std::to_string(k));
            safe_run("Chimp", [&]() { run_columnar_encoder_test<Chimp>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Elf", [&]() { run_columnar_encoder_test<Elf>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Gorilla", [&]() { run_columnar_encoder_test<Gorilla>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Rle", [&]() { run_columnar_encoder_test<Rle>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Huffman", [&]() { run_columnar_encoder_test<Huffman>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("ALP", [&]() { run_columnar_encoder_test<Alp>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Reger", [&]() { run_multidim_reger_test(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("KCluster", [&]() { run_cluster_encoder_test<KCluster>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("ACluster", [&]() { run_cluster_encoder_test<ACluster>(datasets[i], k, rows, reduced_path, output_dir, &none_comp); });
            safe_run("Vortex+Rle", [&]() { run_vortex_test(datasets[i], k, rows, reduced_path, output_dir); });
            std::filesystem::remove(reduced_path);
            log_progress("DIM DONE  dataset=" + datasets[i] + " k=" + std::to_string(k));
        }
    }

    std::ofstream per_dataset(output_dir + "/dimensionality_by_dataset.csv");
    per_dataset << "Dataset,UsedDims,Algorithm,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row\n";

    std::map<std::string, std::vector<LoggedResultRecord>> grouped;
    for (size_t i = start_idx; i < loggedResults.size(); ++i) {
        const auto& r = loggedResults[i];
        if (r.scenario != "dimensionality") continue;
        per_dataset << r.dataset << "," << r.dim << "," << r.algorithm << ","
                    << std::fixed << std::setprecision(4) << r.compression_ratio << ","
                    << std::setprecision(2) << r.enc_avg_ns << "," << r.dec_avg_ns << "\n";
        grouped[std::to_string(r.dim) + "|" + r.algorithm].push_back(r);
    }
    per_dataset.flush();

    std::ofstream cross_avg(output_dir + "/dimensionality_cross_dataset_avg.csv");
    cross_avg << "UsedDims,Algorithm,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row\n";
    for (const auto& kv : grouped) {
        const auto& vec = kv.second;
        double sum_ratio = 0.0, sum_enc = 0.0, sum_dec = 0.0;
        for (const auto& r : vec) {
            sum_ratio += r.compression_ratio;
            sum_enc += r.enc_avg_ns;
            sum_dec += r.dec_avg_ns;
        }
        double n = vec.empty() ? 1.0 : static_cast<double>(vec.size());
        size_t split = kv.first.find('|');
        std::string dim_str = kv.first.substr(0, split);
        std::string algo = kv.first.substr(split + 1);
        cross_avg << dim_str << "," << algo << ","
                  << std::fixed << std::setprecision(4) << (sum_ratio / n) << ","
                  << std::setprecision(2) << (sum_enc / n) << "," << (sum_dec / n) << "\n";
    }
    cross_avg.flush();
    std::filesystem::remove_all(cache_dir);
}

void run_acluster_adjacent_reference_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    std::vector<int> adjacent_values = {1, 2, 3};
    std::map<int, std::vector<std::tuple<double, double, double, double>>> agg;
    std::ofstream detail(output_dir + "/acluster_adjacent_reference_detail.csv");
    detail << "Dataset,AdjacentReference,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row,ReferencePoints,Rows,CompressedBytes\n";

    ACluster compressor;
    int page_size = 10000, pack_size = 10, block_size = 10, k_placeholder = 100;

    for (size_t i = 0; i < datasets.size(); ++i) {
        std::string dataset = datasets[i];
        std::string path = file_root + dataset + ".csv";
        int dim = dims[i];
        for (int adj : adjacent_values) {
            log_progress("ADJ START dataset=" + dataset + " adjacent_reference=" + std::to_string(adj));
            int mutable_dim = dim;
            auto start_enc = std::chrono::high_resolution_clock::now();
            auto encoded_with_diag = compressor.encode_multidim_with_diag(path, mutable_dim, k_placeholder, page_size, pack_size, block_size, adj);
            auto end_enc = std::chrono::high_resolution_clock::now();

            const auto& compressed = encoded_with_diag.first;
            const auto& diag = encoded_with_diag.second;
            auto start_dec = std::chrono::high_resolution_clock::now();
            auto decoded = compressor.decode_multidim(compressed);
            auto end_dec = std::chrono::high_resolution_clock::now();

            auto original_rows = read_csv_rows(path);
            verify_data_order_oblivious(original_rows, decoded, dataset + "_adj_" + std::to_string(adj));

            int rows = static_cast<int>(original_rows.size());
            long long original_size = static_cast<long long>(rows) * mutable_dim * sizeof(double);
            double comp_ratio = compressed.empty() ? 0.0 : static_cast<double>(original_size) / compressed.size();
            double enc_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_enc - start_enc).count()) / rows : 0.0;
            double dec_ns = rows > 0 ? static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(end_dec - start_dec).count()) / rows : 0.0;
            double refs = static_cast<double>(diag.total_reference_points);

            detail << dataset << "," << adj << ","
                   << std::fixed << std::setprecision(4) << comp_ratio << ","
                   << std::setprecision(2) << enc_ns << "," << dec_ns << ","
                   << static_cast<long long>(refs) << "," << rows << "," << compressed.size() << "\n";
            detail.flush();
            agg[adj].push_back({comp_ratio, enc_ns, dec_ns, refs});
            log_progress("ADJ DONE  dataset=" + dataset + " adjacent_reference=" + std::to_string(adj) +
                         " ratio=" + std::to_string(comp_ratio) + " refs=" + std::to_string(static_cast<long long>(refs)));
        }
    }

    std::ofstream summary(output_dir + "/acluster_adjacent_reference_summary_avg.csv");
    summary << "AdjacentReference,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row,AvgReferencePoints\n";
    for (int adj : adjacent_values) {
        const auto& v = agg[adj];
        double s_ratio = 0.0, s_enc = 0.0, s_dec = 0.0, s_ref = 0.0;
        for (const auto& t : v) {
            s_ratio += std::get<0>(t);
            s_enc += std::get<1>(t);
            s_dec += std::get<2>(t);
            s_ref += std::get<3>(t);
        }
        double n = v.empty() ? 1.0 : static_cast<double>(v.size());
        summary << adj << ","
                << std::fixed << std::setprecision(4) << (s_ratio / n) << ","
                << std::setprecision(2) << (s_enc / n) << "," << (s_dec / n) << ","
                << std::setprecision(2) << (s_ref / n) << "\n";
    }
    summary.flush();
}

void run_subspace_grouping_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims,
    int dimension_group_size
) {
    init_results_csv(output_dir + "/subspace_group_raw_results.csv");
    std::ofstream method_cost_init(output_dir + "/method_cost_components_by_dataset.csv");
    method_cost_init << "Dataset,Algorithm,HeaderBitsPerRow,BaseBitsPerRow,ResidualBitsPerRow,IDBitsPerRow\n";
    method_cost_init.flush();
    std::ofstream sub_mode_init(output_dir + "/sub_mode_usage.csv");
    sub_mode_init << "Dataset,Algorithm,Mode\n";
    sub_mode_init.flush();
    std::ofstream sub_strategy_init(output_dir + "/sub_strategy_usage.csv");
    sub_strategy_init << "Dataset,Algorithm,Strategy\n";
    sub_strategy_init.flush();
    currentScenario = "subspace_g" + std::to_string(dimension_group_size);
    size_t start_idx = loggedResults.size();
    NoneCompressor none_comp;

    auto safe_run = [&](const std::string& method_name, const std::function<void()>& fn) {
        try {
            fn();
        } catch (const std::exception& e) {
            std::cerr << "Method failed [" << method_name << "]: " << e.what() << std::endl;
        }
    };

    for (size_t i = 0; i < datasets.size(); ++i) {
        if (dims[i] <= 1) continue;
        std::string path = file_root + datasets[i] + ".csv";
        int dim = dims[i];
        safe_run("KCluster", [&]() { run_cluster_encoder_test<KCluster>(datasets[i], dim, 0, path, output_dir, &none_comp); });
        safe_run("ACluster", [&]() { run_cluster_encoder_test<ACluster>(datasets[i], dim, 0, path, output_dir, &none_comp); });
        safe_run("KClusterSub", [&]() { run_cluster_encoder_test<KClusterSub>(datasets[i], dim, 0, path, output_dir, &none_comp, dimension_group_size); });
        safe_run("AClusterSub", [&]() { run_cluster_encoder_test<AClusterSub>(datasets[i], dim, 0, path, output_dir, &none_comp, dimension_group_size); });
    }

    std::ofstream per_dataset(output_dir + "/subspace_group_by_dataset.csv");
    per_dataset << "Dataset,Algorithm,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row\n";
    std::map<std::string, std::vector<LoggedResultRecord>> grouped;
    for (size_t i = start_idx; i < loggedResults.size(); ++i) {
        const auto& r = loggedResults[i];
        if (r.scenario != currentScenario) continue;
        per_dataset << r.dataset << "," << r.algorithm << ","
                    << std::fixed << std::setprecision(4) << r.compression_ratio << ","
                    << std::setprecision(2) << r.enc_avg_ns << "," << r.dec_avg_ns << "\n";
        grouped[r.algorithm].push_back(r);
    }
    per_dataset.flush();

    std::ofstream cross_avg(output_dir + "/subspace_group_cross_dataset_avg.csv");
    cross_avg << "Algorithm,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row\n";
    for (const auto& kv : grouped) {
        const auto& vec = kv.second;
        double sum_ratio = 0.0, sum_enc = 0.0, sum_dec = 0.0;
        for (const auto& r : vec) {
            sum_ratio += r.compression_ratio;
            sum_enc += r.enc_avg_ns;
            sum_dec += r.dec_avg_ns;
        }
        double n = vec.empty() ? 1.0 : static_cast<double>(vec.size());
        cross_avg << kv.first << ","
                  << std::fixed << std::setprecision(4) << (sum_ratio / n) << ","
                  << std::setprecision(2) << (sum_enc / n) << "," << (sum_dec / n) << "\n";
    }
    cross_avg.flush();

    std::ifstream sub_cost_in(output_dir + "/method_cost_components_by_dataset.csv");
    if (sub_cost_in.is_open()) {
        std::string line;
        std::getline(sub_cost_in, line);
        struct CostAgg { double h=0,b=0,r=0,i=0; };
        std::map<std::string, CostAgg> sums;
        std::map<std::string, int> counts;
        while (std::getline(sub_cost_in, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            std::string dataset, algo, h_s, b_s, r_s, i_s;
            if (!std::getline(ss, dataset, ',')) continue;
            if (!std::getline(ss, algo, ',')) continue;
            if (!std::getline(ss, h_s, ',')) continue;
            if (!std::getline(ss, b_s, ',')) continue;
            if (!std::getline(ss, r_s, ',')) continue;
            if (!std::getline(ss, i_s, ',')) continue;
            double h_v = 0.0, b_v = 0.0, r_v = 0.0, i_v = 0.0;
            try {
                h_v = std::stod(h_s);
                b_v = std::stod(b_s);
                r_v = std::stod(r_s);
                i_v = std::stod(i_s);
            } catch (...) {
                continue;
            }
            sums[algo].h += h_v;
            sums[algo].b += b_v;
            sums[algo].r += r_v;
            sums[algo].i += i_v;
            counts[algo] += 1;
        }
        std::ofstream sub_cost_avg(output_dir + "/method_cost_components_avg.csv");
        sub_cost_avg << "Algorithm,AvgHeaderBitsPerRow,AvgBaseBitsPerRow,AvgResidualBitsPerRow,AvgIDBitsPerRow\n";
        for (const auto& kv : sums) {
            int cnt = counts[kv.first] > 0 ? counts[kv.first] : 1;
            sub_cost_avg << kv.first << ","
                         << std::fixed << std::setprecision(4) << (kv.second.h / cnt) << ","
                         << std::fixed << std::setprecision(4) << (kv.second.b / cnt) << ","
                         << std::fixed << std::setprecision(4) << (kv.second.r / cnt) << ","
                         << std::fixed << std::setprecision(4) << (kv.second.i / cnt) << "\n";
        }
        sub_cost_avg.flush();
    }
}

void run_kcluster_kmeans_all_dataset_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    init_results_csv(output_dir + "/kcluster_kmeans_raw_results.csv");
    currentScenario = "kcluster_kmeans_all";
    size_t start_idx = loggedResults.size();
    NoneCompressor none_comp;

    auto safe_run = [&](const std::string& method_name, const std::function<void()>& fn) {
        try {
            fn();
        } catch (const std::exception& e) {
            std::cerr << "Method failed [" << method_name << "]: " << e.what() << std::endl;
        }
    };

    for (size_t i = 0; i < datasets.size(); ++i) {
        std::string path = file_root + datasets[i] + ".csv";
        int dim = dims[i];
        safe_run("KClusterKMeans", [&]() { run_cluster_encoder_test<KClusterKMeans>(datasets[i], dim, 0, path, output_dir, &none_comp); });
    }

    std::ofstream per_dataset(output_dir + "/kcluster_kmeans_by_dataset.csv");
    per_dataset << "Dataset,CompressionRatio,EncTime_ns_per_row,DecTime_ns_per_row\n";
    double sum_ratio = 0.0, sum_enc = 0.0, sum_dec = 0.0;
    int cnt = 0;
    for (size_t i = start_idx; i < loggedResults.size(); ++i) {
        const auto& r = loggedResults[i];
        if (r.scenario != currentScenario || r.algorithm != "KClusterKMeans") continue;
        per_dataset << r.dataset << ","
                    << std::fixed << std::setprecision(4) << r.compression_ratio << ","
                    << std::setprecision(2) << r.enc_avg_ns << "," << r.dec_avg_ns << "\n";
        sum_ratio += r.compression_ratio;
        sum_enc += r.enc_avg_ns;
        sum_dec += r.dec_avg_ns;
        cnt += 1;
    }
    per_dataset.flush();

    std::ofstream avg_out(output_dir + "/kcluster_kmeans_cross_dataset_avg.csv");
    avg_out << "Algorithm,AvgCompressionRatio,AvgEncTime_ns_per_row,AvgDecTime_ns_per_row\n";
    double n = cnt > 0 ? static_cast<double>(cnt) : 1.0;
    avg_out << "KClusterKMeans,"
            << std::fixed << std::setprecision(4) << (sum_ratio / n) << ","
            << std::setprecision(2) << (sum_enc / n) << "," << (sum_dec / n) << "\n";
    avg_out.flush();
}


int main() {
    try {
        std::cout << std::unitbuf;
        std::string fileRoot = "../data/";
        std::string outputDir = "output";
        std::filesystem::create_directories(outputDir);

        init_results_csv(outputDir + "/benchmark_results.csv");

        std::vector<std::string> datasets = {"SSD-bench", "profile-income", "Blockchain-tr", "Crop", "gas","Wind-Speed", "CT"};
        std::vector<int> data_size = {8926, 14825, 99999, 24000, 13910, 100000, 581012};
        
        std::vector<int> dims = {1, 1, 1, 46, 64, 1, 10};
        int dimension_group_size = 4;
        std::vector<double> reference_ratios = {0.10, 0.20, 0.30, 0.40, 0.50};
        std::vector<int> dim_expand_list = {8, 16, 32, 64};

        // =========================
        // Primary Benchmark Entry
        // =========================
        // Runs the default full benchmark loop and logs compression ratio / time for major methods.
        NoneCompressor none_comp;
        GzipCompressor gzip_comp;
        std::vector<SecondaryCompressor*> strategies = {&gzip_comp};
        // std::vector<SecondaryCompressor*> strategies = {&none_comp, &gzip_comp};

        for (size_t i = 0; i < datasets.size(); ++i) {
            std::string file_path = fileRoot + datasets[i] + ".csv";
            for (auto* strategy : strategies) {
                run_columnar_encoder_test<Chimp>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_columnar_encoder_test<Elf>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_columnar_encoder_test<Gorilla>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_columnar_encoder_test<Rle>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_columnar_encoder_test<Huffman>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_columnar_encoder_test<Alp>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_multidim_reger_test(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_cluster_encoder_test<KCluster>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
                run_cluster_encoder_test<ACluster>(datasets[i], dims[i], data_size[i], file_path, outputDir, strategy);
            }
        }

        // =========================
        // Experiment Call Catalog
        // =========================
        // Runs a single Vortex + RLE multidimensional test on one dataset.
        // run_vortex_test(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir);

        // Runs one scalar columnar compression baseline test (Elf/Gorilla/Chimp/RLE/Huffman/ALP).
        // run_columnar_encoder_test<Chimp>(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir, &gzip_comp);

        // Runs one multi-dimensional Reger test.
        // run_multidim_reger_test(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir, &gzip_comp);

        // Runs one cluster encoder test instance for KCluster with optional secondary compressor.
        // run_cluster_encoder_test<KCluster>(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir, &none_comp);

        // Runs one byte-oriented baseline compression test (Vortex+RLE over raw row-major bytes).
        // run_byte_compression_test(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir);

        // Runs one parameterized cluster experiment instance (k/page_size/pack_size).
        // run_cluster_encoder_experiment_instance<KCluster>(datasets[0], dims[0], data_size[0], fileRoot + datasets[0] + ".csv", outputDir, 100, 10000, 10);

        // Runs k-value sensitivity experiment.
        // run_k_sensitivity_experiment(fileRoot, outputDir);

        // Runs page-size sensitivity experiment.
        // run_page_size_sensitivity_experiment(fileRoot, outputDir);

        // Runs pack-size sensitivity experiment.
        // run_pack_size_sensitivity_experiment(fileRoot, outputDir);

        // Runs cluster reference-ratio sweep for KCluster and ACluster.
        // run_cluster_reference_experiment(fileRoot, outputDir, datasets, dims, reference_ratios);

        // Runs ACluster three-mode reference experiment (NoRef / AdjRef / RandomShuffle).
        // run_acluster_three_mode_reference_experiment(fileRoot, outputDir, datasets, dims, reference_ratios);

        // Runs ACluster outlier-reference experiment with static and dynamic modes.
        // run_acluster_outlier_experiment(fileRoot, outputDir, datasets, dims, reference_ratios);

        // Runs dimensionality impact experiment by expanding selected dimensions.
        // run_dimensionality_impact_experiment(fileRoot, outputDir, datasets.back(), dims.back(), dim_expand_list);

        // Runs adjacent-reference sensitivity experiment for ACluster.
        // run_acluster_adjacent_reference_experiment(fileRoot, outputDir, datasets, dims);

        // Runs subspace grouping experiment for KClusterSub and AClusterSub.
        // run_subspace_grouping_experiment(fileRoot, outputDir, datasets, dims, dimension_group_size);

        // Runs KClusterKMeans full-dataset experiment.
        // run_kcluster_kmeans_all_dataset_experiment(fileRoot, outputDir, datasets, dims);
        
    } catch (const std::exception& e) {
        std::cerr << "A critical error occurred in main: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}

       
