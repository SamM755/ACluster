#pragma once

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cmath>
#include <ctime>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "algorithms/alp.h"
#include "algorithms/chimp.h"
#include "algorithms/elf.h"
#include "algorithms/gorilla.h"
#include "algorithms/huffman.h"
#include "algorithms/rle.h"
#include "algorithms/reger.h"
#include "algorithms/vortex.h"
#include "algorithms/cluster_encoder_logic.h"
using EndToEndRows = std::vector<std::vector<double>>;
struct EndToEndSelection {
    std::vector<int> columns;
    std::vector<int> dimension_path;
    double best_cv_drop = 0;
};
struct EndToEndMetrics {
    uint64_t cluster_bytes = 0, alp_bytes = 0, correction_bytes = 0, metadata_bytes = 0;
    uint64_t baseline_bytes = 0;
    uint64_t correction_values = 0, references = 0, pages = 0;
    long long selection_ns = 0, preprocessing_ns = 0, cluster_ns = 0, alp_ns = 0;
    long long encode_ns = 0, decode_ns = 0;
    long long baseline_ns = 0;
};
struct EndToEndEncoding {
    std::vector<uint8_t> bytes;
    EndToEndSelection selection;
    EndToEndMetrics metrics;
};
EndToEndRows read_end_to_end_table(const std::string& path, size_t expected_rows, int dim);
EndToEndSelection select_end_to_end_columns(const EndToEndRows& rows);
EndToEndEncoding encode_end_to_end_table(const EndToEndRows& rows,
    const std::string& method, uint32_t seed = 20260929, int page_size = 10000);
EndToEndRows decode_end_to_end_table(const std::vector<uint8_t>& bytes);
void verify_end_to_end_table(const EndToEndRows& original, const EndToEndRows& decoded, int page_size = 10000);
void test_end_to_end_pipeline();

struct OptimalPageMeasurement {
    ClusterLogic::ClusterPageEncoding encoding;
    long long encode_ns = 0, decode_ns = 0, preprocess_ns = 0;
    long long codec_ns = 0, read_ns = 0, write_ns = 0;
};
void test_optimal_1d();
OptimalPageMeasurement measure_optimal_comparison_page(
    const std::vector<std::string>& page, const std::string& method, int k,
    int page_size, int pack_size, int block_size, uint32_t seed,
    const std::vector<std::vector<double>>& sorted_expected,
    const std::string& input_path = "", std::streamoff input_offset = 0,
    const std::string& output_path = "");

void write_binary_file(const std::string&, const std::vector<uint8_t>&);
