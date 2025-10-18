#include "reger.h"
#include "../utils/bit_stream.h"
#include "../utils/type_converter.h"
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <array>
#include <stdexcept>
#include <cstring>
#include <limits>

#ifdef _MSC_VER
#include <intrin.h>
#endif

namespace { 
namespace RegerUtils {

inline int getBitWith(uint64_t num) {
    if (num == 0) return 1;
#ifdef _MSC_VER
    unsigned long index; if (_BitScanReverse64(&index, num)) return static_cast<int>(index) + 1; return 1;
#else
    return 64 - __builtin_clzll(num);
#endif
}

inline uint64_t zigzag_encode(int64_t n) { return (static_cast<uint64_t>(n) << 1) ^ (n >> 63); }
inline int64_t zigzag_decode(uint64_t n) { return (n >> 1) ^ -static_cast<int64_t>(n & 1); }

constexpr int SEGMENT_SIZE = 16; 

struct ColumnCostInfo {
    int64_t total_bits = 0;
    std::array<double, 2> theta_v = {0.0, 1.0};
    int64_t min_residual = 0;
    std::vector<int> segment_bit_widths;
    std::vector<int64_t> residuals;
};

// OPTIMIZED: Operates on an indirect view of the data via p_map, for a single column.
ColumnCostInfo calculate_single_col_cost_on_view(
    const std::vector<std::vector<int64_t>>& block, 
    const std::vector<uint32_t>& p_map,
    int col_idx,
    bool store_residuals = false) 
{
    ColumnCostInfo result;
    if (p_map.size() < 2) { result.total_bits = p_map.empty() ? 0 : 64; return result; }
    
    long double sum_x=0, sum_y=0, sum_xx=0, sum_xy=0;
    for (size_t i = 1; i < p_map.size(); ++i) {
        long double x = block[p_map[i-1]][col_idx];
        long double y = block[p_map[i]][col_idx];
        sum_x+=x; sum_y+=y; sum_xx+=x*x; sum_xy+=x*y;
    }
    int m = p_map.size() - 1;
    double den = (double)m * sum_xx - sum_x * sum_x;
    if (std::abs(den) > 1e-9) { result.theta_v[0] = (sum_xx*sum_y-sum_x*sum_xy)/den; result.theta_v[1] = ((double)m*sum_xy-sum_x*sum_y)/den; }
    
    std::vector<int64_t> temp_residuals; temp_residuals.reserve(m);
    result.min_residual = std::numeric_limits<int64_t>::max();
    for (size_t i = 1; i < p_map.size(); ++i) {
        int64_t p = static_cast<int64_t>(std::round(result.theta_v[0] + result.theta_v[1] * block[p_map[i-1]][col_idx]));
        int64_t r = block[p_map[i]][col_idx] - p;
        temp_residuals.push_back(r);
        if (r < result.min_residual) result.min_residual = r;
    }
    if (p_map.size() <= 1) result.min_residual = 0;
    if (store_residuals) result.residuals = temp_residuals;

    int64_t total_res_bits = 0;
    int num_segments = (m + SEGMENT_SIZE - 1) / SEGMENT_SIZE;
    if (store_residuals) result.segment_bit_widths.reserve(num_segments);

    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        uint64_t max_zigzag_res_in_segment = 0;
        size_t start = seg_idx * SEGMENT_SIZE;
        size_t end = std::min(start + SEGMENT_SIZE, temp_residuals.size());
        for (size_t i = start; i < end; ++i) {
            uint64_t zz_res = zigzag_encode(temp_residuals[i] - result.min_residual);
            if (zz_res > max_zigzag_res_in_segment) max_zigzag_res_in_segment = zz_res;
        }
        int bit_width = getBitWith(max_zigzag_res_in_segment);
        if (store_residuals) result.segment_bit_widths.push_back(bit_width);
        total_res_bits += (end - start) * bit_width;
    }
    result.total_bits = total_res_bits + 64 + 64 + 64 + 64;
    return result;
}

// OPTIMIZED: This version now uses the _on_view helper.
int64_t calculate_multidim_cost_on_view(const std::vector<std::vector<int64_t>>& block, const std::vector<uint32_t>& p_map) {
    if (block.empty()) return 0;
    int64_t total_cost = 0;
    for (size_t j = 0; j < block[0].size(); ++j) {
        total_cost += calculate_single_col_cost_on_view(block, p_map, j).total_bits;
    }
    return total_cost;
}

int find_candidate_row_on_view(const std::vector<std::vector<int64_t>>& block, const std::vector<uint32_t>& p_map) {
    if (p_map.size() < 2) return -1;
    std::vector<double> row_error(p_map.size(), 0.0);
    
    for (size_t j = 0; j < block[0].size(); ++j) {
        auto info = calculate_single_col_cost_on_view(block, p_map, j);
        for (size_t i = 1; i < p_map.size(); ++i) {
            int64_t p = static_cast<int64_t>(std::round(info.theta_v[0] + info.theta_v[1] * block[p_map[i-1]][j]));
            row_error[i] += std::abs(static_cast<double>(block[p_map[i]][j]) - p);
        }
    }
    auto max_it = std::max_element(row_error.begin() + 1, row_error.end());
    return (max_it != row_error.end()) ? std::distance(row_error.begin(), max_it) : -1;
}

ColumnCostInfo calculate_single_col_cost(const std::vector<int64_t>& col_data, bool store_residuals = false) {
    ColumnCostInfo result;
    if (col_data.size() < 2) { result.total_bits = col_data.empty() ? 0 : 64; return result; }
    
    long double sum_x=0, sum_y=0, sum_xx=0, sum_xy=0;
    for (size_t i = 1; i < col_data.size(); ++i) { long double x = col_data[i-1], y = col_data[i]; sum_x+=x; sum_y+=y; sum_xx+=x*x; sum_xy+=x*y; }
    int m = col_data.size() - 1;
    double den = (double)m * sum_xx - sum_x * sum_x;
    if (std::abs(den) > 1e-9) { result.theta_v[0] = (sum_xx*sum_y-sum_x*sum_xy)/den; result.theta_v[1] = ((double)m*sum_xy-sum_x*sum_y)/den; }
    
    std::vector<int64_t> temp_residuals; temp_residuals.reserve(m);
    result.min_residual = std::numeric_limits<int64_t>::max();
    for (size_t i = 1; i < col_data.size(); ++i) {
        int64_t p = static_cast<int64_t>(std::round(result.theta_v[0] + result.theta_v[1] * col_data[i-1]));
        int64_t r = col_data[i] - p;
        temp_residuals.push_back(r);
        if (r < result.min_residual) result.min_residual = r;
    }
    if (col_data.size() <= 1) result.min_residual = 0;
    if (store_residuals) result.residuals = temp_residuals;

    int64_t total_res_bits = 0;
    int num_segments = (m + SEGMENT_SIZE - 1) / SEGMENT_SIZE;
    if (store_residuals) result.segment_bit_widths.reserve(num_segments);

    for (int seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
        uint64_t max_zigzag_res_in_segment = 0;
        size_t start = seg_idx * SEGMENT_SIZE;
        size_t end = std::min(start + SEGMENT_SIZE, temp_residuals.size());
        for (size_t i = start; i < end; ++i) {
            uint64_t zz_res = zigzag_encode(temp_residuals[i] - result.min_residual);
            if (zz_res > max_zigzag_res_in_segment) max_zigzag_res_in_segment = zz_res;
        }
        int bit_width = getBitWith(max_zigzag_res_in_segment);
        if (store_residuals) result.segment_bit_widths.push_back(bit_width);
        total_res_bits += (end - start) * bit_width;
    }

    result.total_bits = total_res_bits + 64 + 64 + 64 + 64; // metadata cost
    return result;
}

int64_t calculate_multidim_cost(const std::vector<std::vector<int64_t>>& block) {
    if (block.empty()) return 0;
    int n_rows = block.size(), n_cols = block[0].size();
    int64_t total_cost = 0;
    std::vector<int64_t> col_data(n_rows);
    for (int j = 0; j < n_cols; ++j) {
        for (int i = 0; i < n_rows; ++i) col_data[i] = block[i][j];
        total_cost += calculate_single_col_cost(col_data).total_bits;
    }
    return total_cost;
}

int find_candidate_row(const std::vector<std::vector<int64_t>>& block) {
    if (block.size() < 2) return -1;
    int n_cols = block[0].size();
    std::vector<double> row_error(block.size(), 0.0);
    std::vector<int64_t> col_data(block.size());
    for (int j = 0; j < n_cols; ++j) {
        for (size_t i = 0; i < block.size(); ++i) col_data[i] = block[i][j];
        auto info = calculate_single_col_cost(col_data);
        for (size_t i = 1; i < block.size(); ++i) {
            int64_t p = static_cast<int64_t>(std::round(info.theta_v[0] + info.theta_v[1] * col_data[i-1]));
            row_error[i] += std::abs(static_cast<double>(col_data[i]) - p);
        }
    }
    auto max_it = std::max_element(row_error.begin() + 1, row_error.end());
    return (max_it != row_error.end()) ? std::distance(row_error.begin(), max_it) : -1;
}

// original
// void ReorderingTimeSeries_MultiDim(std::vector<std::vector<int64_t>>& block, bool is_first_block) {
//     if (block.size() < 2) return;
    
//     // Phase 1: Find best initial permutation
//     std::vector<std::vector<int64_t>> best_initial_block = block;
//     int64_t best_initial_cost = calculate_multidim_cost(best_initial_block);
    
//     for (int j = 0; j < block[0].size(); ++j) {
//         std::vector<std::vector<int64_t>> temp_block = block;
//         std::sort(temp_block.begin(), temp_block.end(), [j](const auto& a, const auto& b){ return a[j] < b[j]; });
//         int64_t cost = calculate_multidim_cost(temp_block);
//         if (cost < best_initial_cost) {
//             best_initial_cost = cost;
//             best_initial_block = temp_block;
//         }
//     }
//     block = best_initial_block;

//     // Phase 2: Iterative Refinement
//     int64_t current_best_cost = best_initial_cost;
//     int max_iter = 5;
//     for (int iter = 0; iter < max_iter; ++iter) {
//         int alpha = find_candidate_row(block);
//         if (alpha == -1) break;
//         int best_beta = -1;
//         int64_t cost_after_swap = current_best_cost;
//         for (size_t beta = 0; beta < block.size(); ++beta) {
//             if (beta == alpha) continue;
//             auto temp_block = block;
//             std::swap(temp_block[alpha], temp_block[beta]);
//             int64_t new_cost = calculate_multidim_cost(temp_block);
//             if (new_cost < cost_after_swap) { cost_after_swap = new_cost; best_beta = beta; }
//         }
//         if (best_beta != -1) {
//             std::swap(block[alpha], block[best_beta]);
//             current_best_cost = cost_after_swap;
//         } else {
//             break;
//         }
//     }
// }

// OPTIMIZED: The main reordering logic is heavily optimized to avoid data copies.
void ReorderingTimeSeries_MultiDim(std::vector<std::vector<int64_t>>& block, bool is_first_block) {
    if (block.size() < 2) return;
    
    std::vector<uint32_t> p_map(block.size());
    std::iota(p_map.begin(), p_map.end(), 0);
    
    int64_t best_initial_cost = calculate_multidim_cost_on_view(block, p_map);
    std::vector<uint32_t> best_initial_pmap = p_map;

    for (size_t j = 0; j < block[0].size(); ++j) {
        std::vector<uint32_t> temp_pmap = p_map;
        std::sort(temp_pmap.begin(), temp_pmap.end(), 
            [&](uint32_t a, uint32_t b){ return block[a][j] < block[b][j]; });
        
        int64_t cost = calculate_multidim_cost_on_view(block, temp_pmap);
        if (cost < best_initial_cost) {
            best_initial_cost = cost;
            best_initial_pmap = temp_pmap;
        }
    }
    p_map = best_initial_pmap;
    
    int64_t current_best_cost = best_initial_cost;
    int max_iter = 3;
    for (int iter = 0; iter < max_iter; ++iter) {
        int alpha_idx = find_candidate_row_on_view(block, p_map);
        if (alpha_idx == -1) break;
        
        int best_beta_idx = -1;
        int64_t cost_after_swap = current_best_cost;
        
        size_t search_range = 16;
        size_t start_beta = (alpha_idx > (int)search_range) ? alpha_idx - search_range : 0;
        size_t end_beta = std::min(p_map.size(), (size_t)alpha_idx + search_range);

        for (size_t beta_idx = start_beta; beta_idx < end_beta; ++beta_idx) {
            if (beta_idx == (size_t)alpha_idx) continue;
            std::vector<uint32_t> temp_pmap = p_map;
            std::swap(temp_pmap[alpha_idx], temp_pmap[beta_idx]);
            int64_t new_cost = calculate_multidim_cost_on_view(block, temp_pmap);
            if (new_cost < cost_after_swap) {
                cost_after_swap = new_cost;
                best_beta_idx = beta_idx;
            }
        }
        if (best_beta_idx != -1) {
            std::swap(p_map[alpha_idx], p_map[best_beta_idx]);
            current_best_cost = cost_after_swap;
        } else {
            break;
        }
    }
    
    std::vector<std::vector<int64_t>> final_reordered_block(block.size());
    for(size_t i=0; i<block.size(); ++i) final_reordered_block[i] = block[p_map[i]];
    block = final_reordered_block;
}

} // namespace RegerUtils
}


// original
// std::vector<uint8_t> Reger::encode_multidim(const std::vector<std::vector<int64_t>>& rows, const std::vector<double>& scaling_factors) {
//     if (rows.empty()) return {};
//     const int block_size = 512;
//     OutputBitStream out;
//     uint32_t num_cols = rows[0].size();
//     out.write_bits(rows.size(), 64);
//     out.write_bits(num_cols, 32);
//     out.write_bits(block_size, 32);
//     out.write_bits(static_cast<uint16_t>(RegerUtils::SEGMENT_SIZE), 16);
//     for (double factor : scaling_factors) out.write_bits(double_to_long(factor), 64);
    
//     bool first_block_processed = false;
//     for (size_t i = 0; i < rows.size(); i += block_size) {
//         bool is_first_block = !first_block_processed;
//         size_t end = std::min(i + (size_t)block_size, rows.size());
//         std::vector<std::vector<int64_t>> block;
//         for(size_t j = i; j < end; ++j) block.push_back(rows[j]);
        
//         RegerUtils::ReorderingTimeSeries_MultiDim(block, is_first_block);

//         out.write_bits(static_cast<uint32_t>(block.size()), 32);
        
//         for (size_t j = 0; j < num_cols; ++j) {
//             std::vector<int64_t> col_data;
//             for(const auto& row : block) col_data.push_back(row[j]);
//             auto info = RegerUtils::calculate_single_col_cost(col_data, true);
            
//             out.write_bits(double_to_long(info.theta_v[0]), 64);
//             out.write_bits(double_to_long(info.theta_v[1]), 64);
//             out.write_bits(static_cast<uint64_t>(col_data.empty() ? 0 : col_data[0]), 64);
//             out.write_bits(static_cast<uint64_t>(info.min_residual), 64);
//             out.write_bits(static_cast<uint16_t>(info.segment_bit_widths.size()), 16);
//             for(int width : info.segment_bit_widths) out.write_bits(width, 8);
            
//             for(size_t seg_idx = 0; seg_idx < info.segment_bit_widths.size(); ++seg_idx) {
//                 int bit_width = info.segment_bit_widths[seg_idx];
//                 size_t res_start = seg_idx * RegerUtils::SEGMENT_SIZE;
//                 size_t res_end = std::min(res_start + RegerUtils::SEGMENT_SIZE, info.residuals.size());
//                 if (bit_width > 0) {
//                     for(size_t res_idx = res_start; res_idx < res_end; ++res_idx) {
//                         out.write_bits(RegerUtils::zigzag_encode(info.residuals[res_idx] - info.min_residual), bit_width);
//                     }
//                 }
//             }
//         }
//         first_block_processed = true;
//     }
//     return out.get_bytes();
// }

// std::pair<std::vector<std::vector<int64_t>>, std::vector<double>> Reger::decode_multidim(const std::vector<uint8_t>& compressed_data) {
//     if (compressed_data.empty()) return {};
//     InputBitStream in(compressed_data);
//     uint64_t total_rows = in.read_bits(64);
//     uint32_t num_cols = in.read_bits(32);
//     in.read_bits(32); // block_size_meta
//     uint16_t segment_size = in.read_bits(16);
//     std::vector<double> scaling_factors(num_cols);
//     for (uint32_t j = 0; j < num_cols; ++j) scaling_factors[j] = long_to_double(in.read_bits(64));
//     std::vector<std::vector<int64_t>> decoded_rows;
//     decoded_rows.reserve(total_rows);
    
//     uint64_t rows_decoded_count = 0;
//     while(rows_decoded_count < total_rows) {
//         uint32_t current_block_size = in.read_bits(32);
//         if (current_block_size == 0) break;
//         std::vector<std::vector<int64_t>> decoded_block(current_block_size, std::vector<int64_t>(num_cols));

//         for (uint32_t j = 0; j < num_cols; ++j) {
//             double theta0 = long_to_double(in.read_bits(64));
//             double theta1 = long_to_double(in.read_bits(64));
//             int64_t prev_val = static_cast<int64_t>(in.read_bits(64));
//             int64_t min_res = static_cast<int64_t>(in.read_bits(64));
//             if (current_block_size > 0) decoded_block[0][j] = prev_val;

//             uint16_t num_segments = in.read_bits(16);
//             std::vector<uint8_t> segment_bit_widths(num_segments);
//             for(auto& width : segment_bit_widths) width = in.read_bits(8);
            
//             size_t res_decoded_count = 0;
//             for(uint16_t seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
//                 uint8_t bit_width = segment_bit_widths[seg_idx];
//                 size_t res_to_decode_in_seg = std::min((size_t)segment_size, (current_block_size > 0 ? current_block_size - 1 : 0) - res_decoded_count);
//                 for(size_t k = 0; k < res_to_decode_in_seg; ++k) {
//                     int64_t pred = static_cast<int64_t>(std::round(theta0 + theta1 * prev_val));
//                     uint64_t packed_res = (bit_width > 0) ? in.read_bits(bit_width) : 0;
//                     int64_t unzigzagged_res = RegerUtils::zigzag_decode(packed_res);
//                     int64_t current_val = pred + unzigzagged_res + min_res;
//                     if (res_decoded_count + k + 1 < current_block_size) {
//                         decoded_block[res_decoded_count + k + 1][j] = current_val;
//                     }
//                     prev_val = current_val;
//                 }
//                 res_decoded_count += res_to_decode_in_seg;
//             }
//         }
//         for (const auto& row : decoded_block) {
//             decoded_rows.push_back(row);
//         }
//         rows_decoded_count += current_block_size;
//     }
//     return {decoded_rows, scaling_factors};
// }

// optimized
std::vector<uint8_t> Reger::encode_multidim(const std::vector<std::vector<int64_t>>& rows, const std::vector<double>& scaling_factors) {
    if (rows.empty()) return {};
    const int block_size = 512;
    OutputBitStream out;
    uint32_t num_cols = rows[0].size();
    out.write_bits(rows.size(), 64);
    out.write_bits(num_cols, 32);
    out.write_bits(block_size, 32);
    out.write_bits(static_cast<uint16_t>(RegerUtils::SEGMENT_SIZE), 16);
    for (double factor : scaling_factors) out.write_bits(double_to_long(factor), 64);
    
    for (size_t i = 0; i < rows.size(); i += block_size) {
        size_t end = std::min(i + (size_t)block_size, rows.size());
        std::vector<std::vector<int64_t>> block;
        for(size_t j = i; j < end; ++j) block.push_back(rows[j]);
        
        RegerUtils::ReorderingTimeSeries_MultiDim(block, i==0);
        out.write_bits(static_cast<uint32_t>(block.size()), 32);
        
        for (size_t j = 0; j < num_cols; ++j) {
            std::vector<int64_t> col_data;
            for(const auto& row : block) col_data.push_back(row[j]);
            auto info = RegerUtils::calculate_single_col_cost(col_data, true);
            out.write_bits(double_to_long(info.theta_v[0]), 64);
            out.write_bits(double_to_long(info.theta_v[1]), 64);
            out.write_bits(static_cast<uint64_t>(col_data.empty() ? 0 : col_data[0]), 64);
            out.write_bits(static_cast<uint64_t>(info.min_residual), 64);
            out.write_bits(static_cast<uint16_t>(info.segment_bit_widths.size()), 16);
            for(int width : info.segment_bit_widths) out.write_bits(width, 8);
            for(size_t seg_idx = 0; seg_idx < info.segment_bit_widths.size(); ++seg_idx) {
                int bit_width = info.segment_bit_widths[seg_idx];
                size_t res_start = seg_idx * RegerUtils::SEGMENT_SIZE;
                size_t res_end = std::min(res_start + RegerUtils::SEGMENT_SIZE, info.residuals.size());
                if (bit_width > 0) {
                    for(size_t res_idx = res_start; res_idx < res_end; ++res_idx) {
                        out.write_bits(RegerUtils::zigzag_encode(info.residuals[res_idx] - info.min_residual), bit_width);
                    }
                }
            }
        }
    }
    return out.get_bytes();
}

std::pair<std::vector<std::vector<int64_t>>, std::vector<double>> Reger::decode_multidim(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    InputBitStream in(compressed_data);
    uint64_t total_rows = in.read_bits(64);
    uint32_t num_cols = in.read_bits(32);
    in.read_bits(32); // block_size_meta
    uint16_t segment_size = in.read_bits(16);
    std::vector<double> scaling_factors(num_cols);
    for (uint32_t j = 0; j < num_cols; ++j) scaling_factors[j] = long_to_double(in.read_bits(64));
    std::vector<std::vector<int64_t>> decoded_rows;
    decoded_rows.reserve(total_rows);
    
    uint64_t rows_decoded_count = 0;
    while(rows_decoded_count < total_rows) {
        uint32_t current_block_size = in.read_bits(32);
        if (current_block_size == 0) break;
        std::vector<std::vector<int64_t>> decoded_block(current_block_size, std::vector<int64_t>(num_cols));
        for (uint32_t j = 0; j < num_cols; ++j) {
            double theta0 = long_to_double(in.read_bits(64));
            double theta1 = long_to_double(in.read_bits(64));
            int64_t prev_val = static_cast<int64_t>(in.read_bits(64));
            int64_t min_res = static_cast<int64_t>(in.read_bits(64));
            if (current_block_size > 0) decoded_block[0][j] = prev_val;
            uint16_t num_segments = in.read_bits(16);
            std::vector<uint8_t> segment_bit_widths(num_segments);
            for(auto& width : segment_bit_widths) width = in.read_bits(8);
            size_t res_decoded_count = 0;
            for(uint16_t seg_idx = 0; seg_idx < num_segments; ++seg_idx) {
                uint8_t bit_width = segment_bit_widths[seg_idx];
                size_t res_to_decode_in_seg = std::min((size_t)segment_size, (current_block_size > 0 ? current_block_size - 1 : 0) - res_decoded_count);
                for(size_t k = 0; k < res_to_decode_in_seg; ++k) {
                    int64_t pred = static_cast<int64_t>(std::round(theta0 + theta1 * prev_val));
                    uint64_t packed_res = (bit_width > 0) ? in.read_bits(bit_width) : 0;
                    int64_t unzigzagged_res = RegerUtils::zigzag_decode(packed_res);
                    int64_t current_val = pred + unzigzagged_res + min_res;
                    if (res_decoded_count + k + 1 < current_block_size) decoded_block[res_decoded_count + k + 1][j] = current_val;
                    prev_val = current_val;
                }
                res_decoded_count += res_to_decode_in_seg;
            }
        }
        for (const auto& row : decoded_block) decoded_rows.push_back(row);
        rows_decoded_count += current_block_size;
    }
    return {decoded_rows, scaling_factors};
}