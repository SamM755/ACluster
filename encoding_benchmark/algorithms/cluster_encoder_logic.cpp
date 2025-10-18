#include "cluster_common.h"
#include "cluster_encoder_logic.h"
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdint>
#include <string>
#include <sstream>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <unordered_set>
#include <map>
#include <iostream> // For debug
#include <iomanip>  // For debug
#include <stdexcept>
#include <fstream>


namespace {
namespace GlobalDebug {
    bool DEBUG_PRINT_ENABLED = false;

}
}

namespace ClusterLogic {

struct BitBuffer {
    std::vector<uint8_t> buffer;
    size_t currentBitPos = 0;
    
    void appendBit(bool bit) {
        size_t byteIndex = currentBitPos / 8;
        int bitIndexInByte = currentBitPos % 8;
        if (byteIndex >= buffer.size()) buffer.push_back(0);
        if (bit) buffer[byteIndex] |= (1 << (7 - bitIndexInByte));
        currentBitPos++;
    }

    void merge(const BitBuffer& other) {
        for (size_t i = 0; i < other.currentBitPos; ++i) {
            size_t byteIndex = i / 8;
            int bitOffset = 7 - (i % 8);
            bool bit = (other.buffer[byteIndex] & (1 << bitOffset)) != 0;
            appendBit(bit);
        }
    }
    
    std::vector<uint8_t> toByteArray() const { return buffer; }
    size_t size() const { return currentBitPos; }
};

namespace BitBufferUtils {
    void appendToBitstream(BitBuffer& bitBuffer, long long num, int bits) {
        for (int i = bits - 1; i >= 0; --i) bitBuffer.appendBit((num >> i) & 1);
    }
    int bitsRequiredNoSign(long long value) {
        if (value == 0) return 1;
        unsigned long long u_value = (value > 0) ? value : -value;
        #if defined(__GNUC__) || defined(__clang__)
        return 64 - __builtin_clzll(u_value);
        #else
        // Fallback for MSVC or other compilers
        int length = 0; while(u_value > 0){ u_value >>= 1; length++; } return length == 0 ? 1 : length;
        #endif
    }
}


BitStreamReader::BitStreamReader(const std::vector<uint8_t>& data)
    : buffer(data), bufferSize(data.size()), currentBitPos(0) {}

long long BitStreamReader::readBits(int numBits) {
    if (numBits < 0 || numBits > 64) throw std::invalid_argument("Cannot read " + std::to_string(numBits) + " bits.");
    if (currentBitPos + numBits > bufferSize * 8) throw std::out_of_range("Not enough bits in stream to read " + std::to_string(numBits) + " bits.");
    
    long long value = 0;
    for (int i = 0; i < numBits; ++i) {
        size_t byteIndex = currentBitPos / 8;
        int bitIndex = 7 - (currentBitPos % 8);
        bool bit = (buffer[byteIndex] & (1 << bitIndex)) != 0;
        value = (value << 1) | (bit ? 1LL : 0LL);
        currentBitPos++;
    }
    return value;
}


namespace { // Anonymous helpers for Preprocessor
    void parse_csv_line(const std::string& line, std::vector<std::string>& cells) {
        std::stringstream line_stream(line);
        std::string cell;
        cells.clear();
        while (std::getline(line_stream, cell, ',')) cells.push_back(cell);
    }
    
    long long stringToLong(const std::string& s, int decimal_places_to_scale) {
        if (s.empty()) {
            return 0;
        }

        if (s.find('e') != std::string::npos || s.find('E') != std::string::npos) {
            return 0;
        }

        try {
            bool is_negative = (s[0] == '-');
            std::string abs_s = is_negative ? s.substr(1) : s;

            auto dot_pos = abs_s.find('.');
            long long int_val = 0;
            long long frac_val = 0;

            if (dot_pos == std::string::npos) {
                if (!abs_s.empty()) {
                int_val = std::stoll(abs_s);
                }
            } else {
                std::string integral_part_str = abs_s.substr(0, dot_pos);
                if (!integral_part_str.empty()) {
                    int_val = std::stoll(integral_part_str);
                }

                std::string fractional_part_str = abs_s.substr(dot_pos + 1);
                if (!fractional_part_str.empty()) {
                    if (fractional_part_str.length() > 18) {
                        fractional_part_str = fractional_part_str.substr(0, 18);
                    }
                    frac_val = std::stoll(fractional_part_str);
                    
                    int frac_len = fractional_part_str.length();
                    for (int i = 0; i < decimal_places_to_scale - frac_len; ++i) {
                        frac_val *= 10;
                    }
                    for (int i = 0; i < frac_len - decimal_places_to_scale; ++i) {
                        frac_val /= 10;
                    }
                }
            }
            
            long long total_abs_val = int_val;
            for (int i = 0; i < decimal_places_to_scale; ++i) {
                if (total_abs_val > std::numeric_limits<long long>::max() / 10) {
                    total_abs_val = std::numeric_limits<long long>::max(); 
                    break;
                }
                total_abs_val *= 10;
            }
            
            if (total_abs_val <= std::numeric_limits<long long>::max() - frac_val) {
                total_abs_val += frac_val;
            } else {
                total_abs_val = std::numeric_limits<long long>::max(); 
            }
            return is_negative ? -total_abs_val : total_abs_val;
        } catch (const std::exception& e) {
            return 0;
        }
    }
}
PreprocessorResult preprocessPageFromStrings(const std::vector<std::string>& csv_page, int& dim_ref) {
    if (csv_page.empty()) return {nullptr, {}, {}};
    std::vector<int> max_decimal_places;
    std::vector<std::vector<std::string>> page_as_strings;
    page_as_strings.reserve(csv_page.size());
    int num_rows = 0;
    bool is_dim_known = (dim_ref > 0);
    if (is_dim_known) max_decimal_places.assign(dim_ref, 0);
    for (const auto& line : csv_page) {
        if (line.empty()) continue;
        std::vector<std::string> cells;
        parse_csv_line(line, cells);
        if (!is_dim_known) {
            dim_ref = cells.size();
            if (dim_ref == 0) continue;
            max_decimal_places.assign(dim_ref, 0);
            is_dim_known = true;
        }
        if (cells.size() != static_cast<size_t>(dim_ref)) continue;
        for (int j = 0; j < dim_ref; ++j) {
            const auto& cell = cells[j];
            auto dot_pos = cell.find('.');
            if (dot_pos != std::string::npos) max_decimal_places[j] = std::max(max_decimal_places[j], (int)(cell.length() - dot_pos - 1));
        }
        page_as_strings.push_back(cells);
        num_rows++;
    }
    if (num_rows == 0) return {nullptr, {}, {}};
    auto integerSoA = std::make_unique<SoAData<long long>>();
    integerSoA->dim = dim_ref;
    integerSoA->num_rows = num_rows;
    integerSoA->columns.resize(dim_ref);
    std::vector<long long> minValues(dim_ref);
    for (int j = 0; j < dim_ref; ++j) {
        auto& column = integerSoA->columns[j];
        column.resize(num_rows);
        long long current_min = std::numeric_limits<long long>::max();
        for (int i = 0; i < num_rows; ++i) {
            long long val = stringToLong(page_as_strings[i][j], max_decimal_places[j]);
            column[i] = val;
            current_min = std::min(current_min, val);
        }
        minValues[j] = current_min;
    }
    for (int j = 0; j < dim_ref; ++j) {
        auto& column = integerSoA->columns[j];
        long long min_val = minValues[j];
        for (int i = 0; i < num_rows; ++i) column[i] -= min_val;
    }
    return {std::move(integerSoA), max_decimal_places, minValues};
}


namespace{
    inline int bitLength(long long value) {
        if (value == 0) return 1;
        unsigned long long u_value = (value > 0) ? value : -value;
        #if defined(__GNUC__) || defined(__clang__)
        return 64 - __builtin_clzll(u_value);
        #else
        int length = 0; while (u_value > 0) { u_value >>= 1; length++; } return length;
        #endif
    }

    long long calculateTotalLogCost(const std::vector<long long>& p1, const std::vector<long long>& p2, int dim) {
        long long totalCost = 0;
        for (int i = 0; i < dim; ++i) {
            long long residual = p1[i] - p2[i];
            totalCost += bitLength(residual >= 0 ? residual : -residual) + 1;
        }
        return totalCost;
    }
    long long calculateBasePointStorageCost(const std::vector<long long>& basePoint, int dim) {
        long long storageCost = 0;
        for (int i = 0; i < dim; ++i) {
            storageCost += bitLength(basePoint[i] >= 0 ? basePoint[i] : -basePoint[i]) + 1;
        }
        return storageCost;
    }
    long long manhattanDistance(const std::vector<long long>& p1, const std::vector<long long>& p2) {
        long long sum = 0;
        for (size_t i = 0; i < p1.size(); ++i) { sum += std::abs(p1[i] - p2[i]); }
        return sum;
    }

    // Helper to get a point (row) from SoA data as an AoS vector.
    std::vector<long long> get_point_from_soa(const SoAData<long long>& soa_data, int row_idx) {
        std::vector<long long> point(soa_data.dim);
        for (int d = 0; d < soa_data.dim; ++d) {
            point[d] = soa_data.columns[d][row_idx];
        }
        return point;
    }

    struct VectorHasher {
        std::size_t operator()(const std::vector<long long>& vec) const {
            std::size_t seed = vec.size();
            for(long long i : vec) { seed ^= std::hash<long long>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2); }
            return seed;
        }
    };

    std::vector<std::vector<long long>> acceleratedInitialization(const SoAData<long long>& data, int k, std::mt19937& gen) {
        std::vector<std::vector<long long>> medoids;
        if (data.num_rows == 0 || k <= 0) return medoids;
        medoids.reserve(k);
        
        // Use a set of vectors to ensure VALUE uniqueness, just like Java's Set<String>
        std::unordered_set<std::vector<long long>, VectorHasher> selected_medoids_set;

        std::uniform_int_distribution<> distrib(0, data.num_rows - 1);
        
        // Find the first medoid
        int first_index = distrib(gen);
        auto first_medoid = get_point_from_soa(data, first_index);
        medoids.push_back(first_medoid);
        selected_medoids_set.insert(first_medoid);
        
        std::vector<long long> distances(data.num_rows);
        for (int i = 0; i < data.num_rows; ++i) {
            distances[i] = manhattanDistance(get_point_from_soa(data, i), medoids[0]);
        }

        for (int i = 1; i < k; ++i) {
            if (selected_medoids_set.size() >= (size_t)data.num_rows) break;
            
            long long total_distance = std::accumulate(distances.begin(), distances.end(), 0LL);
            if (total_distance == 0) {
                int next_idx = 0;
                while(selected_medoids_set.count(get_point_from_soa(data, next_idx))) {
                    next_idx = (next_idx + 1) % data.num_rows;
                }
                auto new_medoid = get_point_from_soa(data, next_idx);
                medoids.push_back(new_medoid);
                selected_medoids_set.insert(new_medoid);
                continue;
            }
            
            // Weighted random selection
            std::uniform_real_distribution<double> dist_double(0.0, (double)total_distance);
            long long rand_val = static_cast<long long>(dist_double(gen));
            
            long long current_sum = 0;
            int chosen_idx = -1;
            for(int j = 0; j < data.num_rows; ++j) {
                current_sum += distances[j];
                if (current_sum >= rand_val) {
                    chosen_idx = j;
                    break;
                }
            }
            if (chosen_idx == -1) chosen_idx = data.num_rows - 1;

            // Linearly probe for the next UNSELECTED VALUE, same as Java
            while(selected_medoids_set.count(get_point_from_soa(data, chosen_idx))) {
                chosen_idx = (chosen_idx + 1) % data.num_rows;
            }
            
            auto new_medoid = get_point_from_soa(data, chosen_idx);
            medoids.push_back(new_medoid);
            selected_medoids_set.insert(new_medoid);

            // Update distances
            for (int idx = 0; idx < data.num_rows; ++idx) {
                long long dist_new = manhattanDistance(get_point_from_soa(data, idx), new_medoid);
                distances[idx] = std::min(distances[idx], dist_new);
            }
        }
        return medoids;
    }

    std::vector<std::vector<long long>> updateMedoids(const SoAData<long long>& data, const std::vector<int>& clusters, const std::vector<std::vector<long long>>& currentMedoids, int dim) {
        int k = currentMedoids.size();
        if (k == 0) return {};
        std::vector<std::vector<long long>> newMedoids(k, std::vector<long long>(dim));
        std::vector<std::vector<int>> cluster_point_indices(k);
        for (int i = 0; i < data.num_rows; ++i) {
            if (clusters[i] != -1 && clusters[i] < k) { cluster_point_indices[clusters[i]].push_back(i); }
        }
        for (int medoidIndex = 0; medoidIndex < k; ++medoidIndex) {
            const auto& member_indices = cluster_point_indices[medoidIndex];
            if (member_indices.empty()) { newMedoids[medoidIndex] = currentMedoids[medoidIndex]; continue; }
            
            long long best_total_cost = std::numeric_limits<long long>::max();
            int best_medoid_candidate_idx_in_data = -1;

            for (int candidate_idx_in_data : member_indices) {
                long long current_total_cost = 0;
                auto candidate_medoid = get_point_from_soa(data, candidate_idx_in_data);
                
                for (int member_idx_in_data : member_indices) {
                    if (candidate_idx_in_data == member_idx_in_data) continue;
                    auto other_point = get_point_from_soa(data, member_idx_in_data);
                    for(int d = 0; d < dim; ++d) {
                        long long residual = candidate_medoid[d] - other_point[d];
                        current_total_cost += bitLength(residual >= 0 ? residual : -residual) + 1;
                    }
                }

                // for (int member_idx_in_data : member_indices) {
                //     current_total_cost += calculateTotalLogCost(candidate_medoid, get_point_from_soa(data, member_idx_in_data), dim);
                // }

                if (current_total_cost < best_total_cost) {
                    best_total_cost = current_total_cost;
                    best_medoid_candidate_idx_in_data = candidate_idx_in_data;
                }
            }
            if (best_medoid_candidate_idx_in_data != -1) {
                newMedoids[medoidIndex] = get_point_from_soa(data, best_medoid_candidate_idx_in_data);
            } else { newMedoids[medoidIndex] = currentMedoids[medoidIndex]; }
        }
        return newMedoids;
    }

    struct MedSortHelper {
        std::vector<long long> medoid;
        long long size;
        int original_index;
        bool operator<(const MedSortHelper& other) const { return size < other.size; }
    };

    ClusteringResult sortResults(std::vector<std::vector<long long>>& medoids, std::vector<int>& assignment, std::vector<long long>& sizes) {
        int k = medoids.size();
        if (k == 0) return {{}, {}, {}};
        std::vector<MedSortHelper> sorters;
        sorters.reserve(k);
        for(int i = 0; i < k; ++i) { sorters.push_back({medoids[i], sizes[i], i}); }
        std::sort(sorters.begin(), sorters.end());
        ClusteringResult res;
        res.medoids.resize(k);
        res.cluster_sizes.resize(k);
        std::vector<int> old_to_new_map(k);
        for(int i = 0; i < k; ++i) {
            res.medoids[i] = sorters[i].medoid;
            res.cluster_sizes[i] = sorters[i].size;
            old_to_new_map[sorters[i].original_index] = i;
        }
        res.cluster_assignment.resize(assignment.size());
        for(size_t i = 0; i < assignment.size(); ++i) {
            if (assignment[i] != -1 && static_cast<size_t>(assignment[i]) < old_to_new_map.size()) {
                res.cluster_assignment[i] = old_to_new_map[assignment[i]];
            } else { res.cluster_assignment[i] = -1; }
        }
        return res;
    }
    std::vector<std::vector<long long>> convert_SoA_to_AoS(const SoAData<long long>& soa_data) {
    if (soa_data.num_rows == 0) return {};
    std::vector<std::vector<long long>> aos_data(soa_data.num_rows, std::vector<long long>(soa_data.dim));
    for (int j = 0; j < soa_data.dim; ++j) {
        for (int i = 0; i < soa_data.num_rows; ++i) {
            aos_data[i][j] = soa_data.columns[j][i];
        }
    }
    return aos_data;
}
}

// ablation study version
KMedoidResult kMedoidLogCost_with_diag(const SoAData<long long>& data, int k, int max_iter, double tol, int dim) {
    int n = data.num_rows;
    if (n == 0) return {{}, 0};

    
    auto aos_data = convert_SoA_to_AoS(data);
    std::sort(aos_data.begin(), aos_data.end());
    auto last = std::unique(aos_data.begin(), aos_data.end());
    int distinctCount = std::distance(aos_data.begin(), last);
    if (distinctCount < k) k = distinctCount;

    std::random_device rd; std::mt19937 gen(rd());
    auto medoids = acceleratedInitialization(data, k, gen);
    k = medoids.size();
    if (k==0) return {{}, 0};
    
    std::vector<int> cluster_assignment(n);
    long long previous_total_cost = std::numeric_limits<long long>::max();
    long long final_total_cost = 0; 
    
    for (int iter = 0; iter < max_iter; ++iter) {
        long long total_cost_this_round = 0;
        for (int i = 0; i < n; ++i) {
            long long min_cost = std::numeric_limits<long long>::max();
            int assigned_medoid_idx = -1;
            auto current_point = get_point_from_soa(data, i);
            for (int m = 0; m < k; ++m) {
                long long cost = calculateTotalLogCost(current_point, medoids[m], dim);
                if (cost < min_cost) {
                    min_cost = cost;
                    assigned_medoid_idx = m;
                } else if (cost == min_cost) {
                    if (medoids[m] == current_point) {
                        assigned_medoid_idx = m;
                        break;
                    }
                }
            }
            cluster_assignment[i] = assigned_medoid_idx;
            total_cost_this_round += min_cost;
        }

        final_total_cost = total_cost_this_round; 

        if (iter > 0 && std::abs((double)previous_total_cost - total_cost_this_round) < tol) break;
        previous_total_cost = total_cost_this_round;
        
        auto new_medoids = updateMedoids(data, cluster_assignment, medoids, dim);
        if (medoids == new_medoids) break;
        medoids = new_medoids;
    }
    
    std::vector<long long> cluster_sizes(k, 0);
    for (int assignment : cluster_assignment) { if (assignment != -1 && (size_t)assignment < cluster_sizes.size()) cluster_sizes[assignment]++; }
    
    return {sortResults(medoids, cluster_assignment, cluster_sizes), final_total_cost}; 
}

AClusterResult adaptiveGreedyBasePointSelection_with_diag(const SoAData<long long>& page_data, int dim) {
    auto data_aos = convert_SoA_to_AoS(page_data);
    size_t n = data_aos.size();
    if (n == 0) return {{},0};
    std::vector<std::vector<long long>> medoids_list;
    std::unordered_set<std::vector<long long>, VectorHasher> existing_medoids_set;
    std::vector<std::unordered_set<int>> points_in_cluster_list;
    std::vector<int> point_to_leader_map(n, -1);
    medoids_list.push_back(data_aos[0]);
    existing_medoids_set.insert(data_aos[0]);
    points_in_cluster_list.emplace_back();
    points_in_cluster_list[0].insert(0);
    point_to_leader_map[0] = 0;
    for (size_t i = 1; i < n; ++i) {
        const auto& current_point = data_aos[i];
        if (existing_medoids_set.count(current_point)) {
            int existing_leader_index = -1;
            for (size_t j = 0; j < medoids_list.size(); ++j) { if (medoids_list[j] == current_point) { existing_leader_index = j; break; } }
            if (existing_leader_index != -1) { points_in_cluster_list[existing_leader_index].insert(i); point_to_leader_map[i] = existing_leader_index; }
            continue;
        }
        int best_leader_index = -1;
        long long min_cost_to_existing_leader = std::numeric_limits<long long>::max();
        for (size_t j = 0; j < medoids_list.size(); ++j) {
            long long cost = calculateTotalLogCost(current_point, medoids_list[j], dim);
            if (cost < min_cost_to_existing_leader) { min_cost_to_existing_leader = cost; best_leader_index = j; }
        }
        long long savings_from_current_point = min_cost_to_existing_leader;
        long long savings_from_reassignment = 0;
        const auto& points_in_best_leader_cluster = points_in_cluster_list[best_leader_index];
        for (int point_index_in_cluster : points_in_best_leader_cluster) {
            const auto& p = data_aos[point_index_in_cluster];
            long long cost_to_old_leader = calculateTotalLogCost(p, medoids_list[best_leader_index], dim);
            long long cost_to_new_potential_leader = calculateTotalLogCost(p, current_point, dim);
            if (cost_to_new_potential_leader < cost_to_old_leader) { savings_from_reassignment += (cost_to_old_leader - cost_to_new_potential_leader); }
        }
        long long total_savings = savings_from_current_point + savings_from_reassignment;
        long long storage_cost_for_new_point = calculateBasePointStorageCost(current_point, dim);
        if (total_savings > storage_cost_for_new_point) {
            int new_leader_id = medoids_list.size();
            medoids_list.push_back(current_point);
            existing_medoids_set.insert(current_point);
            points_in_cluster_list.emplace_back();
            points_in_cluster_list.back().insert(i);
            point_to_leader_map[i] = new_leader_id;
            std::vector<int> points_to_reassign;
            for (int point_idx : points_in_cluster_list[best_leader_index]) {
                const auto& p = data_aos[point_idx];
                long long cost_to_old = calculateTotalLogCost(p, medoids_list[best_leader_index], dim);
                long long cost_to_new = calculateTotalLogCost(p, current_point, dim);
                if (cost_to_new < cost_to_old) { points_to_reassign.push_back(point_idx); }
            }
            for (int point_idx : points_to_reassign) {
                points_in_cluster_list[best_leader_index].erase(point_idx);
                points_in_cluster_list[new_leader_id].insert(point_idx);
                point_to_leader_map[point_idx] = new_leader_id;
            }
        } else {
            points_in_cluster_list[best_leader_index].insert(i);
            point_to_leader_map[i] = best_leader_index;
        }
    }
    int k = medoids_list.size();
    auto discovered_medoids = medoids_list;
    std::vector<long long> raw_cluster_sizes(k);
    for (int i = 0; i < k; ++i) { raw_cluster_sizes[i] = points_in_cluster_list[i].size(); }
    
    long long final_total_cost = 0;
    for(size_t i = 0; i < point_to_leader_map.size(); ++i) {
        if(point_to_leader_map[i] != -1) {
            final_total_cost += calculateTotalLogCost(data_aos[i], medoids_list[point_to_leader_map[i]], dim);
        }
    }
    for(const auto& medoid : medoids_list) {
        final_total_cost += calculateBasePointStorageCost(medoid, dim);
    }
    return {sortResults(discovered_medoids, point_to_leader_map, raw_cluster_sizes), final_total_cost};
}


// original version
ClusteringResult kMedoidLogCost(const SoAData<long long>& data, int k, int max_iter, double tol, int dim) {
    int n = data.num_rows;
    if (n == 0) return {{}, {}, {}};

    auto aos_data = convert_SoA_to_AoS(data);
    std::sort(aos_data.begin(), aos_data.end());
    auto last = std::unique(aos_data.begin(), aos_data.end());
    int distinctCount = std::distance(aos_data.begin(), last);
    int original_k = k;
    if (distinctCount < k) {
        k = distinctCount;
    }
    if (GlobalDebug::DEBUG_PRINT_ENABLED) {
        std::cout << "[DEBUG|CLUSTER] Distinct data points: " << distinctCount 
                  << ". Initial k: " << original_k << ". Adjusted k: " << k << "\n";
    }

    std::random_device rd; std::mt19937 gen(rd());
    auto medoids = acceleratedInitialization(data, k, gen);
    k = medoids.size();
    if (k==0) return {};
    
    std::vector<int> cluster_assignment(n);
    long long previous_total_cost = std::numeric_limits<long long>::max();
    
    for (int iter = 0; iter < max_iter; ++iter) {
        long long total_cost_this_round = 0;
        for (int i = 0; i < n; ++i) {
            long long min_cost = std::numeric_limits<long long>::max();
            int assigned_medoid_idx = -1;
            auto current_point = get_point_from_soa(data, i);
            for (int m = 0; m < k; ++m) {
                long long cost = calculateTotalLogCost(current_point, medoids[m], dim);
                
                if (cost < min_cost) {
                    min_cost = cost;
                    assigned_medoid_idx = m;
                } else if (cost == min_cost) {
                    if (medoids[m] == current_point) {
                        assigned_medoid_idx = m;
                        break;
                    }
                }
            }
            cluster_assignment[i] = assigned_medoid_idx;
            total_cost_this_round += min_cost;
        }
        if (iter > 0 && std::abs((double)previous_total_cost - total_cost_this_round) < tol) break;
        previous_total_cost = total_cost_this_round;
        
        auto new_medoids = updateMedoids(data, cluster_assignment, medoids, dim);
        if (medoids == new_medoids) break;
        medoids = new_medoids;
    }
    
    std::vector<long long> cluster_sizes(k, 0);
    for (int assignment : cluster_assignment) { if (assignment != -1 && (size_t)assignment < cluster_sizes.size()) cluster_sizes[assignment]++; }
    return sortResults(medoids, cluster_assignment, cluster_sizes);
}

// original version
ClusteringResult adaptiveGreedyBasePointSelection(const SoAData<long long>& page_data, int dim) {
    auto data_aos = convert_SoA_to_AoS(page_data);
    size_t n = data_aos.size();
    if (n == 0) return {{}, {}, {}};
    std::vector<std::vector<long long>> medoids_list;
    std::unordered_set<std::vector<long long>, VectorHasher> existing_medoids_set;
    std::vector<std::unordered_set<int>> points_in_cluster_list;
    std::vector<int> point_to_leader_map(n, -1);
    medoids_list.push_back(data_aos[0]);
    existing_medoids_set.insert(data_aos[0]);
    points_in_cluster_list.emplace_back();
    points_in_cluster_list[0].insert(0);
    point_to_leader_map[0] = 0;
    for (size_t i = 1; i < n; ++i) {
        const auto& current_point = data_aos[i];
        if (existing_medoids_set.count(current_point)) {
            int existing_leader_index = -1;
            for (size_t j = 0; j < medoids_list.size(); ++j) { if (medoids_list[j] == current_point) { existing_leader_index = j; break; } }
            if (existing_leader_index != -1) { points_in_cluster_list[existing_leader_index].insert(i); point_to_leader_map[i] = existing_leader_index; }
            continue;
        }
        int best_leader_index = -1;
        long long min_cost_to_existing_leader = std::numeric_limits<long long>::max();
        for (size_t j = 0; j < medoids_list.size(); ++j) {
            long long cost = calculateTotalLogCost(current_point, medoids_list[j], dim);
            if (cost < min_cost_to_existing_leader) { min_cost_to_existing_leader = cost; best_leader_index = j; }
        }
        long long savings_from_current_point = min_cost_to_existing_leader;
        long long savings_from_reassignment = 0;
        const auto& points_in_best_leader_cluster = points_in_cluster_list[best_leader_index];
        for (int point_index_in_cluster : points_in_best_leader_cluster) {
            const auto& p = data_aos[point_index_in_cluster];
            long long cost_to_old_leader = calculateTotalLogCost(p, medoids_list[best_leader_index], dim);
            long long cost_to_new_potential_leader = calculateTotalLogCost(p, current_point, dim);
            if (cost_to_new_potential_leader < cost_to_old_leader) { savings_from_reassignment += (cost_to_old_leader - cost_to_new_potential_leader); }
        }
        long long total_savings = savings_from_current_point + savings_from_reassignment;
        long long storage_cost_for_new_point = calculateBasePointStorageCost(current_point, dim);
        if (total_savings > storage_cost_for_new_point) {
            int new_leader_id = medoids_list.size();
            medoids_list.push_back(current_point);
            existing_medoids_set.insert(current_point);
            points_in_cluster_list.emplace_back();
            points_in_cluster_list.back().insert(i);
            point_to_leader_map[i] = new_leader_id;
            std::vector<int> points_to_reassign;
            for (int point_idx : points_in_cluster_list[best_leader_index]) {
                const auto& p = data_aos[point_idx];
                long long cost_to_old = calculateTotalLogCost(p, medoids_list[best_leader_index], dim);
                long long cost_to_new = calculateTotalLogCost(p, current_point, dim);
                if (cost_to_new < cost_to_old) { points_to_reassign.push_back(point_idx); }
            }
            for (int point_idx : points_to_reassign) {
                points_in_cluster_list[best_leader_index].erase(point_idx);
                points_in_cluster_list[new_leader_id].insert(point_idx);
                point_to_leader_map[point_idx] = new_leader_id;
            }
        } else {
            points_in_cluster_list[best_leader_index].insert(i);
            point_to_leader_map[i] = best_leader_index;
        }
    }
    int k = medoids_list.size();
    auto discovered_medoids = medoids_list;
    std::vector<long long> raw_cluster_sizes(k);
    for (int i = 0; i < k; ++i) { raw_cluster_sizes[i] = points_in_cluster_list[i].size(); }
    return sortResults(discovered_medoids, point_to_leader_map, raw_cluster_sizes);
}


namespace {
    inline long long zigzagEncode_scalar(long long n) {
        return (n << 1) ^ (n >> 63);
    }
}

std::vector<std::vector<long long>> residualCalculationZigzag_sorted(
    const SoAData<long long>& data,
    const std::vector<std::vector<long long>>& medoids,
    const std::vector<int>& assignments,
    const std::vector<long long>& sorted_cluster_sizes,
    int dim
) {
    int n = data.num_rows;
    size_t k = medoids.size();
    if (n == 0) return {};

    std::vector<int> writePointers(k);
    int cumulativeCount = 0;
    for (size_t i = 0; i < k; ++i) {
        writePointers[i] = cumulativeCount;
        cumulativeCount += static_cast<int>(sorted_cluster_sizes[i]);
    }
    
    std::vector<std::vector<long long>> sortedResidualSeries(n, std::vector<long long>(dim));
    
    for (int i = 0; i < n; i++) {
        int clusterId = assignments[i];
        if (clusterId < 0 || static_cast<size_t>(clusterId) >= k) continue;

        const auto& base = medoids[clusterId];
        int targetIndex = writePointers[clusterId];

        for (int j = 0; j < dim; j++) {
            long long residual = data.columns[j][i] - base[j];
            sortedResidualSeries[targetIndex][j] = zigzagEncode_scalar(residual);
        }
        writePointers[clusterId]++;
    }
    return sortedResidualSeries;
}


struct BasePointsResult {
    std::vector<long long> min_base;
    std::vector<int> max_base_bit_len;
    BitBuffer base_bitstream;
};
struct EncodedDataResult {
    BitBuffer residual_bitstream;
    std::vector<std::vector<int>> pack_res_metadata;
    int pack_num;
    int total_data_points;
};
struct ClusterEncoderPageResult {
    std::vector<int> max_decimal_placements;
    std::vector<long long> min_values;
    BitBuffer fre_bitstream;
    BitBuffer residual_bitstream;
    std::vector<std::vector<int>> pack_res_metadata;
    int pack_num;
    int total_data_points;
    std::vector<long long> min_base;
    std::vector<int> max_base_bit_len;
    BitBuffer base_bitstream;
    int page_k;
};






namespace { // Anonymous helpers for encoding
    BasePointsResult generateBitstreamEncodedBasePointsOptimized(
        const std::vector<std::vector<long long>>& medoids, int dim) 
    {
        if (medoids.empty()) return {std::vector<long long>(dim, 0), std::vector<int>(dim, 0), BitBuffer()};
        size_t k = medoids.size();
        std::vector<long long> min_base(dim, std::numeric_limits<long long>::max());
        for (const auto& point : medoids) for (int i = 0; i < dim; ++i) min_base[i] = std::min(min_base[i], point[i]);
        std::vector<std::vector<long long>> offset_medoids(k, std::vector<long long>(dim));
        std::vector<int> max_base_bit_len(dim, 0);
        for (size_t i = 0; i < k; ++i) {
            for (int j = 0; j < dim; ++j) {
                long long offset_value = medoids[i][j] - min_base[j];
                offset_medoids[i][j] = offset_value;
                max_base_bit_len[j] = std::max(max_base_bit_len[j], BitBufferUtils::bitsRequiredNoSign(offset_value));
            }
        }
        BitBuffer base_bitstream;
        for (const auto& offset_point : offset_medoids) for (int i = 0; i < dim; ++i) BitBufferUtils::appendToBitstream(base_bitstream, offset_point[i], max_base_bit_len[i]);
        return {min_base, max_base_bit_len, base_bitstream};
    }
    BitBuffer generateBitstreamFrequency(
        const std::vector<long long>& cluster_size, int block_size) 
    {
        if (cluster_size.empty()) return BitBuffer();
        BitBuffer freq_bit, freq_meta;
        std::vector<long long> cluster_delta(cluster_size.size());
        if (!cluster_size.empty()) {
            cluster_delta[0] = cluster_size[0];
            for (size_t i = 1; i < cluster_size.size(); ++i) cluster_delta[i] = cluster_size[i] - cluster_size[i - 1];
        }
        int total_length = cluster_size.size();
        int num_blocks = (total_length + block_size - 1) / block_size;
        for (int block_idx = 0; block_idx < num_blocks; ++block_idx) {
            int start = block_idx * block_size;
            int end = std::min(start + block_size, total_length);
            long long max_frequency = 0;
            for (int i = start; i < end; ++i) max_frequency = std::max(max_frequency, cluster_delta[i]);
            int max_freq_bit = BitBufferUtils::bitsRequiredNoSign(max_frequency);
            BitBufferUtils::appendToBitstream(freq_meta, max_freq_bit, 8);
            for (int i = start; i < end; ++i) BitBufferUtils::appendToBitstream(freq_bit, cluster_delta[i], max_freq_bit);
        }
        freq_meta.merge(freq_bit);
        return freq_meta;
    }
    EncodedDataResult generateBitstreamEncodedData(
        const std::vector<std::vector<long long>>& residual_series, int pack_size, int dim)
    {
        int total_data_points = residual_series.size();
        if (total_data_points == 0) return {BitBuffer(), {}, 0, 0};
        int pack_num = (total_data_points + pack_size - 1) / pack_size;
        BitBuffer residual_bit;
        std::vector<std::vector<int>> pack_res_metadata(pack_num, std::vector<int>(dim));
        int pack_index = 0;
        for (int pack_start = 0; pack_start < total_data_points; pack_start += pack_size) {
            int pack_end = std::min(pack_start + pack_size, total_data_points);
            std::vector<int> max_residual_bits(dim, 0);
            for (int i = pack_start; i < pack_end; ++i) for (int j = 0; j < dim; ++j) max_residual_bits[j] = std::max(max_residual_bits[j], BitBufferUtils::bitsRequiredNoSign(residual_series[i][j]));
            pack_res_metadata[pack_index] = max_residual_bits;
            for (int i = pack_start; i < pack_end; ++i) for (int j = 0; j < dim; ++j) BitBufferUtils::appendToBitstream(residual_bit, residual_series[i][j], max_residual_bits[j]);
            pack_index++;
        }
        return {residual_bit, pack_res_metadata, pack_num, total_data_points};
    }
    struct SaveResult {
        long long total_bytes;
        long long total_header_size;
        long long total_base_size;
        long long total_fre_size;
        long long total_residual_size;
    };
    std::vector<uint8_t> assemble_final_bytestream(
        const std::vector<ClusterEncoderPageResult>& page_results,
        int dim, int pack_size, int block_size, int page_size) 
    {
        if (page_results.empty()) return {};
        BitBuffer final_output_bitstream;
        
        if (GlobalDebug::DEBUG_PRINT_ENABLED) std::cout << "\n--- [DEBUG|WRITE] Writing Global Header ---\n";
        BitBufferUtils::appendToBitstream(final_output_bitstream, dim, 8);
        BitBufferUtils::appendToBitstream(final_output_bitstream, pack_size, 16);
        BitBufferUtils::appendToBitstream(final_output_bitstream, block_size, 16);
        BitBufferUtils::appendToBitstream(final_output_bitstream, (long long)page_results.size(), 32);
        BitBufferUtils::appendToBitstream(final_output_bitstream, page_size, 16);
        if (GlobalDebug::DEBUG_PRINT_ENABLED) std::cout << "  - Global Header: dim=" << dim << ", pack=" << pack_size << ", block=" << block_size << ", page_count=" << page_results.size() << ", page_size=" << page_size << "\n";
        
        bool debug_page = false;
        for (size_t i = 0; i < page_results.size(); ++i) {
            if(i!=2){
                debug_page = false;
            }else{
                debug_page = true;
            }
            const auto& page = page_results[i];
            BitBuffer page_metadata_bitstream;
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "\n--- [DEBUG|WRITE] Writing Page " << i + 1 << " Metadata ---\n";
            
            BitBufferUtils::appendToBitstream(page_metadata_bitstream, (long long)page.page_k, 16);
            BitBufferUtils::appendToBitstream(page_metadata_bitstream, page.total_data_points, 16);
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  - Page Header: k=" << page.page_k << ", points=" << page.total_data_points << "\n";

            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  --- Writing Min Values & Decimals ---\n";
            for (int j = 0; j < dim; ++j) {
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, page.max_decimal_placements[j], 8);
                long long val = page.min_values[j];
                long long abs_val = std::abs(val);
                int bit_len = BitBufferUtils::bitsRequiredNoSign(abs_val);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, bit_len, 8);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, (val >= 0 ? 0 : 1), 1);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, abs_val, bit_len);
                if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": Decimals=" << page.max_decimal_placements[j] << ", MinValLen=" << bit_len << ", Sign=" << (val >= 0 ? 0 : 1) << ", AbsMinVal=" << abs_val << "\n";
            }
            
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  --- Writing Min Bases & Base Bit Lengths ---\n";
            for (int j = 0; j < dim; ++j) {
                long long val = page.min_base[j];
                long long abs_val = std::abs(val);
                int bit_len = BitBufferUtils::bitsRequiredNoSign(abs_val);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, bit_len, 8);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, (val >= 0 ? 0 : 1), 1);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, abs_val, bit_len);
                BitBufferUtils::appendToBitstream(page_metadata_bitstream, page.max_base_bit_len[j], 8);
                if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": MinBaseLen=" << bit_len << ", Sign=" << (val >= 0 ? 0 : 1) << ", AbsMinBase=" << abs_val << ", MaxBaseBitLen=" << page.max_base_bit_len[j] << "\n";
            }

            BitBufferUtils::appendToBitstream(page_metadata_bitstream, (long long)page.pack_res_metadata.size(), 16);
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  - Wrote Pack Metadata Count: " << page.pack_res_metadata.size() << "\n";
            for (const auto& pack : page.pack_res_metadata) for (int val : pack) BitBufferUtils::appendToBitstream(page_metadata_bitstream, val, 8);
            
            final_output_bitstream.merge(page_metadata_bitstream);
            final_output_bitstream.merge(page.base_bitstream);
            final_output_bitstream.merge(page.fre_bitstream);
            final_output_bitstream.merge(page.residual_bitstream);
        }
        return final_output_bitstream.toByteArray();
    }
}

ClusterEncoderPageResult clusterEncoderPage(
    std::unique_ptr<PreprocessorResult> preproc_res_ptr,
    const std::function<ClusteringResult(const SoAData<long long>&, int, int)>& cluster_strategy,
    int k, int pack_size, int dim, int block_size, int page_idx)
{
    if (!preproc_res_ptr) throw std::runtime_error("Preprocessor result is null.");
    auto& preproc_res = *preproc_res_ptr;

    ClusteringResult clustering_res = cluster_strategy(*preproc_res_ptr->data, k, dim);

    // FOR DEBUG
    if (GlobalDebug::DEBUG_PRINT_ENABLED && page_idx == 2) {
        std::cout << "\n--- [DEBUG|ENCODE] Page " << page_idx + 1 << " Clustering Results ---\n";

        std::cout << "--- [DEBUG|ENCODE] Medoids (Top 10) ---\n";
        int medoids_to_print = std::min((int)clustering_res.medoids.size(), 10);
        for (int m = 0; m < medoids_to_print; ++m) {
            std::cout << "  - Medoid " << std::setw(3) << m << ": { ";
            for(long long val : clustering_res.medoids[m]) {
                std::cout << val << " ";
            }
            std::cout << "}\n";
        }

        std::cout << "--- [DEBUG|ENCODE] Frequencies (Top 10) ---\n";
        int freqs_to_print = std::min((int)clustering_res.cluster_sizes.size(), 10);
        for (int m = 0; m < freqs_to_print; ++m) {
            std::cout << "  - Freq for Medoid " << std::setw(3) << m << ": " << clustering_res.cluster_sizes[m] << "\n";
        }
    }
    // DEBUG END

    int page_k = clustering_res.medoids.size();
    BitBuffer fre_bitstream = generateBitstreamFrequency(clustering_res.cluster_sizes, block_size);
    
    std::vector<std::vector<long long>> residuals = residualCalculationZigzag_sorted(
        *preproc_res_ptr->data, clustering_res.medoids, clustering_res.cluster_assignment, 
        clustering_res.cluster_sizes, dim);

    EncodedDataResult bitstream_res = generateBitstreamEncodedData(residuals, pack_size, dim);
    BasePointsResult basepoint_res = generateBitstreamEncodedBasePointsOptimized(clustering_res.medoids, dim);
    
    return {
        preproc_res.max_decimal_places, preproc_res.min_values,
        fre_bitstream, bitstream_res.residual_bitstream,
        bitstream_res.pack_res_metadata, bitstream_res.pack_num,
        bitstream_res.total_data_points, basepoint_res.min_base,
        basepoint_res.max_base_bit_len, basepoint_res.base_bitstream,
        page_k
    };
}

std::pair<ClusterEncoderPageResult, PageDiagnostics> clusterEncoderPage_with_diag(
    std::unique_ptr<PreprocessorResult> preproc_res_ptr,
    const std::string& algorithm_name,
    int k, int pack_size, int dim, int block_size, int page_idx)
{
    if (!preproc_res_ptr) throw std::runtime_error("Preprocessor result is null.");
    auto& preproc_res = *preproc_res_ptr;

    PageDiagnostics diag;
    diag.page_index = page_idx;
    diag.algorithm_name = algorithm_name;
    auto& timings = diag.timings;

    // --- 1. Clustering with timing and cost retrieval ---
    auto start_cluster = std::chrono::high_resolution_clock::now();
    ClusteringResult clustering_res;
    if (algorithm_name == "KCluster") {
        auto k_res = kMedoidLogCost_with_diag(*preproc_res.data, k, 10, 1.0e-3, dim);
        clustering_res = k_res.result;
        diag.kcluster_costs.final_residual_cost_bits = k_res.final_cost;
    } else if (algorithm_name == "ACluster") {
        auto a_res = adaptiveGreedyBasePointSelection_with_diag(*preproc_res.data, dim);
        clustering_res = a_res.result;
        diag.acluster_costs.final_total_cost_bits = a_res.final_cost;
    } else {
        // Fallback for any other algorithm if needed
        clustering_res = kMedoidLogCost(*preproc_res.data, k, 10, 1.0e-3, dim); // Example
    }
    auto end_cluster = std::chrono::high_resolution_clock::now();
    timings.clustering_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_cluster - start_cluster).count();
    
    // --- 2. Residual calculation with timing ---
    auto start_residual = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<long long>> residuals = residualCalculationZigzag_sorted(
        *preproc_res.data, clustering_res.medoids, clustering_res.cluster_assignment, 
        clustering_res.cluster_sizes, dim);
    auto end_residual = std::chrono::high_resolution_clock::now();
    timings.residual_calc_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_residual - start_residual).count();
    
    // --- 3. Bitstream generation with timing and size calculation ---
    auto start_bitstream = std::chrono::high_resolution_clock::now();
    
    int page_k = clustering_res.medoids.size();
    BitBuffer fre_bitstream = generateBitstreamFrequency(clustering_res.cluster_sizes, block_size);
    EncodedDataResult bitstream_res = generateBitstreamEncodedData(residuals, pack_size, dim);
    BasePointsResult basepoint_res = generateBitstreamEncodedBasePointsOptimized(clustering_res.medoids, dim);
    
    auto end_bitstream = std::chrono::high_resolution_clock::now();
    timings.bitstream_gen_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_bitstream - start_bitstream).count();

    // --- 4. Populate component sizes in diagnostics ---
    diag.sizes.frequencies_bits = fre_bitstream.size();
    diag.sizes.residuals_bits = bitstream_res.residual_bitstream.size();
    diag.sizes.medoids_bits = basepoint_res.base_bitstream.size();

    // Construct the page result for assembly
    ClusterEncoderPageResult page_result = {
        preproc_res.max_decimal_places, preproc_res.min_values,
        std::move(fre_bitstream), std::move(bitstream_res.residual_bitstream),
        std::move(bitstream_res.pack_res_metadata), bitstream_res.pack_num,
        bitstream_res.total_data_points, std::move(basepoint_res.min_base),
        std::move(basepoint_res.max_base_bit_len), std::move(basepoint_res.base_bitstream),
        page_k
    };
    
    return {std::move(page_result), std::move(diag)};
}



std::pair<std::vector<uint8_t>, GlobalDiagnostics> encode_multidim_impl_with_diag(
    const std::string& csv_file_path, int& dim,
    const std::string& algorithm_name,
    int k, int pack_size, int block_size, int page_size) 
{
    std::vector<ClusterEncoderPageResult> all_page_results;
    GlobalDiagnostics global_diag;

    std::ifstream file_stream(csv_file_path);
    if (!file_stream.is_open()) throw std::runtime_error("Cannot open data file: " + csv_file_path);
    std::vector<std::string> page_lines;
    page_lines.reserve(page_size);
    std::string line;
    bool file_is_finished = false;

    while (!file_is_finished) {
        page_lines.clear();
        while (page_lines.size() < (size_t)page_size && std::getline(file_stream, line)) if (!line.empty()) page_lines.push_back(line);
        if (file_stream.eof()) file_is_finished = true;
        if (page_lines.empty()) break;

        auto start_preproc = std::chrono::high_resolution_clock::now();
        auto preproc_res_ptr = std::make_unique<PreprocessorResult>(preprocessPageFromStrings(page_lines, dim));
        auto end_preproc = std::chrono::high_resolution_clock::now();
        
        if (!preproc_res_ptr || !preproc_res_ptr->data || preproc_res_ptr->data->num_rows == 0) continue;
        
        int current_dim = preproc_res_ptr->data->dim;
        int page_index_counter = all_page_results.size();

        auto [page_result, page_diag] = clusterEncoderPage_with_diag(
            std::move(preproc_res_ptr), algorithm_name, k, pack_size, current_dim, block_size, page_index_counter);
        
        page_diag.timings.preprocessing_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_preproc - start_preproc).count();
        
        // Aggregate diagnostics
        global_diag.total_timings.preprocessing_ns += page_diag.timings.preprocessing_ns;
        global_diag.total_timings.clustering_ns += page_diag.timings.clustering_ns;
        global_diag.total_timings.residual_calc_ns += page_diag.timings.residual_calc_ns;
        global_diag.total_timings.bitstream_gen_ns += page_diag.timings.bitstream_gen_ns;

        global_diag.total_sizes.medoids_bits += page_diag.sizes.medoids_bits;
        global_diag.total_sizes.frequencies_bits += page_diag.sizes.frequencies_bits;
        global_diag.total_sizes.residuals_bits += page_diag.sizes.residuals_bits;

        if (GlobalDebug::DEBUG_PRINT_ENABLED) {
            if (algorithm_name == "KCluster") {
                std::cout << "[DIAG] Page " << page_index_counter << " KCluster Final Cost: " << page_diag.kcluster_costs.final_residual_cost_bits << "\n";
            } else if (algorithm_name == "ACluster") {
                std::cout << "[DIAG] Page " << page_index_counter << " ACluster Final Cost: " << page_diag.acluster_costs.final_total_cost_bits << "\n";
            }
        }
        
        all_page_results.push_back(std::move(page_result));
    }

    if (all_page_results.empty()) return {{}, {}};

    // Calculate header size
    BitBuffer header_buf;
    BitBufferUtils::appendToBitstream(header_buf, dim, 8);
    BitBufferUtils::appendToBitstream(header_buf, pack_size, 16);
    BitBufferUtils::appendToBitstream(header_buf, block_size, 16);
    BitBufferUtils::appendToBitstream(header_buf, (long long)all_page_results.size(), 32);
    BitBufferUtils::appendToBitstream(header_buf, page_size, 16);
    global_diag.total_sizes.header_bits += header_buf.size();

    for (const auto& page : all_page_results) {
        BitBuffer page_header_buf;
        BitBufferUtils::appendToBitstream(page_header_buf, (long long)page.page_k, 16);
        BitBufferUtils::appendToBitstream(page_header_buf, page.total_data_points, 16);
        for (int j = 0; j < dim; ++j) {
            page_header_buf.currentBitPos += 8; // max_decimal_places
            page_header_buf.currentBitPos += 8 + 1 + 32; // min_values (estimate)
            page_header_buf.currentBitPos += 8 + 1 + 32; // min_base (estimate)
            page_header_buf.currentBitPos += 8; // max_base_bit_len
        }
        page_header_buf.currentBitPos += 16; // pack_res_metadata.size()
        for(const auto& pack : page.pack_res_metadata) page_header_buf.currentBitPos += pack.size() * 8;
        global_diag.total_sizes.header_bits += page_header_buf.size();
    }
    
    return {assemble_final_bytestream(all_page_results, dim, pack_size, block_size, page_size), global_diag};
}

std::vector<uint8_t> encode_multidim_impl(
    const std::string& csv_file_path, int& dim,
    const std::function<ClusteringResult(const SoAData<long long>&, int, int)>& cluster_strategy,
    int k, int pack_size, int block_size, int page_size) 
{
    std::vector<ClusterEncoderPageResult> all_page_results;
    std::ifstream file_stream(csv_file_path);
    if (!file_stream.is_open()) throw std::runtime_error("Cannot open data file: " + csv_file_path);
    std::vector<std::string> page_lines;
    page_lines.reserve(page_size);
    std::string line;
    bool file_is_finished = false;
    while (!file_is_finished) {
        page_lines.clear();
        while (page_lines.size() < (size_t)page_size && std::getline(file_stream, line)) if (!line.empty()) page_lines.push_back(line);
        if (file_stream.eof()) file_is_finished = true;
        if (page_lines.empty()) break;
        auto preproc_res_ptr = std::make_unique<PreprocessorResult>(preprocessPageFromStrings(page_lines, dim));
        if (!preproc_res_ptr || !preproc_res_ptr->data || preproc_res_ptr->data->num_rows == 0) continue;
        int current_dim = preproc_res_ptr->data->dim;
        int page_index_counter = all_page_results.size();
        ClusterEncoderPageResult page_result = clusterEncoderPage(std::move(preproc_res_ptr), cluster_strategy, k, pack_size, current_dim, block_size, page_index_counter);
        all_page_results.push_back(std::move(page_result));
    }
    if (all_page_results.empty()) return {};
    return assemble_final_bytestream(all_page_results, dim, pack_size, block_size, page_size);
}

std::vector<std::vector<double>> decode_multidim_impl(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.empty()) return {};
    BitStreamReader reader(compressed_data);

    if (GlobalDebug::DEBUG_PRINT_ENABLED) std::cout << "\n--- [DEBUG|READ] Starting Decoding Process ---\n";
    int dim = reader.readBits(8);
    int pack_size = reader.readBits(16);
    int block_size = reader.readBits(16);
    int page_count = reader.readBits(32);
    int page_size = reader.readBits(16);
    (void)page_size;
    if (GlobalDebug::DEBUG_PRINT_ENABLED) std::cout << "  - Read Global Header: dim=" << dim << ", pack=" << pack_size << ", block=" << block_size << ", page_count=" << page_count << ", page_size=" << page_size << "\n";

    std::vector<double> pow10_lookup; for (int i = 0; i < 20; ++i) pow10_lookup.push_back(std::pow(10, i));
    std::vector<std::vector<double>> all_data_rows;
    
    bool debug_page = false;
    for (int i = 0; i < page_count; ++i) {
        if(i!=2){
            debug_page = false;
        }else{
            debug_page = true;
        }
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "\n--- [DEBUG|READ] Decoding Page " << i + 1 << " ---\n";
        int k = reader.readBits(16);
        int page_data_points = reader.readBits(16);
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  - Read Page Header: k=" << k << ", points=" << page_data_points << "\n";

        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  --- Reading Min Values & Decimals ---\n";
        std::vector<int> max_decimals(dim);
        std::vector<long long> min_values(dim);
        for (int j = 0; j < dim; ++j) {
            max_decimals[j] = reader.readBits(8);
            int bit_len = reader.readBits(8);
            long long sign = reader.readBits(1);
            long long abs_val = reader.readBits(bit_len);
            min_values[j] = (sign == 0) ? abs_val : -abs_val;
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": Decimals=" << max_decimals[j] << ", MinValLen=" << bit_len << ", Sign=" << sign << ", AbsMinVal=" << abs_val << "\n";
        }
        
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  --- Reading Min Bases & Base Bit Lengths ---\n";
        std::vector<long long> min_bases(dim);
        std::vector<int> max_base_bit_len(dim);
        for (int j = 0; j < dim; ++j) {
            int bit_len = reader.readBits(8);
            long long sign = reader.readBits(1);
            long long abs_val = reader.readBits(bit_len);
            min_bases[j] = (sign == 0) ? abs_val : -abs_val;
            max_base_bit_len[j] = reader.readBits(8);
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": MinBaseLen=" << bit_len << ", Sign=" << sign << ", AbsMinBase=" << abs_val << ", MaxBaseBitLen=" << max_base_bit_len[j] << "\n";
        }
        
        int pack_num = reader.readBits(16);
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  - Read Pack Metadata Count: " << pack_num << "\n";
        std::vector<std::vector<int>> pack_metadata(pack_num, std::vector<int>(dim));
        for (int p = 0; p < pack_num; ++p) for (int j = 0; j < dim; ++j) pack_metadata[p][j] = reader.readBits(8);
        
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "pack_metadata" << "\n";
        std::vector<std::vector<long long>> medoids_long(k, std::vector<long long>(dim));
        for (int m = 0; m < k; ++m) for (int j = 0; j < dim; ++j) medoids_long[m][j] = reader.readBits(max_base_bit_len[j]) + min_bases[j];
        
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
            std::cout << "\n--- [DEBUG|DECODE] Decoded Medoids ---\n";
            for (int m = 0; m < 10; ++m) {
                std::cout << "  - Medoid " << std::setw(3) << m << ": { ";
                for(long long val : medoids_long[m]) {
                    std::cout << val << " ";
                }
                std::cout << "}\n";
            }
        }

        std::vector<long long> cluster_sizes(k, 0);
        if (k > 0) {
            int num_freq_blocks = (k + block_size - 1) / block_size;
            
            // --- Step 1: Read ALL metadata (max_bits) first ---
            std::vector<int> max_bits_per_block(num_freq_blocks);
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
                std::cout << "\n--- [DEBUG|FREQ_READ] Reading all " << num_freq_blocks << " max_bit metadata blocks ---\n";
            }
            for (int b = 0; b < num_freq_blocks; ++b) {
                max_bits_per_block[b] = reader.readBits(8);
                if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
                    std::cout << "  - Block " << b << ": Will read with max_bit=" << max_bits_per_block[b] << "\n";
                }
            }

            // --- Step 2: Now, read ALL delta data using the stored metadata ---
            std::vector<long long> deltas(k);
            int freq_idx = 0;
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
                std::cout << "--- [DEBUG|FREQ_READ] Reading all delta data blocks ---\n";
            }
            for (int b = 0; b < num_freq_blocks && freq_idx < k; ++b) {
                int max_bit_for_this_block = max_bits_per_block[b];
                int end = std::min((b + 1) * block_size, k);
                for (int j = freq_idx; j < end; ++j) {
                    deltas[j] = reader.readBits(max_bit_for_this_block);
                }
                freq_idx = end;
            }

            // --- Step 3: Reconstruct the cumulative sum (This part is correct) ---
            cluster_sizes[0] = deltas[0];
            for (size_t j = 1; j < deltas.size(); ++j) {
                cluster_sizes[j] = cluster_sizes[j - 1] + deltas[j];
            }
        }

        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
            std::cout << "\n--- [DEBUG|DECODE] Decoded Frequencies (Cluster Sizes) ---\n";
            for (int m = 0; m < 10; ++m) {
                std::cout << "  - Freq for Medoid " << std::setw(3) << m << ": " << cluster_sizes[m] << "\n";
            }
        }

        std::vector<std::vector<long long>> residual_series(page_data_points, std::vector<long long>(dim));
        int data_counter = 0;
        if (page_data_points > 0) {
            for (int p = 0; p < pack_num; ++p) {
                const auto& bits = pack_metadata[p];
                int points_in_pack = (p == pack_num - 1) ? (page_data_points - (p * pack_size)) : pack_size;
                for (int pt = 0; pt < points_in_pack; ++pt) {
                    if (data_counter >= page_data_points) break;
                    for (int j = 0; j < dim; ++j) residual_series[data_counter][j] = zigzagDecode(reader.readBits(bits[j]));
                    data_counter++;
                }
            }
        }
        
        if (page_data_points > 0 && k > 0) {
            int current_point_idx = 0;
            for (int medoid_idx = 0; medoid_idx < k; ++medoid_idx) {
                // Get the number of points IN THIS SPECIFIC CLUSTER
                long long points_in_this_cluster = cluster_sizes[medoid_idx];
                const auto& base_point = medoids_long[medoid_idx];

                for (int p_count = 0; p_count < points_in_this_cluster; ++p_count) {
                    if (current_point_idx < page_data_points) {
                        std::vector<double> row(dim);
                        for (int j = 0; j < dim; ++j) {
                            long long int_val = base_point[j] + residual_series[current_point_idx][j] + min_values[j];
                            row[j] = static_cast<double>(int_val) / pow10_lookup[max_decimals[j]];
                        }
                        all_data_rows.push_back(row);
                        current_point_idx++;
                    }
                }
            }
        } else if (page_data_points > 0) {
            for (int d = 0; d < page_data_points; ++d) {
                std::vector<double> row(dim);
                for (int j = 0; j < dim; ++j) row[j] = static_cast<double>(residual_series[d][j] + min_values[j]) / pow10_lookup[max_decimals[j]];
                all_data_rows.push_back(row);
            }
        }
    }
    return all_data_rows;
}


} // namespace ClusterLogic