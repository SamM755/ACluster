#include "cluster_common.h"
#include "cluster_encoder_logic.h"
#include "optimal_1d.h"
#include "cluster_1d.h"
#include <functional>
#include <chrono>
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
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <fstream>
#include <cctype>
#include <cstdlib>

namespace {
namespace GlobalDebug {
    bool DEBUG_PRINT_ENABLED = false;

}
}

namespace ClusterLogic {

namespace {
bool stream_debug_enabled() {
    const char* value = std::getenv("CLUSTER_STREAM_DEBUG");
    return value != nullptr && std::string(value) == "1";
}
}

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

namespace {
    void parse_csv_line(const std::string& line, std::vector<std::string>& cells) {
        std::stringstream line_stream(line);
        std::string cell;
        cells.clear();
        while (std::getline(line_stream, cell, ',')) cells.push_back(cell);
    }

    struct ParsedDecimal {
        bool negative = false;
        std::string digits;
        long long fractional_digits = 0;
    };

    constexpr int kMaxStoredDecimalPlaces = 18;
    constexpr int kSafeScaledMagnitudeDigits = 18;

    ParsedDecimal parse_decimal_string(const std::string& text) {
        size_t begin = 0;
        size_t end = text.size();
        while (begin < end && std::isspace(static_cast<unsigned char>(text[begin]))) ++begin;
        while (end > begin && std::isspace(static_cast<unsigned char>(text[end - 1]))) --end;
        if (begin == end) throw std::invalid_argument("empty numeric cell");

        ParsedDecimal parsed;
        size_t pos = begin;
        if (text[pos] == '+' || text[pos] == '-') parsed.negative = text[pos++] == '-';
        bool saw_dot = false;
        long long fractional_digits = 0;
        while (pos < end && text[pos] != 'e' && text[pos] != 'E') {
            const char ch = text[pos++];
            if (ch == '.') {
                if (saw_dot) throw std::invalid_argument("multiple decimal points");
                saw_dot = true;
            } else if (std::isdigit(static_cast<unsigned char>(ch))) {
                parsed.digits.push_back(ch);
                if (saw_dot) ++fractional_digits;
            } else {
                throw std::invalid_argument("invalid numeric character");
            }
        }
        if (parsed.digits.empty()) throw std::invalid_argument("numeric cell has no digits");

        long long exponent = 0;
        if (pos < end) {
            ++pos;
            bool exponent_negative = false;
            if (pos < end && (text[pos] == '+' || text[pos] == '-')) exponent_negative = text[pos++] == '-';
            const size_t exponent_begin = pos;
            while (pos < end && std::isdigit(static_cast<unsigned char>(text[pos]))) {
                const int digit = text[pos++] - '0';
                if (exponent > 1000000) throw std::out_of_range("scientific exponent is too large");
                exponent = exponent * 10 + digit;
            }
            if (pos != end || exponent_begin == pos) throw std::invalid_argument("invalid scientific exponent");
            if (exponent_negative) exponent = -exponent;
        }
        parsed.fractional_digits = fractional_digits - exponent;
        return parsed;
    }

    int significant_digit_count(const ParsedDecimal& parsed) {
        const size_t first = parsed.digits.find_first_not_of('0');
        return first == std::string::npos ? 0 : static_cast<int>(parsed.digits.size() - first);
    }

    long long stringToLongRounded(const std::string& s, int decimal_places_to_scale) {
        const ParsedDecimal parsed = parse_decimal_string(s);
        const long long shift = static_cast<long long>(decimal_places_to_scale) - parsed.fractional_digits;
        if (shift < -1000000 || shift > 1000000) throw std::out_of_range("invalid decimal scaling shift");

        unsigned long long magnitude = 0;
        const auto max_value = static_cast<unsigned long long>(std::numeric_limits<long long>::max());
        const size_t first = parsed.digits.find_first_not_of('0');
        if (first == std::string::npos) return 0;
        const std::string_view digits(parsed.digits.data() + first, parsed.digits.size() - first);
        const auto append_digit = [&](unsigned digit) {
            if (magnitude > (max_value - digit) / 10) throw std::out_of_range("scaled value exceeds int64 range");
            magnitude = magnitude * 10 + digit;
        };

        if (shift >= 0) {
            for (char ch : digits) append_digit(static_cast<unsigned>(ch - '0'));
            for (long long i = 0; i < shift; ++i) {
                if (magnitude > max_value / 10) throw std::out_of_range("scaled value exceeds int64 range");
                magnitude *= 10;
            }
        } else {
            const long long dropped_digits = -shift;
            if (dropped_digits < static_cast<long long>(digits.size())) {
                const size_t kept_digits = digits.size() - static_cast<size_t>(dropped_digits);
                for (size_t i = 0; i < kept_digits; ++i) append_digit(static_cast<unsigned>(digits[i] - '0'));
                if (digits[kept_digits] >= '5') {
                    if (magnitude == max_value) throw std::out_of_range("scaled value exceeds int64 range");
                    ++magnitude;
                }
            } else if (dropped_digits == static_cast<long long>(digits.size()) && digits.front() >= '5') {
                magnitude = 1;
            }
        }

        const long long value = static_cast<long long>(magnitude);
        return parsed.negative ? -value : value;
    }
}
PreprocessorResult preprocessPageFromStrings(const std::vector<std::string>& csv_page, int& dim_ref) {
    if (csv_page.empty()) return {nullptr, {}, {}};
    std::vector<int> max_decimal_places;
    std::vector<int> max_safe_decimal_places;
    std::vector<std::vector<std::string>> page_as_strings;
    page_as_strings.reserve(csv_page.size());
    int num_rows = 0;
    bool is_dim_known = (dim_ref > 0);
    if (is_dim_known) {
        max_decimal_places.assign(dim_ref, 0);
        max_safe_decimal_places.assign(dim_ref, kMaxStoredDecimalPlaces);
    }
    for (const auto& line : csv_page) {
        if (line.empty()) continue;
        std::vector<std::string> cells;
        parse_csv_line(line, cells);
        if (!is_dim_known) {
            dim_ref = cells.size();
            if (dim_ref == 0) continue;
            max_decimal_places.assign(dim_ref, 0);
            max_safe_decimal_places.assign(dim_ref, kMaxStoredDecimalPlaces);
            is_dim_known = true;
        }
        if (cells.size() != static_cast<size_t>(dim_ref)) continue;
        for (int j = 0; j < dim_ref; ++j) {
            const ParsedDecimal parsed = parse_decimal_string(cells[j]);
            if (parsed.fractional_digits > std::numeric_limits<int>::max()) {
                throw std::out_of_range("effective decimal precision exceeds supported range");
            }
            max_decimal_places[j] = std::max(max_decimal_places[j], static_cast<int>(std::max(0LL, parsed.fractional_digits)));
            const int significant_digits = significant_digit_count(parsed);
            if (significant_digits > 0) {
                const long long safe_scale = parsed.fractional_digits
                    + static_cast<long long>(kSafeScaledMagnitudeDigits - significant_digits);
                if (safe_scale < 0) {
                    throw std::out_of_range("value exceeds int64 range even without decimal scaling");
                }
                max_safe_decimal_places[j] = std::min(
                    max_safe_decimal_places[j],
                    static_cast<int>(std::min<long long>(safe_scale, kMaxStoredDecimalPlaces))
                );
            }
        }
        page_as_strings.push_back(cells);
        num_rows++;
    }
    if (num_rows == 0) return {nullptr, {}, {}};
    for (int j = 0; j < dim_ref; ++j) {
        const int requested_scale = max_decimal_places[j];
        const int safe_scale = max_safe_decimal_places[j];
        if (requested_scale > safe_scale) {
            std::cerr << "Cluster preprocessing: column " << j
                      << " reduced decimal scale from " << requested_scale
                      << " to " << safe_scale
                      << " to keep scaled values within int64." << std::endl;
            max_decimal_places[j] = safe_scale;
        }
    }
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
            long long val = stringToLongRounded(page_as_strings[i][j], max_decimal_places[j]);
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
    inline long long calculateScalarResidualCost(long long x, long long b) {
        long long residual = x - b;
        return bitLength(residual >= 0 ? residual : -residual) + 1;
    }
    inline long long calculateScalarBaseCost(long long b) {
        return bitLength(b >= 0 ? b : -b) + 1;
    }
    long long manhattanDistance(const std::vector<long long>& p1, const std::vector<long long>& p2) {
        long long sum = 0;
        for (size_t i = 0; i < p1.size(); ++i) { sum += std::abs(p1[i] - p2[i]); }
        return sum;
    }

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

        std::unordered_set<std::vector<long long>, VectorHasher> selected_medoids_set;

        std::uniform_int_distribution<> distrib(0, data.num_rows - 1);

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

            while(selected_medoids_set.count(get_point_from_soa(data, chosen_idx))) {
                chosen_idx = (chosen_idx + 1) % data.num_rows;
            }

            auto new_medoid = get_point_from_soa(data, chosen_idx);
            medoids.push_back(new_medoid);
            selected_medoids_set.insert(new_medoid);

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

KMedoidResult kMedoidLogCost_with_diag(const SoAData<long long>& data, int k, int max_iter, double tol, int dim,
                                    const uint32_t* seed) {
    int n = data.num_rows;
    if (n == 0) return {{}, 0};

    auto aos_data = convert_SoA_to_AoS(data);
    std::sort(aos_data.begin(), aos_data.end());
    auto last = std::unique(aos_data.begin(), aos_data.end());
    int distinctCount = std::distance(aos_data.begin(), last);
    if (distinctCount < k) k = distinctCount;

    std::random_device rd; std::mt19937 gen(seed ? *seed : rd());
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

AClusterResult adaptiveGreedyBasePointSelection_with_diag(const SoAData<long long>& page_data, int dim, int adjacent_reference) {
    auto data_aos = convert_SoA_to_AoS(page_data);
    size_t n = data_aos.size();
    if (n == 0) return {{},0};
    adjacent_reference = std::max(1, adjacent_reference);
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
        std::vector<std::pair<long long, int>> cost_rank;
        cost_rank.reserve(medoids_list.size());
        int best_leader_index = -1;
        long long min_cost_to_existing_leader = std::numeric_limits<long long>::max();
        for (size_t j = 0; j < medoids_list.size(); ++j) {
            long long cost = calculateTotalLogCost(current_point, medoids_list[j], dim);
            cost_rank.push_back({cost, static_cast<int>(j)});
            if (cost < min_cost_to_existing_leader) { min_cost_to_existing_leader = cost; best_leader_index = static_cast<int>(j); }
        }
        std::sort(cost_rank.begin(), cost_rank.end(), [](const auto& a, const auto& b) {
            return a.first != b.first ? a.first < b.first : a.second < b.second;
        });
        int considered_clusters = std::min(static_cast<int>(cost_rank.size()), adjacent_reference);
        std::vector<int> candidate_cluster_indices;
        candidate_cluster_indices.reserve(considered_clusters);
        for (int ci = 0; ci < considered_clusters; ++ci) {
            candidate_cluster_indices.push_back(cost_rank[ci].second);
        }
        long long savings_from_current_point = min_cost_to_existing_leader;
        long long savings_from_reassignment = 0;
        for (int cluster_idx : candidate_cluster_indices) {
            const auto& points_in_cluster = points_in_cluster_list[cluster_idx];
            for (int point_index_in_cluster : points_in_cluster) {
                const auto& p = data_aos[point_index_in_cluster];
                long long cost_to_old_leader = calculateTotalLogCost(p, medoids_list[cluster_idx], dim);
                long long cost_to_new_potential_leader = calculateTotalLogCost(p, current_point, dim);
                if (cost_to_new_potential_leader < cost_to_old_leader) { savings_from_reassignment += (cost_to_old_leader - cost_to_new_potential_leader); }
            }
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
            for (int cluster_idx : candidate_cluster_indices) {
                std::vector<int> points_to_reassign;
                for (int point_idx : points_in_cluster_list[cluster_idx]) {
                    const auto& p = data_aos[point_idx];
                    long long cost_to_old = calculateTotalLogCost(p, medoids_list[cluster_idx], dim);
                    long long cost_to_new = calculateTotalLogCost(p, current_point, dim);
                    if (cost_to_new < cost_to_old) { points_to_reassign.push_back(point_idx); }
                }
                for (int point_idx : points_to_reassign) {
                    points_in_cluster_list[cluster_idx].erase(point_idx);
                    points_in_cluster_list[new_leader_id].insert(point_idx);
                    point_to_leader_map[point_idx] = new_leader_id;
                }
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

ReferenceRedundancyDiagnostics analyzeReferenceRedundancy(
    const SoAData<long long>& page_data,
    const ClusteringResult& clustering_result,
    int dim)
{
    ReferenceRedundancyDiagnostics diag;
    const auto data_aos = convert_SoA_to_AoS(page_data);
    const auto& references = clustering_result.medoids;
    const int reference_count = static_cast<int>(references.size());
    diag.reference_count = reference_count;

    if (reference_count == 0 || data_aos.empty()) {
        return diag;
    }

    std::vector<long long> reference_costs(reference_count, 0);
    std::vector<long long> replacement_penalties(reference_count, 0);
    for (int reference_index = 0; reference_index < reference_count; ++reference_index) {
        reference_costs[reference_index] = calculateBasePointStorageCost(
            references[reference_index], dim
        );
        diag.reference_cost_bits += reference_costs[reference_index];
    }

    for (const auto& point : data_aos) {
        long long best_cost = std::numeric_limits<long long>::max();
        long long second_best_cost = std::numeric_limits<long long>::max();
        int best_reference_index = -1;

        for (int reference_index = 0; reference_index < reference_count; ++reference_index) {
            const long long cost = calculateTotalLogCost(
                point, references[reference_index], dim
            );
            if (cost < best_cost) {
                second_best_cost = best_cost;
                best_cost = cost;
                best_reference_index = reference_index;
            } else if (cost < second_best_cost) {
                second_best_cost = cost;
            }
        }

        if (best_reference_index < 0) {
            throw std::runtime_error("Reference redundancy analysis found no reference for a data point.");
        }
        diag.residual_cost_bits += best_cost;
        if (reference_count > 1) {
            replacement_penalties[best_reference_index] += second_best_cost - best_cost;
        }
    }

    if (reference_count == 1) {
        return diag;
    }

    for (int reference_index = 0; reference_index < reference_count; ++reference_index) {
        const long long removal_delta = replacement_penalties[reference_index]
            - reference_costs[reference_index];
        if (removal_delta <= 0) {
            ++diag.individually_removable_reference_count;
            diag.removable_reference_cost_bits += reference_costs[reference_index];
        }
    }
    return diag;
}

AOptimalResult adaptiveOptimal1D_with_diag(const SoAData<long long>& page_data, int dim) {
    if (dim != 1) throw std::invalid_argument("AOptimal only supports 1D data.");
    auto optimal = optimal1D(page_data);
    return {std::move(optimal.result), optimal.objective};
}

ClusteringResult sortClusteringBySize(std::vector<std::vector<long long>>& medoids,
    std::vector<int>& assignment, std::vector<long long>& sizes) {
    return sortResults(medoids, assignment, sizes);
}

ClusteringResult adaptiveOptimal1D(const SoAData<long long>& page_data, int dim) {
    return adaptiveOptimal1D_with_diag(page_data, dim).result;
}

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

ClusteringResult kMeansEuclidean(const SoAData<long long>& data, int k, int max_iter, double tol, int dim) {
    int n = data.num_rows;
    if (n == 0) return {{}, {}, {}};

    auto aos_data = convert_SoA_to_AoS(data);
    std::sort(aos_data.begin(), aos_data.end());
    auto last = std::unique(aos_data.begin(), aos_data.end());
    int distinctCount = std::distance(aos_data.begin(), last);
    if (distinctCount < k) k = distinctCount;

    std::random_device rd;
    std::mt19937 gen(rd());
    auto all_points = convert_SoA_to_AoS(data);
    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), gen);
    std::vector<std::vector<long long>> centroids;
    centroids.reserve(k);
    std::unordered_set<std::vector<long long>, VectorHasher> seen_centroids;
    for (int idx : indices) {
        const auto& p = all_points[idx];
        if (seen_centroids.insert(p).second) {
            centroids.push_back(p);
            if (static_cast<int>(centroids.size()) >= k) break;
        }
    }
    k = static_cast<int>(centroids.size());
    if (k == 0) return {{}, {}, {}};

    std::vector<int> assignment(n, -1);
    auto squared_l2 = [&](const std::vector<long long>& p, const std::vector<long long>& c) {
        long double sum = 0.0L;
        for (int d = 0; d < dim; ++d) {
            long double diff = static_cast<long double>(p[d]) - static_cast<long double>(c[d]);
            sum += diff * diff;
        }
        return sum;
    };

    for (int iter = 0; iter < max_iter; ++iter) {
        for (int i = 0; i < n; ++i) {
            const auto& p = all_points[i];
            long double best_dist = std::numeric_limits<long double>::infinity();
            int best_idx = 0;
            for (int ci = 0; ci < k; ++ci) {
                long double d2 = squared_l2(p, centroids[ci]);
                if (d2 < best_dist) {
                    best_dist = d2;
                    best_idx = ci;
                } else if (d2 == best_dist && centroids[ci] == p) {
                    best_idx = ci;
                    break;
                }
            }
            assignment[i] = best_idx;
        }

        std::vector<std::vector<long double>> sums(k, std::vector<long double>(dim, 0.0L));
        std::vector<int> counts(k, 0);
        for (int i = 0; i < n; ++i) {
            int ci = assignment[i];
            if (ci < 0 || ci >= k) continue;
            counts[ci] += 1;
            const auto& p = all_points[i];
            for (int d = 0; d < dim; ++d) sums[ci][d] += static_cast<long double>(p[d]);
        }

        std::vector<std::vector<long long>> new_centroids = centroids;
        long double total_shift = 0.0L;
        for (int ci = 0; ci < k; ++ci) {
            if (counts[ci] == 0) continue;
            for (int d = 0; d < dim; ++d) {
                long double mean = sums[ci][d] / static_cast<long double>(counts[ci]);
                long long rounded = static_cast<long long>(std::llround(mean));
                total_shift += std::abs(static_cast<long double>(rounded) - static_cast<long double>(centroids[ci][d]));
                new_centroids[ci][d] = rounded;
            }
        }

        centroids.swap(new_centroids);
        if (total_shift <= tol) break;
    }

    std::vector<long long> cluster_sizes(k, 0);
    for (int ci : assignment) if (ci >= 0 && ci < k) cluster_sizes[ci] += 1;
    return sortResults(centroids, assignment, cluster_sizes);
}

ClusteringResult kMeansEuclideanQuantizedReferences(
    const PreprocessorResult& preprocessed,
    int k,
    int max_iter,
    long double tol,
    int dim,
    KMeansReferenceQuantization quantization
) {
    if (!preprocessed.data || preprocessed.data->num_rows == 0) return {{}, {}, {}};
    const auto& data = *preprocessed.data;
    const int n = data.num_rows;
    if (dim <= 0 || dim != data.dim) {
        throw std::invalid_argument("KMeans quantization received an invalid dimension.");
    }
    if (preprocessed.max_decimal_places.size() != static_cast<size_t>(dim) ||
        preprocessed.min_values.size() != static_cast<size_t>(dim)) {
        throw std::invalid_argument("KMeans quantization metadata does not match the input dimension.");
    }

    auto scaled_points = convert_SoA_to_AoS(data);
    auto distinct_points = scaled_points;
    std::sort(distinct_points.begin(), distinct_points.end());
    distinct_points.erase(std::unique(distinct_points.begin(), distinct_points.end()), distinct_points.end());
    k = std::max(1, std::min(k, static_cast<int>(distinct_points.size())));

    std::vector<long double> scales(dim, 1.0L);
    for (int d = 0; d < dim; ++d) {
        for (int p = 0; p < preprocessed.max_decimal_places[d]; ++p) scales[d] *= 10.0L;
    }

    std::vector<std::vector<long double>> points(
        n, std::vector<long double>(dim, 0.0L)
    );
    for (int i = 0; i < n; ++i) {
        for (int d = 0; d < dim; ++d) {
            points[i][d] = (
                static_cast<long double>(scaled_points[i][d]) +
                static_cast<long double>(preprocessed.min_values[d])
            ) / scales[d];
        }
    }

    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::mt19937 generator(42);
    std::shuffle(indices.begin(), indices.end(), generator);

    std::vector<std::vector<long double>> centroids;
    centroids.reserve(k);
    std::unordered_set<std::vector<long long>, VectorHasher> seen_centroids;
    for (int index : indices) {
        if (seen_centroids.insert(scaled_points[index]).second) {
            centroids.push_back(points[index]);
            if (static_cast<int>(centroids.size()) == k) break;
        }
    }
    k = static_cast<int>(centroids.size());
    if (k == 0) return {{}, {}, {}};

    std::vector<int> assignments(n, 0);
    const auto assign_points = [&]() {
        for (int i = 0; i < n; ++i) {
            long double best_distance = std::numeric_limits<long double>::infinity();
            int best_index = 0;
            for (int ci = 0; ci < k; ++ci) {
                long double squared_distance = 0.0L;
                for (int d = 0; d < dim; ++d) {
                    const long double difference = points[i][d] - centroids[ci][d];
                    squared_distance += difference * difference;
                }
                if (squared_distance < best_distance) {
                    best_distance = squared_distance;
                    best_index = ci;
                }
            }
            assignments[i] = best_index;
        }
    };

    for (int iteration = 0; iteration < max_iter; ++iteration) {
        assign_points();
        std::vector<std::vector<long double>> sums(
            k, std::vector<long double>(dim, 0.0L)
        );
        std::vector<int> counts(k, 0);
        for (int i = 0; i < n; ++i) {
            const int ci = assignments[i];
            ++counts[ci];
            for (int d = 0; d < dim; ++d) sums[ci][d] += points[i][d];
        }

        long double max_shift = 0.0L;
        auto updated_centroids = centroids;
        for (int ci = 0; ci < k; ++ci) {
            if (counts[ci] == 0) continue;
            for (int d = 0; d < dim; ++d) {
                const long double updated = sums[ci][d] / static_cast<long double>(counts[ci]);
                max_shift = std::max(max_shift, std::abs(updated - centroids[ci][d]));
                updated_centroids[ci][d] = updated;
            }
        }
        centroids.swap(updated_centroids);
        if (max_shift <= tol) break;
    }
    assign_points();

    std::vector<std::vector<long long>> quantized_references(
        k, std::vector<long long>(dim, 0)
    );
    for (int ci = 0; ci < k; ++ci) {
        for (int d = 0; d < dim; ++d) {
            const long double quantized_value =
                quantization == KMeansReferenceQuantization::RoundTwoDecimals
                ? std::round(centroids[ci][d] * 100.0L) / 100.0L
                : std::trunc(centroids[ci][d]);
            const long double scaled_absolute = std::round(quantized_value * scales[d]);
            const long double normalized =
                scaled_absolute - static_cast<long double>(preprocessed.min_values[d]);
            if (normalized < static_cast<long double>(std::numeric_limits<long long>::min()) ||
                normalized > static_cast<long double>(std::numeric_limits<long long>::max())) {
                throw std::out_of_range("Quantized KMeans reference exceeds int64 range.");
            }
            quantized_references[ci][d] = static_cast<long long>(normalized);
        }
    }

    std::vector<long long> cluster_sizes(k, 0);
    for (int assignment : assignments) ++cluster_sizes[assignment];
    return sortResults(quantized_references, assignments, cluster_sizes);
}

ClusteringResult adaptiveGreedyBasePointSelection(const SoAData<long long>& page_data, int dim, int adjacent_reference) {
    auto data_aos = convert_SoA_to_AoS(page_data);
    size_t n = data_aos.size();
    if (n == 0) return {{}, {}, {}};
    adjacent_reference = std::max(1, adjacent_reference);
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
        std::vector<std::pair<long long, int>> cost_rank;
        cost_rank.reserve(medoids_list.size());
        int best_leader_index = -1;
        long long min_cost_to_existing_leader = std::numeric_limits<long long>::max();
        for (size_t j = 0; j < medoids_list.size(); ++j) {
            long long cost = calculateTotalLogCost(current_point, medoids_list[j], dim);
            cost_rank.push_back({cost, static_cast<int>(j)});
            if (cost < min_cost_to_existing_leader) { min_cost_to_existing_leader = cost; best_leader_index = static_cast<int>(j); }
        }
        std::sort(cost_rank.begin(), cost_rank.end(), [](const auto& a, const auto& b) {
            return a.first != b.first ? a.first < b.first : a.second < b.second;
        });
        int considered_clusters = std::min(static_cast<int>(cost_rank.size()), adjacent_reference);
        std::vector<int> candidate_cluster_indices;
        candidate_cluster_indices.reserve(considered_clusters);
        for (int ci = 0; ci < considered_clusters; ++ci) {
            candidate_cluster_indices.push_back(cost_rank[ci].second);
        }
        long long savings_from_current_point = min_cost_to_existing_leader;
        long long savings_from_reassignment = 0;
        for (int cluster_idx : candidate_cluster_indices) {
            const auto& points_in_cluster = points_in_cluster_list[cluster_idx];
            for (int point_index_in_cluster : points_in_cluster) {
                const auto& p = data_aos[point_index_in_cluster];
                long long cost_to_old_leader = calculateTotalLogCost(p, medoids_list[cluster_idx], dim);
                long long cost_to_new_potential_leader = calculateTotalLogCost(p, current_point, dim);
                if (cost_to_new_potential_leader < cost_to_old_leader) { savings_from_reassignment += (cost_to_old_leader - cost_to_new_potential_leader); }
            }
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
            for (int cluster_idx : candidate_cluster_indices) {
                std::vector<int> points_to_reassign;
                for (int point_idx : points_in_cluster_list[cluster_idx]) {
                    const auto& p = data_aos[point_idx];
                    long long cost_to_old = calculateTotalLogCost(p, medoids_list[cluster_idx], dim);
                    long long cost_to_new = calculateTotalLogCost(p, current_point, dim);
                    if (cost_to_new < cost_to_old) { points_to_reassign.push_back(point_idx); }
                }
                for (int point_idx : points_to_reassign) {
                    points_in_cluster_list[cluster_idx].erase(point_idx);
                    points_in_cluster_list[new_leader_id].insert(point_idx);
                    point_to_leader_map[point_idx] = new_leader_id;
                }
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

namespace {
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
                if (stream_debug_enabled() && i == 0 && (j < 2 || j % 100 == 0 || j + 1 == dim)) {
                    std::cerr << "[STREAM|WRITE] page=0 min_value dim=" << j
                              << " bit_pos=" << page_metadata_bitstream.currentBitPos
                              << " decimals=" << page.max_decimal_placements[j]
                              << " bit_len=" << bit_len
                              << " sign=" << (val >= 0 ? 0 : 1)
                              << " abs=" << abs_val << "\n";
                }
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
                if (stream_debug_enabled() && i == 0 && j < 2) {
                    std::cerr << "[STREAM|WRITE] page=0 min_base dim=" << j
                              << " bit_pos=" << page_metadata_bitstream.currentBitPos
                              << " bit_len=" << bit_len
                              << " sign=" << (val >= 0 ? 0 : 1)
                              << " abs=" << abs_val
                              << " base_width=" << page.max_base_bit_len[j] << "\n";
                }
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

ClusterEncoderPageResult clusterEncoderPageFromClustering(
    std::unique_ptr<PreprocessorResult> preproc_res_ptr,
    ClusteringResult clustering_res,
    int pack_size, int dim, int block_size, int page_idx)
{
    if (!preproc_res_ptr) throw std::runtime_error("Preprocessor result is null.");
    auto& preproc_res = *preproc_res_ptr;

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

    int page_k = clustering_res.medoids.size();
    const long long encoded_point_count = std::accumulate(
        clustering_res.cluster_sizes.begin(), clustering_res.cluster_sizes.end(), 0LL
    );
    if (encoded_point_count != preproc_res_ptr->data->num_rows) {
        throw std::runtime_error(
            "Cluster-size invariant failed while encoding page " + std::to_string(page_idx) +
            ": sizes sum to " + std::to_string(encoded_point_count) +
            ", expected " + std::to_string(preproc_res_ptr->data->num_rows)
        );
    }
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

ClusterEncoderPageResult clusterEncoderPage(
    std::unique_ptr<PreprocessorResult> preproc_res_ptr,
    const std::function<ClusteringResult(const SoAData<long long>&, int, int)>& cluster_strategy,
    int k, int pack_size, int dim, int block_size, int page_idx)
{
    if (!preproc_res_ptr || !preproc_res_ptr->data) {
        throw std::runtime_error("Preprocessor result is null.");
    }
    ClusteringResult clustering_res = cluster_strategy(*preproc_res_ptr->data, k, dim);
    return clusterEncoderPageFromClustering(
        std::move(preproc_res_ptr), std::move(clustering_res),
        pack_size, dim, block_size, page_idx
    );
}

std::pair<ClusterEncoderPageResult, PageDiagnostics> clusterEncoderPage_with_diag(
    std::unique_ptr<PreprocessorResult> preproc_res_ptr,
    const std::string& algorithm_name,
    int k, int pack_size, int dim, int block_size, int page_idx, int adjacent_reference,
    int kcluster_max_iterations)
{
    if (!preproc_res_ptr) throw std::runtime_error("Preprocessor result is null.");
    auto& preproc_res = *preproc_res_ptr;

    PageDiagnostics diag;
    diag.page_index = page_idx;
    diag.algorithm_name = algorithm_name;
    auto& timings = diag.timings;

    auto start_cluster = std::chrono::high_resolution_clock::now();
    ClusteringResult clustering_res;
    if (algorithm_name == "KCluster") {
        auto k_res = kMedoidLogCost_with_diag(
            *preproc_res.data, k, kcluster_max_iterations, 1.0e-3, dim
        );
        clustering_res = k_res.result;
        diag.kcluster_costs.final_residual_cost_bits = k_res.final_cost;
    } else if (algorithm_name == "ACluster") {
        auto a_res = adaptiveGreedyBasePointSelection_with_diag(*preproc_res.data, dim, adjacent_reference);
        clustering_res = a_res.result;
        diag.acluster_costs.final_total_cost_bits = a_res.final_cost;
    } else {

        clustering_res = kMedoidLogCost(
            *preproc_res.data, k, kcluster_max_iterations, 1.0e-3, dim
        );
    }
    auto end_cluster = std::chrono::high_resolution_clock::now();
    timings.clustering_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_cluster - start_cluster).count();

    auto start_residual = std::chrono::high_resolution_clock::now();
    std::vector<std::vector<long long>> residuals = residualCalculationZigzag_sorted(
        *preproc_res.data, clustering_res.medoids, clustering_res.cluster_assignment,
        clustering_res.cluster_sizes, dim);
    auto end_residual = std::chrono::high_resolution_clock::now();
    timings.residual_calc_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_residual - start_residual).count();

    auto start_bitstream = std::chrono::high_resolution_clock::now();

    int page_k = clustering_res.medoids.size();
    const long long encoded_point_count = std::accumulate(
        clustering_res.cluster_sizes.begin(), clustering_res.cluster_sizes.end(), 0LL
    );
    if (encoded_point_count != preproc_res.data->num_rows) {
        throw std::runtime_error(
            "Cluster-size invariant failed while encoding diagnostic page " + std::to_string(page_idx) +
            ": sizes sum to " + std::to_string(encoded_point_count) +
            ", expected " + std::to_string(preproc_res.data->num_rows)
        );
    }
    diag.page_reference_count = page_k;
    for (const auto& medoid : clustering_res.medoids) {
        diag.raw_base_cost_bits += calculateBasePointStorageCost(medoid, dim);
    }
    for (int i = 0; i < preproc_res.data->num_rows; ++i) {
        int cid = clustering_res.cluster_assignment[i];
        if (cid < 0 || cid >= page_k) continue;
        diag.raw_residual_cost_bits += calculateTotalLogCost(get_point_from_soa(*preproc_res.data, i), clustering_res.medoids[cid], dim);
    }
    BitBuffer fre_bitstream = generateBitstreamFrequency(clustering_res.cluster_sizes, block_size);
    EncodedDataResult bitstream_res = generateBitstreamEncodedData(residuals, pack_size, dim);
    BasePointsResult basepoint_res = generateBitstreamEncodedBasePointsOptimized(clustering_res.medoids, dim);

    auto end_bitstream = std::chrono::high_resolution_clock::now();
    timings.bitstream_gen_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_bitstream - start_bitstream).count();

    diag.sizes.frequencies_bits = fre_bitstream.size();
    diag.sizes.residuals_bits = bitstream_res.residual_bitstream.size();
    diag.sizes.medoids_bits = basepoint_res.base_bitstream.size();

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
    int k, int pack_size, int block_size, int page_size, int adjacent_reference,
    int kcluster_max_iterations)
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
            std::move(preproc_res_ptr), algorithm_name, k, pack_size, current_dim, block_size,
            page_index_counter, adjacent_reference, kcluster_max_iterations);

        page_diag.timings.preprocessing_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end_preproc - start_preproc).count();

        global_diag.total_timings.preprocessing_ns += page_diag.timings.preprocessing_ns;
        global_diag.total_timings.clustering_ns += page_diag.timings.clustering_ns;
        global_diag.total_timings.residual_calc_ns += page_diag.timings.residual_calc_ns;
        global_diag.total_timings.bitstream_gen_ns += page_diag.timings.bitstream_gen_ns;

        global_diag.total_sizes.medoids_bits += page_diag.sizes.medoids_bits;
        global_diag.total_sizes.frequencies_bits += page_diag.sizes.frequencies_bits;
        global_diag.total_sizes.residuals_bits += page_diag.sizes.residuals_bits;
        global_diag.total_reference_points += page_diag.page_reference_count;
        global_diag.total_raw_base_cost_bits += page_diag.raw_base_cost_bits;
        global_diag.total_raw_residual_cost_bits += page_diag.raw_residual_cost_bits;

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

            page_header_buf.currentBitPos += 8;
            page_header_buf.currentBitPos += 8 + 1 + 32;
            page_header_buf.currentBitPos += 8 + 1 + 32;
            page_header_buf.currentBitPos += 8;
        }
        page_header_buf.currentBitPos += 16;
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

std::vector<uint8_t> encode_multidim_kmeans_quantized_impl(
    const std::string& csv_file_path,
    int& dim,
    int k,
    int pack_size,
    int block_size,
    int page_size,
    KMeansReferenceQuantization quantization
) {
    std::vector<ClusterEncoderPageResult> all_page_results;
    std::ifstream file_stream(csv_file_path);
    if (!file_stream.is_open()) throw std::runtime_error("Cannot open data file: " + csv_file_path);
    std::vector<std::string> page_lines;
    page_lines.reserve(page_size);
    std::string line;
    bool file_is_finished = false;
    while (!file_is_finished) {
        page_lines.clear();
        while (page_lines.size() < static_cast<size_t>(page_size) && std::getline(file_stream, line)) {
            if (!line.empty()) page_lines.push_back(line);
        }
        if (file_stream.eof()) file_is_finished = true;
        if (page_lines.empty()) break;

        auto preprocessed = std::make_unique<PreprocessorResult>(
            preprocessPageFromStrings(page_lines, dim)
        );
        if (!preprocessed || !preprocessed->data || preprocessed->data->num_rows == 0) continue;
        const int current_dim = preprocessed->data->dim;
        const int page_index = static_cast<int>(all_page_results.size());
        ClusteringResult clustering = kMeansEuclideanQuantizedReferences(
            *preprocessed, k, 20, 1.0e-9L, current_dim, quantization
        );
        all_page_results.push_back(clusterEncoderPageFromClustering(
            std::move(preprocessed), std::move(clustering),
            pack_size, current_dim, block_size, page_index
        ));
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
    const auto require_bit_width = [](int bits, const std::string& context) {
        if (bits < 1 || bits > 64) {
            throw std::runtime_error("Invalid bit width for " + context + ": " + std::to_string(bits));
        }
    };

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
            require_bit_width(bit_len, "min value at page " + std::to_string(i) + ", dimension " + std::to_string(j));
            long long sign = reader.readBits(1);
            long long abs_val = reader.readBits(bit_len);
            min_values[j] = (sign == 0) ? abs_val : -abs_val;
            if (stream_debug_enabled() && i == 0 && (j < 2 || j % 100 == 0 || j + 1 == dim)) {
                std::cerr << "[STREAM|READ] page=0 min_value dim=" << j
                          << " bit_pos=" << reader.position()
                          << " decimals=" << max_decimals[j]
                          << " bit_len=" << bit_len
                          << " sign=" << sign
                          << " abs=" << abs_val << "\n";
            }
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": Decimals=" << max_decimals[j] << ", MinValLen=" << bit_len << ", Sign=" << sign << ", AbsMinVal=" << abs_val << "\n";
        }

        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  --- Reading Min Bases & Base Bit Lengths ---\n";
        std::vector<long long> min_bases(dim);
        std::vector<int> max_base_bit_len(dim);
        for (int j = 0; j < dim; ++j) {
            int bit_len = reader.readBits(8);
            require_bit_width(bit_len, "min base at page " + std::to_string(i) + ", dimension " + std::to_string(j));
            long long sign = reader.readBits(1);
            long long abs_val = reader.readBits(bit_len);
            min_bases[j] = (sign == 0) ? abs_val : -abs_val;
            max_base_bit_len[j] = reader.readBits(8);
            require_bit_width(max_base_bit_len[j], "base offset at page " + std::to_string(i) + ", dimension " + std::to_string(j));
            if (stream_debug_enabled() && i == 0 && j < 2) {
                std::cerr << "[STREAM|READ] page=0 min_base dim=" << j
                          << " bit_pos=" << reader.position()
                          << " bit_len=" << bit_len
                          << " sign=" << sign
                          << " abs=" << abs_val
                          << " base_width=" << max_base_bit_len[j] << "\n";
            }
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "    - Dim " << j << ": MinBaseLen=" << bit_len << ", Sign=" << sign << ", AbsMinBase=" << abs_val << ", MaxBaseBitLen=" << max_base_bit_len[j] << "\n";
        }

        int pack_num = reader.readBits(16);
        if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) std::cout << "  - Read Pack Metadata Count: " << pack_num << "\n";
        std::vector<std::vector<int>> pack_metadata(pack_num, std::vector<int>(dim));
        for (int p = 0; p < pack_num; ++p) {
            for (int j = 0; j < dim; ++j) {
                pack_metadata[p][j] = reader.readBits(8);
                if (pack_metadata[p][j] < 1 || pack_metadata[p][j] > 64) {
                    throw std::runtime_error(
                        "Invalid residual bit width in page " + std::to_string(i) +
                        ", pack " + std::to_string(p) + ", dimension " + std::to_string(j) +
                        ": " + std::to_string(pack_metadata[p][j])
                    );
                }
            }
        }

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

            std::vector<int> max_bits_per_block(num_freq_blocks);
            if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
                std::cout << "\n--- [DEBUG|FREQ_READ] Reading all " << num_freq_blocks << " max_bit metadata blocks ---\n";
            }
            for (int b = 0; b < num_freq_blocks; ++b) {
                max_bits_per_block[b] = reader.readBits(8);
                require_bit_width(max_bits_per_block[b], "cluster frequency at page " + std::to_string(i) + ", block " + std::to_string(b));
                if (GlobalDebug::DEBUG_PRINT_ENABLED && debug_page) {
                    std::cout << "  - Block " << b << ": Will read with max_bit=" << max_bits_per_block[b] << "\n";
                }
            }

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

            cluster_sizes[0] = deltas[0];
            for (size_t j = 1; j < deltas.size(); ++j) {
                cluster_sizes[j] = cluster_sizes[j - 1] + deltas[j];
            }
            const long long decoded_point_count = std::accumulate(cluster_sizes.begin(), cluster_sizes.end(), 0LL);
            if (decoded_point_count != page_data_points) {
                throw std::runtime_error(
                    "Cluster-frequency invariant failed while decoding page " + std::to_string(i) +
                    ": sizes sum to " + std::to_string(decoded_point_count) +
                    ", expected " + std::to_string(page_data_points)
                );
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

ClusterPageEncoding encodePreprocessedPage(
    PreprocessorResult preprocessed, const std::string& method,
    int k, int page_size, int pack_size, int block_size, uint32_t seed,
    bool measure_objective)
{
    if (!preprocessed.data || preprocessed.data->num_rows <= 0 ||
        preprocessed.data->dim <= 0 || preprocessed.data->dim > 255 ||
        preprocessed.data->num_rows > page_size || page_size > 10000 ||
        k <= 0 || pack_size <= 0 || pack_size > 65535 || block_size <= 0 || block_size > 65535) {
        throw std::invalid_argument("Invalid in-memory cluster page.");
    }
    const int dim = preprocessed.data->dim;
    const int n = preprocessed.data->num_rows;
    ClusteringResult clusters;
    ClusterPageEncoding result;
    using StageClock = std::chrono::steady_clock;
    const auto selection_start = StageClock::now();
    if (method == "ACluster1D") {
        clusters = adaptiveGreedy1D(*preprocessed.data).result;
    } else if (method == "KCluster1D") {
        clusters = kMedoid1D(*preprocessed.data,k,10,1e-3,seed).result;
    } else if (method == "ACluster") {
        clusters = adaptiveGreedyBasePointSelection_with_diag(*preprocessed.data, dim, 1).result;
    } else if (method == "KCluster") {
        clusters = kMedoidLogCost_with_diag(*preprocessed.data, k, 10, 1e-3, dim, &seed).result;
    } else if (method == "KOptimal" || method == "AOptimal") {
        auto optimal = optimal1D(*preprocessed.data, method == "KOptimal" ? k : 0);
        result.optimal_dp_states = optimal.dp_states;
        result.optimal_reference_candidates = optimal.reference_candidates;
        result.optimal_sort_ns = optimal.sort_ns;
        result.optimal_prefix_ns = optimal.prefix_ns;
        result.optimal_search_ns = optimal.search_ns;
        result.optimal_dp_ns = optimal.dp_ns;
        result.residual_objective_bits = optimal.residual_cost;
        result.reference_objective_bits = optimal.reference_cost;
        clusters = std::move(optimal.result);
    } else {
        throw std::invalid_argument("Unknown cluster page method: " + method);
    }
    const auto selection_end = StageClock::now();
    result.selection_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(selection_end-selection_start).count();
    if (measure_objective && result.residual_objective_bits < 0) {
        result.residual_objective_bits = 0;
        result.reference_objective_bits = 0;
        for (const auto& ref : clusters.medoids)
            for (auto value : ref) result.reference_objective_bits += paperScalarCost(value);
        for (int i = 0; i < n; ++i)
            for (int d = 0; d < dim; ++d)
                result.residual_objective_bits += paperScalarCost(
                    preprocessed.data->columns[d][i], clusters.medoids[clusters.cluster_assignment[i]][d]);
    }
    const auto diagnostics_end = StageClock::now();
    result.diagnostics_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(diagnostics_end-selection_end).count();
    result.reference_count = clusters.medoids.size();
    std::vector<size_t> offsets(clusters.cluster_sizes.size());
    size_t total = 0;
    for (size_t c = 0; c < offsets.size(); ++c) {
        offsets[c] = total;
        total += clusters.cluster_sizes[c];
    }
    if (total != static_cast<size_t>(n) || clusters.cluster_assignment.size() != total)
        throw std::runtime_error("Invalid cluster page assignment size.");
    result.row_order.resize(n);
    std::vector<size_t> counts(offsets.size(), 0);
    for (int r = 0; r < n; ++r) {
        const int c = clusters.cluster_assignment[r];
        if (c < 0 || static_cast<size_t>(c) >= offsets.size() ||
            counts[c]++ >= static_cast<size_t>(clusters.cluster_sizes[c]))
            throw std::runtime_error("Invalid cluster page assignment.");
        result.row_order[offsets[c]++] = r;
    }
    std::vector<ClusterEncoderPageResult> pages;
    pages.push_back(clusterEncoderPageFromClustering(
        std::make_unique<PreprocessorResult>(std::move(preprocessed)), std::move(clusters),
        pack_size, dim, block_size, 0));
    result.bytes = assemble_final_bytestream(pages, dim, pack_size, block_size, page_size);
    result.sizes.medoids_bits = pages.front().base_bitstream.size();
    result.sizes.frequencies_bits = pages.front().fre_bitstream.size();
    result.sizes.residuals_bits = pages.front().residual_bitstream.size();

    result.sizes.header_bits = result.bytes.size() * 8 - result.sizes.medoids_bits
        - result.sizes.frequencies_bits - result.sizes.residuals_bits;
    result.serialization_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(StageClock::now()-diagnostics_end).count();
    return result;
}

ReferenceCandidateEncoding encodePageWithReferences(
    const PreprocessorResult& preprocessed,
    std::vector<std::vector<long long>> references,
    int page_size, int pack_size, int block_size)
{
    if (!preprocessed.data || preprocessed.data->num_rows <= 0 || references.empty()) {
        throw std::invalid_argument("Reference candidate requires a nonempty page and references.");
    }
    const auto& data = *preprocessed.data;
    if (page_size < data.num_rows || page_size > 10000 || pack_size <= 0 ||
        pack_size > 65535 || block_size <= 0 || block_size > 65535) {
        throw std::invalid_argument("Invalid reference candidate page/packing parameters.");
    }
    auto rows = convert_SoA_to_AoS(data);
    std::set<std::vector<long long>> observed(rows.begin(), rows.end());
    std::sort(references.begin(), references.end());
    if (std::adjacent_find(references.begin(), references.end()) != references.end()) {
        throw std::invalid_argument("Reference candidates must be distinct observed tuples.");
    }
    for (const auto& reference : references) {
        if (observed.count(reference) == 0) {
            throw std::invalid_argument("Reference candidate is not an observed tuple.");
        }
    }

    auto paper_component_cost = [](long long v) -> int {
        return v == 0 ? 1 : bitLength(v) + 1;
    };
    ReferenceCandidateEncoding result;
    result.reference_count = static_cast<int>(references.size());
    for (const auto& reference : references) {
        for (long long value : reference) {
            result.reference_objective_bits += paper_component_cost(value);
        }
    }
    std::vector<int> assignments(data.num_rows);
    std::vector<long long> sizes(references.size(), 0);
    for (int i = 0; i < data.num_rows; ++i) {
        long long best_cost = std::numeric_limits<long long>::max();
        int best = -1;
        for (int j = 0; j < result.reference_count; ++j) {
            long long cost = 0;
            for (int m = 0; m < data.dim; ++m) {
                cost += paper_component_cost(rows[i][m] - references[j][m]);
                if (cost >= best_cost) break;
            }

            if (cost < best_cost) { best_cost = cost; best = j; }
        }
        assignments[i] = best;
        ++sizes[best];
        result.residual_objective_bits += best_cost;
    }
    auto clustering = sortResults(references, assignments, sizes);
    auto copy = std::make_unique<PreprocessorResult>();
    copy->data = std::make_unique<SoAData<long long>>(data);
    copy->max_decimal_places = preprocessed.max_decimal_places;
    copy->min_values = preprocessed.min_values;
    std::vector<ClusterEncoderPageResult> pages;
    pages.push_back(clusterEncoderPageFromClustering(
        std::move(copy), std::move(clustering), pack_size, data.dim, block_size, 0));
    result.bytes = assemble_final_bytestream(pages, data.dim, pack_size, block_size, page_size);
    return result;
}

}
