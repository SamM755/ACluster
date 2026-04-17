#include "acluster_sub.h"
#include "acluster.h"
#include "../utils/csv_reader.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <limits>

namespace {

constexpr double kQuant = 1e6;

void append_u32(std::vector<uint8_t>& out, uint32_t v) {
    out.push_back(static_cast<uint8_t>(v & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 8) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 16) & 0xFF));
    out.push_back(static_cast<uint8_t>((v >> 24) & 0xFF));
}

void append_u64(std::vector<uint8_t>& out, uint64_t v) {
    append_u32(out, static_cast<uint32_t>(v & 0xFFFFFFFFULL));
    append_u32(out, static_cast<uint32_t>((v >> 32) & 0xFFFFFFFFULL));
}

uint32_t read_u32(const std::vector<uint8_t>& in, size_t& pos) {
    if (pos + 4 > in.size()) throw std::runtime_error("AClusterSub invalid stream");
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

void append_uvarint(std::vector<uint8_t>& out, uint64_t v) {
    while (v >= 0x80ULL) {
        out.push_back(static_cast<uint8_t>((v & 0x7FULL) | 0x80ULL));
        v >>= 7;
    }
    out.push_back(static_cast<uint8_t>(v));
}

uint64_t read_uvarint(const std::vector<uint8_t>& in, size_t& pos) {
    uint64_t v = 0;
    int shift = 0;
    while (true) {
        if (pos >= in.size() || shift > 63) throw std::runtime_error("AClusterSub invalid varint");
        uint8_t b = in[pos++];
        v |= static_cast<uint64_t>(b & 0x7F) << shift;
        if ((b & 0x80) == 0) break;
        shift += 7;
    }
    return v;
}

void append_svarint(std::vector<uint8_t>& out, int64_t v) {
    uint64_t zz = (static_cast<uint64_t>(v) << 1) ^ static_cast<uint64_t>(v >> 63);
    append_uvarint(out, zz);
}

int64_t read_svarint(const std::vector<uint8_t>& in, size_t& pos) {
    uint64_t zz = read_uvarint(in, pos);
    return static_cast<int64_t>((zz >> 1) ^ (~(zz & 1) + 1));
}

int bits_required(uint32_t v) {
    int bits = 0;
    do { bits++; v >>= 1; } while (v);
    return bits;
}

std::vector<uint8_t> pack_bits(const std::vector<uint32_t>& vals, int bw) {
    std::vector<uint8_t> out;
    if (vals.empty() || bw <= 0) return out;
    uint64_t acc = 0;
    int acc_bits = 0;
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1ULL);
    for (uint32_t v : vals) {
        acc |= (static_cast<uint64_t>(v) & mask) << acc_bits;
        acc_bits += bw;
        while (acc_bits >= 8) {
            out.push_back(static_cast<uint8_t>(acc & 0xFF));
            acc >>= 8;
            acc_bits -= 8;
        }
    }
    if (acc_bits > 0) out.push_back(static_cast<uint8_t>(acc & 0xFF));
    return out;
}

std::vector<uint32_t> unpack_bits(const std::vector<uint8_t>& bytes, size_t expected, int bw) {
    std::vector<uint32_t> out;
    out.reserve(expected);
    uint64_t acc = 0;
    int acc_bits = 0;
    size_t p = 0;
    uint64_t mask = (bw == 64) ? ~0ULL : ((1ULL << bw) - 1ULL);
    while (out.size() < expected) {
        while (acc_bits < bw) {
            if (p >= bytes.size()) throw std::runtime_error("AClusterSub unpack underflow");
            acc |= static_cast<uint64_t>(bytes[p++]) << acc_bits;
            acc_bits += 8;
        }
        out.push_back(static_cast<uint32_t>(acc & mask));
        acc >>= bw;
        acc_bits -= bw;
    }
    return out;
}

std::vector<uint8_t> encode_ids(const std::vector<uint32_t>& ids, uint32_t block_size) {
    std::vector<uint8_t> out;
    append_u32(out, static_cast<uint32_t>(ids.size()));
    append_u32(out, block_size);
    size_t p = 0;
    while (p < ids.size()) {
        size_t cnt = std::min<size_t>(block_size, ids.size() - p);
        uint32_t mn = ids[p], mx = ids[p];
        for (size_t i = 1; i < cnt; ++i) {
            mn = std::min(mn, ids[p + i]);
            mx = std::max(mx, ids[p + i]);
        }
        std::vector<uint32_t> deltas(cnt);
        for (size_t i = 0; i < cnt; ++i) deltas[i] = ids[p + i] - mn;
        int bw = bits_required(mx - mn);
        auto packed = pack_bits(deltas, bw);
        append_u32(out, static_cast<uint32_t>(cnt));
        append_u32(out, mn);
        out.push_back(static_cast<uint8_t>(bw));
        append_u32(out, static_cast<uint32_t>(packed.size()));
        out.insert(out.end(), packed.begin(), packed.end());
        p += cnt;
    }
    return out;
}

std::vector<uint32_t> decode_ids(const std::vector<uint8_t>& in) {
    size_t pos = 0;
    uint32_t total = read_u32(in, pos);
    uint32_t block_size = read_u32(in, pos);
    std::vector<uint32_t> ids;
    ids.reserve(total);
    while (ids.size() < total) {
        uint32_t cnt = read_u32(in, pos);
        uint32_t mn = read_u32(in, pos);
        int bw = static_cast<int>(in[pos++]);
        uint32_t packed_sz = read_u32(in, pos);
        if (pos + packed_sz > in.size()) throw std::runtime_error("AClusterSub decode ids truncated");
        std::vector<uint8_t> packed(in.begin() + pos, in.begin() + pos + packed_sz);
        pos += packed_sz;
        auto deltas = unpack_bits(packed, cnt, bw);
        for (uint32_t d : deltas) ids.push_back(mn + d);
        if (cnt > block_size) throw std::runtime_error("AClusterSub decode ids bad block");
    }
    return ids;
}

uint32_t zigzag_encode_32(int64_t v) {
    return static_cast<uint32_t>((v << 1) ^ (v >> 63));
}

int64_t zigzag_decode_32(uint32_t v) {
    return static_cast<int64_t>((v >> 1) ^ static_cast<uint32_t>(-static_cast<int32_t>(v & 1)));
}

std::vector<uint8_t> encode_positions(const std::vector<uint32_t>& pos) {
    auto abs_stream = encode_ids(pos, 256);
    std::vector<uint32_t> deltas;
    deltas.reserve(pos.size());
    for (size_t i = 0; i < pos.size(); ++i) {
        int64_t d = static_cast<int64_t>(pos[i]) - static_cast<int64_t>(i);
        deltas.push_back(zigzag_encode_32(d));
    }
    auto delta_stream = encode_ids(deltas, 256);
    std::vector<uint8_t> out;
    if (delta_stream.size() < abs_stream.size()) {
        out.push_back(1);
        out.insert(out.end(), delta_stream.begin(), delta_stream.end());
    } else {
        out.push_back(0);
        out.insert(out.end(), abs_stream.begin(), abs_stream.end());
    }
    return out;
}

std::vector<uint32_t> decode_positions(const std::vector<uint8_t>& in, uint32_t row_count) {
    if (in.empty()) throw std::runtime_error("AClusterSub decode positions empty");
    uint8_t mode = in[0];
    std::vector<uint8_t> body(in.begin() + 1, in.end());
    if (mode == 0) {
        auto pos = decode_ids(body);
        for (uint32_t p : pos) if (p >= row_count) throw std::runtime_error("AClusterSub position out of range");
        return pos;
    }
    if (mode == 1) {
        auto deltas = decode_ids(body);
        std::vector<uint32_t> pos;
        pos.reserve(deltas.size());
        for (size_t i = 0; i < deltas.size(); ++i) {
            int64_t d = zigzag_decode_32(deltas[i]);
            int64_t p = static_cast<int64_t>(i) + d;
            if (p < 0 || p >= static_cast<int64_t>(row_count)) throw std::runtime_error("AClusterSub delta position out of range");
            pos.push_back(static_cast<uint32_t>(p));
        }
        return pos;
    }
    throw std::runtime_error("AClusterSub decode positions mode unknown");
}

long long quant_key(double v) {
    return static_cast<long long>(std::llround(v * kQuant));
}

std::vector<uint8_t> encode_ambiguous_positions(
    const std::vector<std::vector<double>>& dec0,
    int anchor_idx_in_g0,
    const std::vector<std::vector<double>>& dec1,
    int anchor_idx_in_g1plus,
    const std::vector<uint32_t>& pos1
) {
    std::unordered_map<long long, std::vector<uint32_t>> bucket;
    bucket.reserve(dec0.size() * 2);
    for (uint32_t i = 0; i < static_cast<uint32_t>(dec0.size()); ++i) bucket[quant_key(dec0[i][anchor_idx_in_g0])].push_back(i);
    std::vector<uint32_t> local_idx;
    for (size_t i = 0; i < dec1.size(); ++i) {
        long long k = quant_key(dec1[i][anchor_idx_in_g1plus]);
        const auto& b = bucket[k];
        if (b.size() <= 1) continue;
        auto it = std::find(b.begin(), b.end(), pos1[i]);
        if (it == b.end()) throw std::runtime_error("AClusterSub ambiguous encode mapping miss");
        local_idx.push_back(static_cast<uint32_t>(std::distance(b.begin(), it)));
    }
    auto idx_stream = encode_ids(local_idx, 256);
    std::vector<uint8_t> out;
    append_u32(out, static_cast<uint32_t>(idx_stream.size()));
    out.insert(out.end(), idx_stream.begin(), idx_stream.end());
    return out;
}

std::vector<uint32_t> decode_ambiguous_positions(
    const std::vector<std::vector<double>>& dec0,
    int anchor_idx_in_g0,
    const std::vector<std::vector<double>>& dec1,
    int anchor_idx_in_g1plus,
    const std::vector<uint8_t>& ambiguous_stream,
    uint32_t row_count
) {
    size_t p = 0;
    uint32_t idx_sz = read_u32(ambiguous_stream, p);
    if (p + idx_sz > ambiguous_stream.size()) throw std::runtime_error("AClusterSub ambiguous stream truncated");
    std::vector<uint8_t> idx_bytes(ambiguous_stream.begin() + p, ambiguous_stream.begin() + p + idx_sz);
    auto idx = decode_ids(idx_bytes);

    std::unordered_map<long long, std::vector<uint32_t>> bucket;
    bucket.reserve(dec0.size() * 2);
    for (uint32_t i = 0; i < static_cast<uint32_t>(dec0.size()); ++i) bucket[quant_key(dec0[i][anchor_idx_in_g0])].push_back(i);
    std::vector<uint32_t> pos1;
    pos1.reserve(dec1.size());
    size_t idx_ptr = 0;
    for (size_t i = 0; i < dec1.size(); ++i) {
        long long k = quant_key(dec1[i][anchor_idx_in_g1plus]);
        const auto& b = bucket[k];
        if (b.empty()) throw std::runtime_error("AClusterSub ambiguous decode bucket empty");
        if (b.size() == 1) {
            pos1.push_back(b[0]);
            continue;
        }
        if (idx_ptr >= idx.size()) throw std::runtime_error("AClusterSub ambiguous decode idx underflow");
        uint32_t li = idx[idx_ptr++];
        if (li >= b.size()) throw std::runtime_error("AClusterSub ambiguous decode local idx out of range");
        pos1.push_back(b[li]);
    }
    if (idx_ptr != idx.size()) throw std::runtime_error("AClusterSub ambiguous decode idx overflow");
    for (uint32_t p1 : pos1) if (p1 >= row_count) throw std::runtime_error("AClusterSub ambiguous decode pos out of range");
    return pos1;
}

double abs_corr(const std::vector<std::vector<double>>& rows, int c1, int c2) {
    if (rows.empty()) return 0.0;
    double sx = 0.0, sy = 0.0, sxx = 0.0, syy = 0.0, sxy = 0.0, n = 0.0;
    for (const auto& r : rows) {
        double x = r[c1], y = r[c2];
        if (std::isnan(x) || std::isnan(y) || std::isinf(x) || std::isinf(y)) continue;
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; n += 1.0;
    }
    if (n <= 1.0) return 0.0;
    double vx = sxx - sx * sx / n;
    double vy = syy - sy * sy / n;
    if (vx <= 1e-12 || vy <= 1e-12) return 0.0;
    return std::abs((sxy - sx * sy / n) / std::sqrt(vx * vy));
}

std::vector<std::vector<double>> sample_rows(const std::vector<std::vector<double>>& rows, size_t cap = 6000) {
    if (rows.size() <= cap) return rows;
    std::vector<std::vector<double>> s;
    s.reserve(cap);
    std::mt19937 rng(20260408);
    std::uniform_int_distribution<size_t> dis(0, rows.size() - 1);
    for (size_t i = 0; i < cap; ++i) s.push_back(rows[dis(rng)]);
    return s;
}

void build_corr_groups(const std::vector<std::vector<double>>& rows, std::vector<int>& g0, std::vector<int>& g1) {
    int dim = static_cast<int>(rows[0].size());
    g0.clear(); g1.clear();
    if (dim <= 2) {
        g0.push_back(0);
        for (int c = 1; c < dim; ++c) g1.push_back(c);
        return;
    }
    auto sample = sample_rows(rows);
    std::vector<std::vector<double>> corr(dim, std::vector<double>(dim, 0.0));
    for (int i = 0; i < dim; ++i) for (int j = i + 1; j < dim; ++j) corr[i][j] = corr[j][i] = abs_corr(sample, i, j);
    int s0 = 0, s1 = 1;
    double best_sep = -1.0;
    for (int i = 0; i < dim; ++i) for (int j = i + 1; j < dim; ++j) {
        double sep = 1.0 - corr[i][j];
        if (sep > best_sep) { best_sep = sep; s0 = i; s1 = j; }
    }
    g0.push_back(s0); g1.push_back(s1);
    std::vector<int> rem;
    for (int c = 0; c < dim; ++c) if (c != s0 && c != s1) rem.push_back(c);
    int target0 = dim / 2;
    for (int c : rem) {
        double a0 = 0.0, a1 = 0.0;
        for (int x : g0) a0 += corr[c][x];
        for (int x : g1) a1 += corr[c][x];
        a0 /= std::max(1, static_cast<int>(g0.size()));
        a1 /= std::max(1, static_cast<int>(g1.size()));
        if (static_cast<int>(g0.size()) < target0 && (a0 >= a1 || static_cast<int>(g1.size()) >= dim - target0)) g0.push_back(c);
        else g1.push_back(c);
    }
    if (g0.empty()) { g0.push_back(g1.back()); g1.pop_back(); }
    if (g1.empty()) { g1.push_back(g0.back()); g0.pop_back(); }
    std::sort(g0.begin(), g0.end());
    std::sort(g1.begin(), g1.end());
}

double column_variance(const std::vector<std::vector<double>>& rows, int c) {
    if (rows.empty()) return 0.0;
    double s = 0.0, ss = 0.0, n = 0.0;
    for (const auto& r : rows) {
        double v = r[c];
        if (std::isnan(v) || std::isinf(v)) continue;
        s += v;
        ss += v * v;
        n += 1.0;
    }
    if (n <= 1.0) return 0.0;
    return std::max(0.0, ss - s * s / n);
}

double unique_ratio(const std::vector<std::vector<double>>& rows, int c) {
    if (rows.empty()) return 0.0;
    std::unordered_map<long long, int> freq;
    freq.reserve(rows.size() * 2);
    int finite_count = 0;
    for (const auto& r : rows) {
        double v = r[c];
        if (!std::isfinite(v)) continue;
        freq[quant_key(v)] += 1;
        finite_count += 1;
    }
    if (finite_count <= 0) return 0.0;
    return static_cast<double>(freq.size()) / static_cast<double>(finite_count);
}

std::vector<std::pair<std::vector<int>, std::vector<int>>> build_candidate_splits(const std::vector<std::vector<double>>& rows) {
    int dim = static_cast<int>(rows[0].size());
    std::vector<std::pair<std::vector<int>, std::vector<int>>> cands;
    std::unordered_map<std::string, int> seen;
    auto add_cand = [&](std::vector<int> g0, std::vector<int> g1) {
        if (g0.empty() || g1.empty()) return;
        std::sort(g0.begin(), g0.end());
        std::sort(g1.begin(), g1.end());
        std::string k;
        for (int x : g0) { k += std::to_string(x); k.push_back('_'); }
        k.push_back('|');
        for (int x : g1) { k += std::to_string(x); k.push_back('_'); }
        if (seen.count(k)) return;
        seen[k] = 1;
        cands.push_back({std::move(g0), std::move(g1)});
    };

    std::vector<int> g0, g1;
    build_corr_groups(rows, g0, g1);
    add_cand(g0, g1);

    std::vector<std::pair<double, int>> by_corr, by_var;
    auto sample = sample_rows(rows);
    for (int c = 0; c < dim; ++c) {
        double ac = 0.0;
        for (int j = 0; j < dim; ++j) if (j != c) ac += abs_corr(sample, c, j);
        ac /= std::max(1, dim - 1);
        by_corr.push_back({ac, c});
        by_var.push_back({column_variance(sample, c), c});
    }
    std::sort(by_corr.begin(), by_corr.end());
    std::sort(by_var.begin(), by_var.end());

    auto add_tail = [&](const std::vector<std::pair<double, int>>& order, int tail, bool take_small) {
        if (tail <= 0 || tail >= dim) return;
        std::vector<int> tset;
        if (take_small) {
            for (int i = 0; i < tail; ++i) tset.push_back(order[i].second);
        } else {
            for (int i = 0; i < tail; ++i) tset.push_back(order[dim - 1 - i].second);
        }
        std::sort(tset.begin(), tset.end());
        std::vector<int> a, b;
        for (int c = 0; c < dim; ++c) {
            if (std::binary_search(tset.begin(), tset.end(), c)) b.push_back(c);
            else a.push_back(c);
        }
        add_cand(a, b);
    };

    add_tail(by_corr, 1, true);
    add_tail(by_corr, 2, true);
    add_tail(by_corr, 3, true);
    add_tail(by_var, 1, true);
    add_tail(by_var, 2, true);
    add_tail(by_var, 3, true);
    add_tail(by_corr, 1, false);
    add_tail(by_corr, 2, false);
    add_tail(by_corr, 3, false);
    add_tail(by_var, 1, false);
    add_tail(by_var, 2, false);
    add_tail(by_var, 3, false);

    return cands;
}

int best_anchor_col(const std::vector<std::vector<double>>& rows, const std::vector<int>& g0, const std::vector<int>& g1) {
    if (g0.empty()) return -1;
    if (g1.empty()) return g0[0];
    auto sample = sample_rows(rows);
    int best_col = g0[0];
    double best_score = -1.0;
    for (int c0 : g0) {
        double corr_s = 0.0;
        for (int c1 : g1) corr_s += abs_corr(sample, c0, c1);
        corr_s /= static_cast<double>(g1.size());
        double uniq_s = unique_ratio(sample, c0);
        double score = 0.7 * corr_s + 0.3 * uniq_s;
        if (score > best_score) {
            best_score = score;
            best_col = c0;
        }
    }
    return best_col;
}

bool candidate_has_finite_range(const std::vector<std::vector<double>>& rows, const std::vector<int>& cols) {
    for (int c : cols) {
        bool has = false;
        double mn = 0.0, mx = 0.0;
        for (const auto& r : rows) {
            double v = r[c];
            if (!std::isfinite(v)) continue;
            if (!has) {
                mn = mx = v;
                has = true;
            } else {
                mn = std::min(mn, v);
                mx = std::max(mx, v);
            }
        }
        if (!has || mn > mx) return false;
    }
    return true;
}

std::string make_group_csv(
    const std::vector<std::vector<double>>& rows,
    const std::vector<int>& cols,
    const std::string& tag,
    int extra_anchor_col = -1
) {
    std::filesystem::path tmp = std::filesystem::temp_directory_path() / ("aclustersub_g_" + tag + ".csv");
    std::ofstream out(tmp.string());
    if (!out.is_open()) throw std::runtime_error("AClusterSub cannot create temp file");
    out << std::fixed << std::setprecision(6);
    for (const auto& row : rows) {
        for (size_t i = 0; i < cols.size(); ++i) {
            if (i) out << ",";
            out << row[cols[i]];
        }
        if (extra_anchor_col >= 0) out << "," << row[extra_anchor_col];
        out << "\n";
    }
    return tmp.string();
}

std::vector<std::vector<int>> build_fixed_groups(int dim, int group_size) {
    std::vector<std::vector<int>> groups;
    if (dim <= 0) return groups;
    int gs = std::max(1, group_size);
    for (int s = 0; s < dim; s += gs) {
        std::vector<int> g;
        for (int c = s; c < std::min(dim, s + gs); ++c) g.push_back(c);
        groups.push_back(std::move(g));
    }
    if (groups.size() == 1 && dim > 1) {
        groups[0].pop_back();
        groups.push_back({dim - 1});
    }
    return groups;
}

std::vector<std::vector<int>> build_offset_groups(int dim, int group_size, int offset) {
    std::vector<std::vector<int>> groups;
    if (dim <= 0) return groups;
    int gs = std::max(1, group_size);
    std::vector<int> order(dim);
    for (int i = 0; i < dim; ++i) order[i] = (i + offset) % dim;
    for (int s = 0; s < dim; s += gs) {
        std::vector<int> g;
        for (int i = s; i < std::min(dim, s + gs); ++i) g.push_back(order[i]);
        std::sort(g.begin(), g.end());
        groups.push_back(std::move(g));
    }
    if (groups.size() == 1 && dim > 1) {
        groups[0].pop_back();
        groups.push_back({dim - 1});
    }
    return groups;
}

std::vector<std::vector<int>> build_variance_groups(const std::vector<std::vector<double>>& rows, int group_size) {
    int dim = static_cast<int>(rows[0].size());
    std::vector<std::pair<double, int>> ord;
    ord.reserve(dim);
    for (int c = 0; c < dim; ++c) ord.push_back({column_variance(rows, c), c});
    std::sort(ord.begin(), ord.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    int gs = std::max(1, group_size);
    std::vector<std::vector<int>> groups;
    for (int s = 0; s < dim; s += gs) {
        std::vector<int> g;
        for (int i = s; i < std::min(dim, s + gs); ++i) g.push_back(ord[i].second);
        std::sort(g.begin(), g.end());
        groups.push_back(std::move(g));
    }
    if (groups.size() == 1 && dim > 1) {
        groups[0].pop_back();
        groups.push_back({dim - 1});
    }
    return groups;
}

std::vector<std::vector<int>> build_corr_groups_chunked(const std::vector<std::vector<double>>& rows, int group_size) {
    int dim = static_cast<int>(rows[0].size());
    auto sample = sample_rows(rows);
    std::vector<std::pair<double, int>> ord;
    ord.reserve(dim);
    for (int c = 0; c < dim; ++c) {
        double ac = 0.0;
        for (int j = 0; j < dim; ++j) if (j != c) ac += abs_corr(sample, c, j);
        ac /= std::max(1, dim - 1);
        ord.push_back({ac, c});
    }
    std::sort(ord.begin(), ord.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    int gs = std::max(1, group_size);
    std::vector<std::vector<int>> groups;
    for (int s = 0; s < dim; s += gs) {
        std::vector<int> g;
        for (int i = s; i < std::min(dim, s + gs); ++i) g.push_back(ord[i].second);
        std::sort(g.begin(), g.end());
        groups.push_back(std::move(g));
    }
    if (groups.size() == 1 && dim > 1) {
        groups[0].pop_back();
        groups.push_back({dim - 1});
    }
    return groups;
}

std::string groups_signature(const std::vector<std::vector<int>>& groups) {
    std::string s;
    for (const auto& g : groups) {
        s.push_back('[');
        for (int c : g) {
            s += std::to_string(c);
            s.push_back('_');
        }
        s.push_back(']');
    }
    return s;
}

std::vector<std::vector<std::vector<int>>> build_grouping_candidates(const std::vector<std::vector<double>>& rows, int group_size) {
    int dim = static_cast<int>(rows[0].size());
    std::vector<std::vector<std::vector<int>>> cands;
    std::unordered_set<std::string> seen;
    auto add = [&](std::vector<std::vector<int>> g) {
        if (g.empty()) return;
        std::string sig = groups_signature(g);
        if (seen.count(sig)) return;
        seen.insert(sig);
        cands.push_back(std::move(g));
    };
    std::vector<int> group_sizes;
    group_sizes.push_back(8);
    if (dim <= 32) group_sizes.push_back(1);
    if (dim <= 48) group_sizes.push_back(2);
    if (dim <= 64) group_sizes.push_back(4);
    group_sizes.push_back(std::max(1, group_size));
    for (int gs : group_sizes) {
        if (gs <= 0 || gs >= dim) continue;
        add(build_fixed_groups(dim, gs));
    }
    return cands;
}

int choose_group_k(int base_k, int full_dim, int group_dim) {
    if (group_dim <= 0) return base_k;
    double scale = static_cast<double>(full_dim) / static_cast<double>(group_dim);
    int kg = static_cast<int>(std::llround(base_k * scale));
    kg = std::max(32, std::min(kg, 1024));
    return kg;
}

struct GlobalKMeansResult {
    std::vector<std::vector<double>> centroids;
    std::vector<uint32_t> assignment;
};

GlobalKMeansResult run_global_kmeans(const std::vector<std::vector<double>>& rows, int k, int max_iter) {
    int n = static_cast<int>(rows.size());
    int dim = static_cast<int>(rows[0].size());
    k = std::max(1, std::min(k, n));
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::mt19937 rng(42);
    std::shuffle(idx.begin(), idx.end(), rng);
    std::vector<std::vector<double>> centroids;
    centroids.reserve(k);
    for (int i = 0; i < k; ++i) centroids.push_back(rows[idx[i]]);
    std::vector<uint32_t> assign(n, 0);
    for (int it = 0; it < max_iter; ++it) {
        for (int i = 0; i < n; ++i) {
            double best = std::numeric_limits<double>::infinity();
            int best_j = 0;
            for (int j = 0; j < k; ++j) {
                double d2 = 0.0;
                for (int c = 0; c < dim; ++c) {
                    double d = rows[i][c] - centroids[j][c];
                    d2 += d * d;
                }
                if (d2 < best) {
                    best = d2;
                    best_j = j;
                }
            }
            assign[i] = static_cast<uint32_t>(best_j);
        }
        std::vector<std::vector<double>> sums(k, std::vector<double>(dim, 0.0));
        std::vector<int> cnt(k, 0);
        for (int i = 0; i < n; ++i) {
            int j = static_cast<int>(assign[i]);
            cnt[j] += 1;
            for (int c = 0; c < dim; ++c) sums[j][c] += rows[i][c];
        }
        for (int j = 0; j < k; ++j) {
            if (cnt[j] == 0) continue;
            for (int c = 0; c < dim; ++c) centroids[j][c] = sums[j][c] / static_cast<double>(cnt[j]);
        }
    }
    return {centroids, assign};
}

std::string key_from_values(const std::vector<double>& row, const std::vector<int>& cols) {
    std::string k;
    k.reserve(cols.size() * 24);
    for (int c : cols) {
        long long q = static_cast<long long>(std::llround(row[c] * kQuant));
        k += std::to_string(q);
        k.push_back('|');
    }
    return k;
}

std::string key_from_decoded(const std::vector<double>& row, size_t take_cols) {
    std::string k;
    k.reserve(take_cols * 24);
    for (size_t i = 0; i < take_cols; ++i) {
        double v = row[i];
        long long q = static_cast<long long>(std::llround(v * kQuant));
        k += std::to_string(q);
        k.push_back('|');
    }
    return k;
}

std::vector<uint32_t> map_rows(
    const std::vector<std::vector<double>>& original,
    const std::vector<std::vector<double>>& decoded,
    const std::vector<int>& cols,
    size_t decoded_take_cols
) {
    std::unordered_map<std::string, std::vector<uint32_t>> buckets;
    buckets.reserve(original.size() * 2);
    for (uint32_t i = 0; i < static_cast<uint32_t>(original.size()); ++i) {
        buckets[key_from_values(original[i], cols)].push_back(i);
    }
    for (auto& kv : buckets) std::reverse(kv.second.begin(), kv.second.end());
    std::vector<uint32_t> ids;
    ids.reserve(decoded.size());
    for (const auto& r : decoded) {
        auto k = key_from_decoded(r, decoded_take_cols);
        auto it = buckets.find(k);
        if (it == buckets.end() || it->second.empty()) throw std::runtime_error("AClusterSub mapping failed");
        ids.push_back(it->second.back());
        it->second.pop_back();
    }
    return ids;
}

}

std::vector<uint8_t> AClusterSub::encode_multidim(
    const std::string& csv_file_path,
    int& dim,
    int k,
    int page_size,
    int pack_size,
    int block_size
) {
    auto rows = read_csv_rows(csv_file_path);
    if (rows.empty()) return {};
    int full_dim = static_cast<int>(rows[0].size());
    dim = full_dim;

    auto groups = build_fixed_groups(full_dim, dimension_group_size_);
    if (groups.empty()) throw std::runtime_error("AClusterSub empty groups");
    int k_global = std::max(8, std::min(k, 256));
    auto km = run_global_kmeans(rows, k_global, 10);

    std::vector<uint8_t> centroid_blob;
    centroid_blob.reserve(static_cast<size_t>(k_global) * static_cast<size_t>(full_dim) * 3);
    for (int j = 0; j < k_global; ++j) {
        for (int c = 0; c < full_dim; ++c) {
            int64_t q = static_cast<int64_t>(std::llround(km.centroids[j][c] * kQuant));
            append_svarint(centroid_blob, q);
        }
    }
    std::vector<uint8_t> assign_blob;
    assign_blob.reserve(rows.size());
    for (uint32_t id : km.assignment) append_uvarint(assign_blob, id);

    std::vector<std::vector<uint8_t>> residual_streams(groups.size());
    uint64_t residual_bits = 0;
    for (size_t gi = 0; gi < groups.size(); ++gi) {
        auto& rs = residual_streams[gi];
        rs.reserve(rows.size() * groups[gi].size());
        for (size_t r = 0; r < rows.size(); ++r) {
            int cid = static_cast<int>(km.assignment[r]);
            for (int c : groups[gi]) {
                int64_t q = static_cast<int64_t>(std::llround((rows[r][c] - km.centroids[cid][c]) * kQuant));
                append_svarint(rs, q);
            }
        }
        residual_bits += static_cast<uint64_t>(rs.size()) * 8ULL;
    }

    uint64_t base_bits = static_cast<uint64_t>(centroid_blob.size()) * 8ULL;
    uint64_t id_bits = static_cast<uint64_t>(assign_blob.size()) * 8ULL;
    uint64_t header_bytes = 6 + 4 * 6 + 8 * 4 + 4 * 2 + static_cast<uint64_t>(groups.size()) * 8ULL + static_cast<uint64_t>(full_dim) * 4ULL;
    uint64_t header_bits = header_bytes * 8ULL;

    std::vector<uint8_t> out;
    out.push_back('A'); out.push_back('S'); out.push_back('U'); out.push_back('B');
    out.push_back(1);
    out.push_back(1);
    append_u32(out, static_cast<uint32_t>(full_dim));
    append_u32(out, static_cast<uint32_t>(rows.size()));
    append_u32(out, static_cast<uint32_t>(groups.size()));
    append_u32(out, static_cast<uint32_t>(k_global));
    append_u32(out, static_cast<uint32_t>(centroid_blob.size()));
    append_u32(out, static_cast<uint32_t>(assign_blob.size()));
    append_u64(out, header_bits);
    append_u64(out, base_bits);
    append_u64(out, residual_bits);
    append_u64(out, id_bits);
    for (const auto& g : groups) append_u32(out, static_cast<uint32_t>(g.size()));
    for (const auto& rs : residual_streams) append_u32(out, static_cast<uint32_t>(rs.size()));
    for (const auto& g : groups) for (int c : g) append_u32(out, static_cast<uint32_t>(c));
    out.insert(out.end(), centroid_blob.begin(), centroid_blob.end());
    out.insert(out.end(), assign_blob.begin(), assign_blob.end());
    for (const auto& rs : residual_streams) out.insert(out.end(), rs.begin(), rs.end());
    return out;
}

std::vector<std::vector<double>> AClusterSub::decode_multidim(const std::vector<uint8_t>& compressed_data) {
    if (compressed_data.size() < 7) throw std::runtime_error("AClusterSub decode stream too small");
    if (!(compressed_data[0] == 'A' && compressed_data[1] == 'S' && compressed_data[2] == 'U' && compressed_data[3] == 'B')) {
        throw std::runtime_error("AClusterSub decode bad magic");
    }
    uint8_t mode = compressed_data[4];
    size_t pos = 5;
    if (mode != 1) throw std::runtime_error("AClusterSub decode unsupported mode");
    uint8_t strategy = compressed_data[pos++];
    if (strategy == 1) {
        int full_dim = static_cast<int>(read_u32(compressed_data, pos));
        int row_count = static_cast<int>(read_u32(compressed_data, pos));
        int group_count = static_cast<int>(read_u32(compressed_data, pos));
        int k_global = static_cast<int>(read_u32(compressed_data, pos));
        int centroid_sz = static_cast<int>(read_u32(compressed_data, pos));
        int assign_sz = static_cast<int>(read_u32(compressed_data, pos));
        read_u64(compressed_data, pos);
        read_u64(compressed_data, pos);
        read_u64(compressed_data, pos);
        read_u64(compressed_data, pos);
        std::vector<int> group_dims(group_count), residual_sz(group_count);
        for (int i = 0; i < group_count; ++i) group_dims[i] = static_cast<int>(read_u32(compressed_data, pos));
        for (int i = 0; i < group_count; ++i) residual_sz[i] = static_cast<int>(read_u32(compressed_data, pos));
        std::vector<std::vector<int>> groups(group_count);
        for (int i = 0; i < group_count; ++i) {
            groups[i].resize(group_dims[i]);
            for (int j = 0; j < group_dims[i]; ++j) groups[i][j] = static_cast<int>(read_u32(compressed_data, pos));
        }
        if (pos + centroid_sz + assign_sz > compressed_data.size()) throw std::runtime_error("AClusterSub decode strategy1 truncated");
        std::vector<uint8_t> centroid_blob(compressed_data.begin() + pos, compressed_data.begin() + pos + centroid_sz); pos += centroid_sz;
        std::vector<uint8_t> assign_blob(compressed_data.begin() + pos, compressed_data.begin() + pos + assign_sz); pos += assign_sz;
        std::vector<std::vector<uint8_t>> residual_streams(group_count);
        for (int i = 0; i < group_count; ++i) {
            if (pos + residual_sz[i] > compressed_data.size()) throw std::runtime_error("AClusterSub decode residual truncated");
            residual_streams[i] = std::vector<uint8_t>(compressed_data.begin() + pos, compressed_data.begin() + pos + residual_sz[i]);
            pos += residual_sz[i];
        }

        size_t cp = 0;
        std::vector<std::vector<double>> centroids(k_global, std::vector<double>(full_dim, 0.0));
        for (int j = 0; j < k_global; ++j) for (int c = 0; c < full_dim; ++c) centroids[j][c] = static_cast<double>(read_svarint(centroid_blob, cp)) / kQuant;
        size_t ap = 0;
        std::vector<uint32_t> assign(row_count, 0);
        for (int i = 0; i < row_count; ++i) assign[i] = static_cast<uint32_t>(read_uvarint(assign_blob, ap));

        std::vector<std::vector<double>> out(row_count, std::vector<double>(full_dim, 0.0));
        for (int gi = 0; gi < group_count; ++gi) {
            size_t rp = 0;
            for (int r = 0; r < row_count; ++r) {
                int cid = static_cast<int>(assign[r]);
                if (cid < 0 || cid >= k_global) throw std::runtime_error("AClusterSub decode bad centroid id");
                for (int j = 0; j < group_dims[gi]; ++j) {
                    double res = static_cast<double>(read_svarint(residual_streams[gi], rp)) / kQuant;
                    int c = groups[gi][j];
                    out[r][c] = centroids[cid][c] + res;
                }
            }
        }
        return out;
    }
    if (strategy != 0) throw std::runtime_error("AClusterSub decode unsupported strategy");

    int full_dim = static_cast<int>(read_u32(compressed_data, pos));
    int row_count = static_cast<int>(read_u32(compressed_data, pos));
    int group_count = static_cast<int>(read_u32(compressed_data, pos));
    read_u64(compressed_data, pos);
    read_u64(compressed_data, pos);
    read_u64(compressed_data, pos);
    if (group_count <= 0) throw std::runtime_error("AClusterSub decode bad group count");
    std::vector<int> group_dims(group_count);
    std::vector<int> blob_sizes(group_count);
    std::vector<int> pos_sizes(std::max(0, group_count - 1), 0);
    for (int i = 0; i < group_count; ++i) group_dims[i] = static_cast<int>(read_u32(compressed_data, pos));
    for (int i = 0; i < group_count; ++i) blob_sizes[i] = static_cast<int>(read_u32(compressed_data, pos));
    for (int i = 0; i < group_count - 1; ++i) pos_sizes[i] = static_cast<int>(read_u32(compressed_data, pos));
    std::vector<std::vector<int>> groups(group_count);
    for (int i = 0; i < group_count; ++i) {
        groups[i].resize(group_dims[i]);
        for (int j = 0; j < group_dims[i]; ++j) groups[i][j] = static_cast<int>(read_u32(compressed_data, pos));
    }

    std::vector<std::vector<uint8_t>> blobs(group_count);
    for (int i = 0; i < group_count; ++i) {
        if (pos + blob_sizes[i] > compressed_data.size()) throw std::runtime_error("AClusterSub decode blob truncated");
        blobs[i] = std::vector<uint8_t>(compressed_data.begin() + pos, compressed_data.begin() + pos + blob_sizes[i]);
        pos += blob_sizes[i];
    }
    std::vector<std::vector<uint8_t>> pos_streams(group_count > 0 ? group_count - 1 : 0);
    for (int i = 0; i < group_count - 1; ++i) {
        if (pos + pos_sizes[i] > compressed_data.size()) throw std::runtime_error("AClusterSub decode pos truncated");
        pos_streams[i] = std::vector<uint8_t>(compressed_data.begin() + pos, compressed_data.begin() + pos + pos_sizes[i]);
        pos += pos_sizes[i];
    }

    ACluster base;
    std::vector<std::vector<std::vector<double>>> dec_groups(group_count);
    for (int i = 0; i < group_count; ++i) dec_groups[i] = base.decode_multidim(blobs[i]);
    if (dec_groups[0].size() != static_cast<size_t>(row_count)) throw std::runtime_error("AClusterSub decode row mismatch");

    std::vector<std::vector<double>> out(row_count, std::vector<double>(full_dim, 0.0));
    for (size_t r = 0; r < dec_groups[0].size(); ++r) {
        for (int j = 0; j < group_dims[0]; ++j) out[r][groups[0][j]] = dec_groups[0][r][j];
    }
    for (int gi = 1; gi < group_count; ++gi) {
        auto pos_map = decode_positions(pos_streams[gi - 1], static_cast<uint32_t>(row_count));
        if (dec_groups[gi].size() != pos_map.size()) throw std::runtime_error("AClusterSub decode group size mismatch");
        std::vector<uint8_t> seen(row_count, 0);
        for (size_t r = 0; r < dec_groups[gi].size(); ++r) {
            uint32_t p = pos_map[r];
            if (p >= static_cast<uint32_t>(row_count) || seen[p]) throw std::runtime_error("AClusterSub decode mapped row invalid");
            seen[p] = 1;
            for (int j = 0; j < group_dims[gi]; ++j) out[p][groups[gi][j]] = dec_groups[gi][r][j];
        }
    }
    return out;
}
