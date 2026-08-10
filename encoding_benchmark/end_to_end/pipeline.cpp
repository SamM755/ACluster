#include "pipeline.h"

#include <charconv>
#include <cstring>

namespace {
using Clock = std::chrono::steady_clock;
long long elapsed(Clock::time_point start) {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(Clock::now() - start).count();
}
uint64_t bits(double v) { uint64_t u; std::memcpy(&u, &v, 8); return u; }
double value(uint64_t u) { double v; std::memcpy(&v, &u, 8); return v; }
void put(std::vector<uint8_t>& out, uint64_t v, int width) {
    for (int b = 0; b < width; ++b) out.push_back(static_cast<uint8_t>(v >> (8 * b)));
}
void payload(std::vector<uint8_t>& out, const std::vector<uint8_t>& data) {
    put(out, data.size(), 8);
    out.insert(out.end(), data.begin(), data.end());
}
struct Reader {
    const std::vector<uint8_t>& data;
    size_t pos = 0;
    uint64_t get(int width) {
        if (width < 1 || width > 8 || data.size() - pos < static_cast<size_t>(width))
            throw std::runtime_error("Truncated end-to-end stream.");
        uint64_t v = 0;
        for (int b = 0; b < width; ++b) v |= uint64_t(data[pos++]) << (8 * b);
        return v;
    }
    std::vector<uint8_t> block() {
        const uint64_t n = get(8);
        if (n > data.size() - pos) throw std::runtime_error("Invalid end-to-end stream length.");
        std::vector<uint8_t> out(data.begin() + pos, data.begin() + pos + n);
        pos += n;
        return out;
    }
    void finish() const {
        if (pos != data.size()) throw std::runtime_error("Unconsumed end-to-end bytes.");
    }
};
constexpr uint32_t table_magic = 0x31425445;
constexpr uint32_t page_magic = 0x31505445;
constexpr uint32_t baseline_magic = 0x32505445;
const std::vector<std::string> baseline_methods{"Chimp", "Elf", "Gorilla", "RLE", "Huffman", "REGER", "Vortex"};

std::vector<uint8_t> scalar_encode(const std::string& method, const std::vector<double>& column) {
    if (method == "Chimp") return Chimp().encode(column);
    if (method == "Elf") return Elf().encode(column);
    if (method == "Gorilla") return Gorilla().encode(column);
    if (method == "RLE") return Rle().encode(column);
    if (method == "Huffman") return Huffman().encode(column);
    throw std::invalid_argument("Unknown scalar baseline.");
}
std::vector<double> scalar_decode(const std::string& method, const std::vector<uint8_t>& bytes) {
    if (method == "Chimp") return Chimp().decode(bytes);
    if (method == "Elf") return Elf().decode(bytes);
    if (method == "Gorilla") return Gorilla().decode(bytes);
    if (method == "RLE") return Rle().decode(bytes);
    if (method == "Huffman") return Huffman().decode(bytes);
    throw std::invalid_argument("Unknown scalar baseline.");
}
using Corrections = std::vector<std::vector<std::pair<int, uint64_t>>>;

void scaling_corrections(const std::vector<double>& column, int c, Corrections& patches) {
    int decimals = 0;
    for (double v : column) {
        const auto s = std::to_string(v);
        const auto dot = s.find('.'), last = s.find_last_not_of('0');
        if (dot != std::string::npos && last != std::string::npos && last > dot)
            decimals = std::max(decimals, static_cast<int>(last-dot));
    }
    const double factor = std::pow(10.0, decimals);
    for (size_t r = 0; r < column.size(); ++r) {
        const double scaled = std::round(column[r]*factor);
        if (!std::isfinite(scaled) || scaled < -0x1p63 || scaled >= 0x1p63)
            throw std::runtime_error("Baseline int64 scaling overflow.");
        const double restored = static_cast<int64_t>(scaled)/factor;
        if (bits(restored) != bits(column[r])) patches[r].push_back({c, bits(column[r])});
    }
}
void write_corrections(std::vector<uint8_t>& out, const Corrections& patches,
                       const std::vector<size_t>& order, int dim, EndToEndMetrics& metrics) {
    uint64_t count = 0;
    for (const auto& row : patches) count += row.size();
    put(out, count, 8);
    for (size_t r = 0; r < order.size(); ++r) for (const auto& patch : patches[order[r]]) {
        put(out, r*dim+patch.first, 4); put(out, patch.second, 8);
    }
    metrics.correction_values += count;
    metrics.correction_bytes += count*12;
}
void read_corrections(Reader& reader, EndToEndRows& rows, int dim) {
    const uint64_t count = reader.get(8);
    if (count > rows.size()*dim) throw std::runtime_error("Invalid correction count.");
    for (uint64_t p = 0; p < count; ++p) {
        const uint64_t index = reader.get(4);
        if (index >= rows.size()*dim) throw std::runtime_error("Invalid correction position.");
        rows[index/dim][index%dim] = value(reader.get(8));
    }
}

std::vector<uint8_t> encode_baseline_page(const EndToEndRows& rows, size_t begin, size_t end,
                                         const std::string& method, EndToEndMetrics& metrics) {
    const auto method_it = std::find(baseline_methods.begin(), baseline_methods.end(), method);
    const int dim = static_cast<int>(rows[0].size());
    const size_t n = end-begin;
    std::vector<uint8_t> out;
    put(out, baseline_magic, 4); put(out, n, 4); put(out, dim, 4);
    put(out, method_it-baseline_methods.begin(), 4);
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    Corrections patches(n);
    auto start = Clock::now();
    if (method == "REGER") {
        std::vector<std::string> lines(n);
        char buffer[128];
        for (size_t r = 0; r < n; ++r) for (int c = 0; c < dim; ++c) {
            const auto converted = std::to_chars(buffer, buffer+sizeof(buffer), rows[begin+r][c]);
            if (converted.ec != std::errc{}) throw std::runtime_error("Binary64 conversion failed.");
            if (c) lines[r].push_back(',');
            lines[r].append(buffer, converted.ptr);
        }
        int parsed_dim = dim;
        auto preprocessed = ClusterLogic::preprocessPageFromStrings(lines, parsed_dim);
        std::vector<double> factors(dim);
        std::vector<std::vector<int64_t>> scaled(n, std::vector<int64_t>(dim));
        for (int c = 0; c < dim; ++c) {
            factors[c] = std::pow(10.0, preprocessed.max_decimal_places[c]);
            for (size_t r = 0; r < n; ++r) {
                const __int128 integer = static_cast<__int128>(preprocessed.data->columns[c][r])+preprocessed.min_values[c];
                if (integer < INT64_MIN || integer > INT64_MAX) throw std::runtime_error("REGER scaling overflow.");
                scaled[r][c] = static_cast<int64_t>(integer);
                if (bits(scaled[r][c]/factors[c]) != bits(rows[begin+r][c]))
                    patches[r].push_back({c, bits(rows[begin+r][c])});
            }
        }
        metrics.preprocessing_ns += elapsed(start);
        start = Clock::now();
        auto bytes = Reger().encode_multidim(scaled, factors, &order);
        metrics.baseline_bytes += bytes.size();
        payload(out, bytes);
        metrics.baseline_ns += elapsed(start);
    } else {

        payload(out, {});
        if (method == "Vortex") {
            EndToEndRows page(rows.begin()+begin, rows.begin()+end);
            auto permutation = Vortex().compute_permutation(page);
            order.assign(permutation.begin(), permutation.end());
            std::vector<double> permutation_values(order.begin(), order.end());
            auto bytes = Rle().encode(permutation_values);
            metrics.baseline_bytes += bytes.size();
            payload(out, bytes);
        }
        metrics.baseline_ns += elapsed(start);
        for (int c = 0; c < dim; ++c) {
            start = Clock::now();
            std::vector<double> column(n);
            for (size_t r = 0; r < n; ++r) column[r] = rows[begin+order[r]][c];

            if (method == "RLE" || method == "Huffman" || method == "Vortex")
                scaling_corrections(column, c, patches);
            metrics.preprocessing_ns += elapsed(start);
            start = Clock::now();
            auto bytes = scalar_encode(method == "Vortex" ? "RLE" : method, column);
            metrics.baseline_bytes += bytes.size();
            payload(out, bytes);
            metrics.baseline_ns += elapsed(start);
        }
        std::iota(order.begin(), order.end(), 0);
    }
    write_corrections(out, patches, order, dim, metrics);
    ++metrics.pages;
    return out;
}

EndToEndRows decode_baseline_page(Reader& reader, size_t n, int dim, uint64_t id) {
    if (id >= baseline_methods.size()) throw std::runtime_error("Unknown baseline page method.");
    const auto& method = baseline_methods[id];
    EndToEndRows rows(n, std::vector<double>(dim));
    std::vector<double> permutation;
    if (method == "REGER") {
        const auto decoded = Reger().decode_multidim(reader.block());
        if (decoded.first.size() != n || decoded.second.size() != static_cast<size_t>(dim))
            throw std::runtime_error("REGER page shape mismatch.");
        for (size_t r = 0; r < n; ++r) {
            if (decoded.first[r].size() != static_cast<size_t>(dim)) throw std::runtime_error("REGER row shape mismatch.");
            for (int c = 0; c < dim; ++c) rows[r][c] = decoded.first[r][c]/decoded.second[c];
        }
    } else {
        if (!reader.block().empty()) throw std::runtime_error("Unexpected scalar multi-column payload.");
        if (method == "Vortex") permutation = Rle().decode(reader.block());
        for (int c = 0; c < dim; ++c) {
            auto column = scalar_decode(method == "Vortex" ? "RLE" : method, reader.block());
            if (column.size() != n) throw std::runtime_error("Baseline row count mismatch.");
            for (size_t r = 0; r < n; ++r) rows[r][c] = column[r];
        }
    }
    read_corrections(reader, rows, dim);
    reader.finish();
    if (method == "Vortex") {
        if (permutation.size() != n) throw std::runtime_error("Vortex permutation length mismatch.");
        EndToEndRows restored(n);
        for (size_t r = 0; r < n; ++r) {
            const double p = permutation[r];
            if (!std::isfinite(p) || p < 0 || p >= n || std::floor(p) != p || !restored[static_cast<size_t>(p)].empty())
                throw std::runtime_error("Invalid Vortex permutation.");
            restored[static_cast<size_t>(p)] = std::move(rows[r]);
        }
        return restored;
    }
    return rows;
}

std::vector<int> complement(int dim, const std::vector<int>& selected) {
    std::vector<bool> used(dim, false);
    for (int c : selected) {
        if (c < 0 || c >= dim || used[c]) throw std::runtime_error("Invalid selected column ID.");
        used[c] = true;
    }
    std::vector<int> remaining;
    for (int c = 0; c < dim; ++c) if (!used[c]) remaining.push_back(c);
    return remaining;
}

std::vector<uint8_t> encode_page(const EndToEndRows& rows, size_t begin, size_t end,
    const std::vector<int>& selected, const std::string& method, uint32_t seed,
    int page_size, EndToEndMetrics& metrics)
{
    const int dim = static_cast<int>(rows[0].size());
    const size_t n = end - begin;
    auto start = Clock::now();
    const auto remaining = complement(dim, selected);

    std::vector<std::vector<std::pair<int, uint64_t>>> corrections(n);
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    PreprocessorResult preprocessed;
    if (!selected.empty()) {
        std::vector<std::string> lines(n);
        char buffer[128];
        for (size_t r = 0; r < n; ++r) {
            for (int c : selected) {
                const auto converted = std::to_chars(buffer, buffer + sizeof(buffer), rows[begin+r][c]);
                if (converted.ec != std::errc{}) throw std::runtime_error("Binary64 conversion failed.");
                if (!lines[r].empty()) lines[r].push_back(',');
                lines[r].append(buffer, converted.ptr);
            }
        }
        int selected_dim = static_cast<int>(selected.size());
        preprocessed = ClusterLogic::preprocessPageFromStrings(lines, selected_dim);
        if (!preprocessed.data || preprocessed.data->num_rows != static_cast<int>(n))
            throw std::runtime_error("Preprocessing changed the page shape.");
        for (size_t j = 0; j < selected.size(); ++j) {
            const double scale = std::pow(10.0, preprocessed.max_decimal_places[j]);
            for (size_t r = 0; r < n; ++r) {
                const double restored = static_cast<double>(preprocessed.data->columns[j][r] +
                    preprocessed.min_values[j]) / scale;
                if (bits(restored) != bits(rows[begin+r][selected[j]]))
                    corrections[r].push_back({selected[j], bits(rows[begin+r][selected[j]])});
            }
        }
    }

    for (size_t r = 0; r < n; ++r)
        for (int c : remaining)
            if (rows[begin+r][c] == 0 && std::signbit(rows[begin+r][c]))
                corrections[r].push_back({c, bits(rows[begin+r][c])});
    metrics.preprocessing_ns += elapsed(start);

    start = Clock::now();
    std::vector<uint8_t> cluster_bytes;
    if (!selected.empty()) {
        auto clustered = ClusterLogic::encodePreprocessedPage(std::move(preprocessed), method,
            100, page_size, 10, 10, seed);
        order = std::move(clustered.row_order);
        cluster_bytes = std::move(clustered.bytes);
        metrics.references += clustered.reference_count;
    }
    metrics.cluster_ns += elapsed(start);

    std::vector<uint8_t> out;
    put(out, page_magic, 4); put(out, n, 4); put(out, dim, 4); put(out, selected.size(), 4);
    for (int c : selected) put(out, c, 4);
    payload(out, cluster_bytes);
    metrics.cluster_bytes += cluster_bytes.size();
    start = Clock::now();
    for (int c : remaining) {
        std::vector<double> column(n);
        for (size_t r = 0; r < n; ++r) column[r] = rows[begin+order[r]][c];
        Alp encoder;
        auto bytes = encoder.encode(column);
        metrics.alp_bytes += bytes.size();
        payload(out, bytes);
    }
    metrics.alp_ns += elapsed(start);
    uint64_t count = 0;
    for (const auto& row : corrections) count += row.size();
    put(out, count, 8);
    for (size_t r = 0; r < n; ++r) {
        for (const auto& patch : corrections[order[r]]) {
            put(out, r * dim + patch.first, 4);
            put(out, patch.second, 8);
        }
    }
    metrics.correction_values += count;
    metrics.correction_bytes += count * 12;
    ++metrics.pages;
    return out;
}

EndToEndRows decode_page(const std::vector<uint8_t>& bytes, int expected_dim) {
    Reader reader{bytes};
    const auto magic = reader.get(4);
    if (magic != page_magic && magic != baseline_magic) throw std::runtime_error("Invalid page magic/version.");
    const size_t n = reader.get(4);
    const int dim = static_cast<int>(reader.get(4));
    const size_t selected_dim = reader.get(4);
    if (magic == baseline_magic) {
        if (dim != expected_dim || !n || n > 10000) throw std::runtime_error("Invalid baseline page shape.");
        return decode_baseline_page(reader, n, dim, selected_dim);
    }
    if (dim != expected_dim || !n || n > 10000 || selected_dim > static_cast<size_t>(dim) || selected_dim > 255)
        throw std::runtime_error("Invalid page shape.");
    std::vector<int> selected(selected_dim);
    for (auto& c : selected) c = static_cast<int>(reader.get(4));
    const auto remaining = complement(dim, selected);
    EndToEndRows rows(n, std::vector<double>(dim));
    auto cluster = reader.block();
    if (!selected.empty()) {
        const auto decoded = ClusterLogic::decode_multidim_impl(cluster);
        if (decoded.size() != n) throw std::runtime_error("Cluster row count mismatch.");
        for (size_t r = 0; r < n; ++r) {
            if (decoded[r].size() != selected_dim) throw std::runtime_error("Cluster dimension mismatch.");
            for (size_t j = 0; j < selected_dim; ++j) rows[r][selected[j]] = decoded[r][j];
        }
    } else if (!cluster.empty()) throw std::runtime_error("Unexpected cluster payload.");
    for (int c : remaining) {
        Alp decoder;
        auto column = decoder.decode(reader.block());
        if (column.size() != n) throw std::runtime_error("ALP row count mismatch.");
        for (size_t r = 0; r < n; ++r) rows[r][c] = column[r];
    }
    const uint64_t count = reader.get(8);
    if (count > n * dim) throw std::runtime_error("Invalid correction count.");
    for (uint64_t p = 0; p < count; ++p) {
        const uint64_t index = reader.get(4);
        if (index >= n * dim) throw std::runtime_error("Invalid correction position.");
        rows[index / dim][index % dim] = value(reader.get(8));
    }
    reader.finish();
    return rows;
}
}

EndToEndRows read_end_to_end_table(const std::string& path, size_t expected_rows, int dim) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("Cannot open full table: " + path);
    EndToEndRows rows;
    rows.reserve(expected_rows);
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) throw std::runtime_error("Empty row in full table.");
        std::vector<double> row;
        size_t begin = 0;
        while (true) {
            size_t end = line.find(',', begin);
            if (end == std::string::npos) end = line.size();
            double v;
            auto parsed = std::from_chars(line.data()+begin, line.data()+end, v);
            if (parsed.ec != std::errc{} || parsed.ptr != line.data()+end || !std::isfinite(v) || std::abs(v) >= 1e17)
                throw std::runtime_error("Invalid/unsupported finite numeric cell at row " + std::to_string(rows.size()+1));
            row.push_back(v);
            if (end == line.size()) break;
            begin = end+1;
        }
        if (row.size() != static_cast<size_t>(dim)) throw std::runtime_error("Full-table dimension mismatch.");
        rows.push_back(std::move(row));
    }
    if (!file.eof() || rows.size() != expected_rows || rows.empty())
        throw std::runtime_error("Full-table row count mismatch: " + path);
    return rows;
}

EndToEndSelection select_end_to_end_columns(const EndToEndRows& rows) {
    if (rows.empty() || rows[0].empty()) throw std::invalid_argument("Empty selector input.");
    const int dim = static_cast<int>(rows[0].size());
    const size_t n = std::min<size_t>(10000, rows.size());
    EndToEndSelection result;
    result.columns.resize(dim);
    std::iota(result.columns.begin(), result.columns.end(), 0);
    result.dimension_path.push_back(dim);
    if (dim == 1) return result;

    std::vector<std::vector<double>> centered(dim, std::vector<double>(n));
    std::vector<double> norm(dim, 0);
    for (int c = 0; c < dim; ++c) {
        double mean = 0;
        for (size_t r = 0; r < n; ++r) mean += rows[r][c];
        mean /= n;
        for (size_t r = 0; r < n; ++r) {
            centered[c][r] = rows[r][c] - mean;
            norm[c] += centered[c][r] * centered[c][r];
        }
        norm[c] = std::sqrt(norm[c]);
    }
    std::vector<std::vector<double>> corr(dim, std::vector<double>(dim, 0));
    for (int c = 0; c < dim; ++c) for (int j = c+1; j < dim; ++j) {
        const double denominator = norm[c] * norm[j];
        double dot = 0;
        if (denominator > 1e-18) {
            for (size_t r = 0; r < n; ++r) dot += centered[c][r] * centered[j][r];
            corr[c][j] = corr[j][c] = std::abs(dot / denominator);
        }
    }
    auto active = result.columns;
    double previous_cv = 0, best_drop = -std::numeric_limits<double>::infinity();
    bool has_best = false;
    for (int step = 0; step < 64 && active.size() > 1; ++step) {
        std::vector<double> scores(active.size(), 0);
        for (size_t c = 0; c < active.size(); ++c) {
            for (int j : active) scores[c] += corr[active[c]][j];
            scores[c] /= active.size()-1;
        }
        double mean = std::accumulate(scores.begin(), scores.end(), 0.0) / scores.size();
        double variance = 0;
        for (double score : scores) variance += (score-mean)*(score-mean);
        const double deviation = std::sqrt(variance / scores.size());
        const double cv = deviation / std::max(std::abs(mean), 1e-18);
        if (step && previous_cv-cv > best_drop+1e-15) {
            best_drop = previous_cv-cv;
            result.columns = active;
            has_best = true;
        }
        previous_cv = cv;
        const double threshold = mean-0.5*deviation;
        std::vector<int> next;
        for (size_t c = 0; c < active.size(); ++c) if (scores[c] >= threshold) next.push_back(active[c]);

        if (next.size() == active.size()) {
            const size_t weakest = std::min_element(scores.begin(), scores.end())-scores.begin();
            next.erase(next.begin()+weakest);
        }
        if (next.empty()) throw std::runtime_error("Selector produced an empty subspace.");
        active = std::move(next);
        result.dimension_path.push_back(static_cast<int>(active.size()));
    }
    if (!has_best) result.columns = active;
    result.best_cv_drop = has_best ? best_drop : 0;
    return result;
}

EndToEndEncoding encode_end_to_end_table(const EndToEndRows& rows,
    const std::string& method, uint32_t seed, int page_size)
{
    const auto start = Clock::now();
    if (rows.empty() || rows[0].empty() || rows[0].size() > 255 || page_size < 1 || page_size > 10000)
        throw std::invalid_argument("Unsupported end-to-end table shape.");
    const bool cluster_method = method == "ACluster" || method == "KCluster";
    const bool baseline = std::find(baseline_methods.begin(), baseline_methods.end(), method) != baseline_methods.end();
    if (method != "ALP" && !cluster_method && !baseline)
        throw std::invalid_argument("Unknown end-to-end method.");
    EndToEndEncoding result;
    const size_t dim = rows[0].size();
    for (const auto& row : rows) {
        if (row.size() != dim) throw std::invalid_argument("Ragged end-to-end table.");
        for (double v : row) if (!std::isfinite(v) || std::abs(v) >= 1e17)
            throw std::invalid_argument("End-to-end input requires finite binary64 values with abs(value)<1e17.");
    }
    auto select_start = Clock::now();
    if (cluster_method) result.selection = select_end_to_end_columns(rows);
    result.metrics.selection_ns = cluster_method ? elapsed(select_start) : 0;
    put(result.bytes, table_magic, 4); put(result.bytes, dim, 4);
    put(result.bytes, rows.size(), 8); put(result.bytes, (rows.size()+page_size-1)/page_size, 4);
    for (size_t begin = 0; begin < rows.size(); begin += page_size) {
        auto page = baseline ? encode_baseline_page(rows, begin, std::min(rows.size(), begin+page_size), method, result.metrics)
            : encode_page(rows, begin, std::min(rows.size(), begin+page_size),
            result.selection.columns, method, seed+static_cast<uint32_t>(begin/page_size), page_size, result.metrics);
        payload(result.bytes, page);
    }
    result.metrics.metadata_bytes = result.bytes.size()-result.metrics.cluster_bytes-
        result.metrics.alp_bytes-result.metrics.baseline_bytes-result.metrics.correction_bytes;
    result.metrics.encode_ns = elapsed(start);
    return result;
}

EndToEndRows decode_end_to_end_table(const std::vector<uint8_t>& bytes) {
    Reader reader{bytes};
    if (reader.get(4) != table_magic) throw std::runtime_error("Invalid table magic/version.");
    const int dim = static_cast<int>(reader.get(4));
    const uint64_t n = reader.get(8), pages = reader.get(4);
    if (dim < 1 || dim > 255 || !n || !pages || pages > n || n > bytes.size()*10000ULL)
        throw std::runtime_error("Invalid table shape.");
    EndToEndRows rows;
    rows.reserve(n);
    for (uint64_t p = 0; p < pages; ++p) {
        auto page = decode_page(reader.block(), dim);
        if (page.size() > n-rows.size()) throw std::runtime_error("Excess decoded rows.");
        rows.insert(rows.end(), std::make_move_iterator(page.begin()), std::make_move_iterator(page.end()));
    }
    if (rows.size() != n) throw std::runtime_error("Decoded table row count mismatch.");
    reader.finish();
    return rows;
}

void verify_end_to_end_table(const EndToEndRows& original, const EndToEndRows& decoded, int page_size) {
    if (original.size() != decoded.size()) throw std::runtime_error("Round-trip row count mismatch.");
    for (size_t begin = 0; begin < original.size(); begin += page_size) {
        const size_t end = std::min(original.size(), begin+page_size);
        std::vector<size_t> a(end-begin), b(end-begin);
        std::iota(a.begin(), a.end(), begin); b = a;
        auto comparator = [](const EndToEndRows& rows, size_t x, size_t y) {
            for (size_t c = 0; c < rows[x].size(); ++c) {
                if (bits(rows[x][c]) != bits(rows[y][c])) return bits(rows[x][c]) < bits(rows[y][c]);
            }
            return false;
        };
        for (size_t r = begin; r < end; ++r) if (original[r].size() != decoded[r].size())
            throw std::runtime_error("Round-trip column count mismatch.");
        std::sort(a.begin(), a.end(), [&](size_t x, size_t y) { return comparator(original, x, y); });
        std::sort(b.begin(), b.end(), [&](size_t x, size_t y) { return comparator(decoded, x, y); });
        for (size_t r = 0; r < a.size(); ++r) for (size_t c = 0; c < original[a[r]].size(); ++c)
            if (bits(original[a[r]][c]) != bits(decoded[b[r]][c]))
                throw std::runtime_error("Exact full-tuple multiset mismatch at page " + std::to_string(begin/page_size));
    }
}

void test_end_to_end_pipeline() {
    EndToEndRows rows;

    for (int r = 0; r < 37; ++r) {
        const double x = (r%7)*0.1;
        rows.push_back({x, 2*x, 3*x, double((r*19)%37), r%3 ? -0.0 : 1e-22});
    }
    rows.push_back(rows[3]);
    const std::vector<std::string> methods{"ACluster", "KCluster", "ALP", "Chimp", "Elf", "Gorilla", "RLE", "Huffman", "REGER", "Vortex"};
    for (const auto& method : methods) {
        std::cout << "Testing " << method << "...\n";
        auto encoded = encode_end_to_end_table(rows, method, 42, 13);
        auto decoded = decode_end_to_end_table(encoded.bytes);
        verify_end_to_end_table(rows, decoded, 13);
        if ((method == "ACluster" || method == "KCluster") && (encoded.selection.columns.empty() || encoded.selection.columns.size() == rows[0].size()))
            throw std::runtime_error("Test failed to exercise split-column path.");
        if (encoded.metrics.cluster_bytes+encoded.metrics.alp_bytes+encoded.metrics.baseline_bytes+encoded.metrics.metadata_bytes+
            encoded.metrics.correction_bytes != encoded.bytes.size()) throw std::runtime_error("Byte accounting mismatch.");
        auto truncated = encoded.bytes; truncated.pop_back();
        bool rejected = false;
        try { decode_end_to_end_table(truncated); } catch (const std::exception&) { rejected = true; }
        if (!rejected) throw std::runtime_error("Truncated stream accepted.");
    }
    for (const EndToEndRows small : {EndToEndRows{{-0.0}}, EndToEndRows{{0.1},{1e-22},{-0.0},{0.1}},
             EndToEndRows{{0.22999999999999998},{0.11764705882352941},{-0.22999999999999998},{1.0000000000000002}},
             EndToEndRows(8, std::vector<double>(4, 7.0))}) {
        for (const auto& method : methods) {
            auto encoded = encode_end_to_end_table(small, method, 42);
            verify_end_to_end_table(small, decode_end_to_end_table(encoded.bytes));
        }
    }

    EndToEndRows wrong = rows;
    std::swap(wrong[0][3], wrong[1][3]);
    bool rejected = false;
    try { verify_end_to_end_table(rows, wrong); } catch (const std::exception&) { rejected = true; }
    if (!rejected) throw std::runtime_error("Verification missed broken tuple associations.");
    std::cout << "End-to-end tests passed: split columns, duplicate tuples, short pages, signed zero,"
                 " precision corrections, constant/1D inputs, byte accounting, truncated streams.\n";
}
