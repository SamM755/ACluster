#include "pipeline.h"
#include <cstring>
OptimalPageMeasurement measure_optimal_comparison_page(
    const std::vector<std::string>& page, const std::string& method, int k,
    int page_size, int pack_size, int block_size, uint32_t seed,
    const std::vector<std::vector<double>>& sorted_expected,
    const std::string& input_path, std::streamoff input_offset,
    const std::string& output_path)
{
    using Clock = std::chrono::steady_clock;
    int dim = 1;
    const auto start = Clock::now();
    std::vector<std::string> fresh_page;
    if (!input_path.empty()) {
        std::ifstream input(input_path, std::ios::binary);
        if (!input || !input.seekg(input_offset)) throw std::runtime_error("Cannot read benchmark input page.");
        std::string line;
        fresh_page.reserve(page.size());
        while (fresh_page.size() < page.size() && std::getline(input, line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty()) fresh_page.push_back(line);
        }
        if (fresh_page.size() != page.size()) throw std::runtime_error("Truncated benchmark input page.");

    }
    const auto read_at = Clock::now();
    auto data = ClusterLogic::preprocessPageFromStrings(input_path.empty() ? page : fresh_page, dim);
    if (!data.data || data.data->num_rows != static_cast<int>(page.size()) || dim != 1)
        throw std::runtime_error("Optimal benchmark preprocessing changed the input row count/dimension.");
    const auto preprocessed = Clock::now();
    auto encoded = ClusterLogic::encodePreprocessedPage(std::move(data), method,
        k, page_size, pack_size, block_size, seed, true);
    const auto encoded_at = Clock::now();
    if (!output_path.empty()) write_binary_file(output_path, encoded.bytes);
    const auto written_at = Clock::now();
    auto decoded = ClusterLogic::decode_multidim_impl(encoded.bytes);
    const auto decoded_at = Clock::now();
    std::sort(decoded.begin(), decoded.end());
    if (decoded.size() != sorted_expected.size()) throw std::runtime_error("Optimal benchmark decoded row-count mismatch.");
    for (size_t i = 0; i < decoded.size(); ++i) {
        if (decoded[i].size() != 1 || std::memcmp(&decoded[i][0], &sorted_expected[i][0], sizeof(double)))
            throw std::runtime_error("Optimal benchmark failed exact binary64 multiset verification.");
    }
    auto elapsed = [](auto a, auto b) { return std::chrono::duration_cast<std::chrono::nanoseconds>(b-a).count(); };
    OptimalPageMeasurement result;
    result.encode_ns = elapsed(start, written_at);
    result.codec_ns = elapsed(read_at, encoded_at);
    result.read_ns = elapsed(start, read_at);
    result.write_ns = elapsed(encoded_at, written_at);
    result.preprocess_ns = elapsed(read_at, preprocessed);
    result.decode_ns = elapsed(written_at, decoded_at);
    result.encoding = std::move(encoded);
    return result;
}
