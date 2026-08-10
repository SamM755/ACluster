#include "experiments.h"
void run_optimal_1d_experiment(const std::string& manifest, const std::string& file_root,
    const std::string& output_dir, int runs, int k, const std::string& dataset,
    int page_size, int pack_size, int block_size)
{
    if (runs < 1 || k < 1 || page_size < 1 || page_size > 10000 || pack_size < 1 ||
        pack_size > 65535 || block_size < 1 || block_size > 65535)
        throw std::invalid_argument("Invalid optimal-1d parameters.");
    auto configs = load_selected_dataset_manifest(manifest, "optimal-1d");
    configs.erase(std::remove_if(configs.begin(), configs.end(), [&](const auto& c) {
        return c.dim != 1 || (dataset == "all" ? c.name != "SSD" && c.name != "BC" && c.name != "PI" : c.name != dataset);
    }), configs.end());
    if (configs.empty()) throw std::invalid_argument("No matching 1D datasets.");
    const std::filesystem::path out(output_dir);
    if (std::filesystem::exists(out / "pages.csv")) throw std::runtime_error("Use a fresh optimal result directory.");
    std::filesystem::create_directories(out);
    std::ofstream pages(out / "pages.csv"), totals(out / "runs.csv"), protocol(out / "protocol.txt");
    if (!pages || !totals || !protocol) throw std::runtime_error("Cannot create optimal results.");
    protocol << "KOptimal,KCluster,AOptimal,ACluster; independent pages of at most 10000 rows.\n"
        << "runs=" << runs << "; K=" << k << "; page_size=" << page_size << "; pack_size=" << pack_size << "; block_size=" << block_size << "\n"
        << "Encoding: CSV read, preprocessing, selection, diagnostics, serialization and cached file write.\n"
        << "Decoding excludes verification. Every page is checked as a binary64 multiset.\n"
        << "No warmups. Rotating method order. Seed=20261001+page.\n";
    protocol.flush();
    pages << "Dataset,Method,Run,Page,Rows,DistinctValues,References,Bytes,ResidualBits,ReferenceBits,EncodeNs,DecodeNs,PreprocessNs,"
          << "CodecNs,ReadNs,WriteNs,SelectionNs,DiagnosticsNs,SerializationNs,OptimalSortNs,OptimalPrefixNs,OptimalSearchNs,OptimalDPNs,DPStates,ReferenceCandidates,Verified\n";
    totals << "Dataset,Method,Run,Rows,Pages,RequestedK,References,Bytes,CR,EncodeNsPerRow,DecodeNsPerRow,ResidualBits,ReferenceBits,"
           << "CodecNsPerRow,ReadNsPerRow,WriteNsPerRow,PreprocessNsPerRow,SelectionNsPerRow,DiagnosticsNsPerRow,SerializationNsPerRow,"
           << "OptimalSortNsPerRow,OptimalPrefixNsPerRow,OptimalSearchNsPerRow,OptimalDPNsPerRow,Verified\n";
    totals << std::setprecision(12);
    const std::vector<std::string> methods{"KOptimal","KCluster","AOptimal","ACluster"};
    for (const auto& config : configs) {
        const auto input = std::filesystem::path(file_root) / (config.name + ".csv");
        std::ifstream source(input, std::ios::binary);
        if (!source) throw std::runtime_error("Cannot read " + input.string());
        std::vector<std::vector<std::string>> page_lines;
        std::vector<std::vector<std::vector<double>>> expected;
        std::vector<size_t> distinct;
        std::vector<std::streamoff> input_offsets;
        std::string line;
        size_t rows = 0;
        while (true) {
            const auto offset = source.tellg();
            if (!std::getline(source, line)) break;
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (rows % page_size == 0) {
                page_lines.emplace_back(); expected.emplace_back(); input_offsets.push_back(offset);
            }
            size_t consumed = 0;
            const double value = std::stod(line, &consumed);
            if (consumed != line.size() || !std::isfinite(value))
                throw std::runtime_error("Invalid scalar input in " + config.name);
            page_lines.back().push_back(line); expected.back().push_back({value}); ++rows;
        }
        if (rows != static_cast<size_t>(config.rows)) throw std::runtime_error("Manifest row mismatch: " + config.name);
        for (auto& page : expected) {
            std::sort(page.begin(), page.end());
            distinct.push_back(std::set<std::vector<double>>(page.begin(), page.end()).size());
        }
        for (int run = 0; run < runs; ++run) {
            for (size_t mi = 0; mi < methods.size(); ++mi) {
                const auto& method = methods[(mi + run) % methods.size()];
                const bool is_optimal = method == "KOptimal" || method == "AOptimal";
                const auto implementation = is_optimal ? method : method+"1D";
                uint64_t bytes = 0, refs = 0;
                long long enc_ns = 0, dec_ns = 0, residual = 0, reference = 0;
                std::vector<long long> stage_totals(11, 0);
                const auto stream_dir = out / "streams" / config.name / method;
                std::filesystem::create_directories(stream_dir);
                for (size_t p = 0; p < page_lines.size(); ++p) {
                    std::cout << config.name << ' ' << method << " run " << run+1 << '/' << runs
                              << " page " << p+1 << '/' << page_lines.size() << " (u=" << distinct[p] << ")\n";
                    auto measured = measure_optimal_comparison_page(page_lines[p], implementation, k, page_size,
                        pack_size, block_size, 20261001u + static_cast<uint32_t>(p), expected[p],
                        input.string(), input_offsets[p], (stream_dir / ("page_" + std::to_string(p) + ".bin")).string());
                    const auto& e = measured.encoding;
                    const std::vector<long long> stages{measured.codec_ns, measured.read_ns, measured.write_ns,
                        measured.preprocess_ns, e.selection_ns, e.diagnostics_ns, e.serialization_ns,
                        e.optimal_sort_ns, e.optimal_prefix_ns, e.optimal_search_ns, e.optimal_dp_ns};
                    for (size_t s = 0; s < stages.size(); ++s) stage_totals[s] += stages[s];
                    if (measured.encode_ns != measured.read_ns + measured.codec_ns + measured.write_ns ||
                        (is_optimal && e.optimal_dp_states != static_cast<int>(distinct[p])))
                        throw std::runtime_error("Phase or state accounting failure.");
                    bytes += e.bytes.size(); refs += e.reference_count;
                    enc_ns += measured.encode_ns; dec_ns += measured.decode_ns;
                    residual += e.residual_objective_bits; reference += e.reference_objective_bits;
                    pages << config.name << ',' << method << ',' << run+1 << ',' << p << ',' << page_lines[p].size()
                          << ',' << distinct[p] << ',' << e.reference_count << ',' << e.bytes.size() << ','
                          << e.residual_objective_bits << ',' << e.reference_objective_bits << ',' << measured.encode_ns
                          << ',' << measured.decode_ns << ',' << measured.preprocess_ns
                          << ',' << measured.codec_ns << ',' << measured.read_ns << ',' << measured.write_ns
                          << ',' << e.selection_ns << ',' << e.diagnostics_ns << ',' << e.serialization_ns
                          << ',' << e.optimal_sort_ns << ',' << e.optimal_prefix_ns << ',' << e.optimal_search_ns << ',' << e.optimal_dp_ns
                          << ',' << e.optimal_dp_states << ',' << e.optimal_reference_candidates
                          << ",true\n";
                    pages.flush();
                    std::cout << "  page encoded: " << measured.encode_ns / 1e9 << " s; search="
                              << e.optimal_search_ns / 1e9 << " s; DP=" << e.optimal_dp_ns / 1e9 << " s\n";
                }
                totals << config.name << ',' << method << ',' << run+1 << ',' << rows << ',' << page_lines.size()
                       << ',' << k << ',' << refs << ',' << bytes << ',' << double(rows)*8/bytes
                       << ',' << double(enc_ns)/rows << ',' << double(dec_ns)/rows << ',' << residual << ',' << reference;
                for (auto ns : stage_totals) totals << ',' << double(ns)/rows;
                totals << ",true\n";
                totals.flush();
                if (!pages || !totals) throw std::runtime_error("Failed to persist optimal measurements.");
                std::cout << "  Verified CR=" << double(rows)*8/bytes << ", enc=" << double(enc_ns)/rows
                          << " ns/row, dec=" << double(dec_ns)/rows << " ns/row\n";
            }
        }
    }
    std::ofstream(out / "COMPLETE") << "All requested dataset/method/repeat/page round trips verified.\n";
}
