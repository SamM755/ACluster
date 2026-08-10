#include "experiments.h"
void run_end_to_end_experiment(const std::string& manifest, const std::string& file_root,
    const std::string& output_dir, int runs, const std::string& requested_dataset,
    const std::string& requested_method, int warmups)
{
    if (runs < 1 || warmups < 0) throw std::invalid_argument("Invalid measurement/warmup counts.");
    std::vector<std::string> methods{"ACluster", "KCluster", "ALP"};
    const std::vector<std::string> baselines{"Chimp", "Elf", "Gorilla", "RLE", "Huffman", "REGER", "Vortex"};
    if (requested_method == "baselines") methods = baselines;
    else if (requested_method == "all") methods.insert(methods.end(), baselines.begin(), baselines.end());
    else {
        methods.insert(methods.end(), baselines.begin(), baselines.end());
        if (std::find(methods.begin(), methods.end(), requested_method) == methods.end())
            throw std::invalid_argument("Unsupported end-to-end method.");
        methods = {requested_method};
    }
    std::ifstream input(manifest);
    if (!input) throw std::runtime_error("Cannot open end-to-end manifest.");
    auto split = [](const std::string& line, char delimiter) {
        std::vector<std::string> fields;
        std::stringstream stream(line);
        std::string field;
        while (std::getline(stream, field, delimiter)) fields.push_back(field);
        return fields;
    };
    std::string line;
    std::getline(input, line);
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line != "Dataset,Rows,OriginalDim,SelectedDim,SelectedColumns")
        throw std::runtime_error("Unexpected end-to-end manifest header.");
    bool found = false;
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const auto fields = split(line, ',');
        if (fields.size() != 5) throw std::runtime_error("Malformed end-to-end manifest row.");
        const auto& name = fields[0];
        if (std::stoi(fields[2]) <= 1) continue;
        if (requested_dataset != "all" && requested_dataset != name) continue;
        found = true;
        const size_t expected_rows = std::stoull(fields[1]);
        const int dim = std::stoi(fields[2]);
        std::vector<int> expected_columns;
        for (const auto& c : split(fields[4], ';')) expected_columns.push_back(std::stoi(c));
        const auto path = std::filesystem::path(file_root) / (name+".csv");
        auto load_start = std::chrono::steady_clock::now();
        const auto rows = read_end_to_end_table(path.string(), expected_rows, dim);
        const double load_ms = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now()-load_start).count();
        const auto dir = std::filesystem::path(output_dir) / name;
        std::filesystem::create_directories(dir);
        const auto runs_path = dir / (requested_method+"_runs.csv");
        if (std::filesystem::exists(runs_path))
            throw std::runtime_error("Result already exists; choose a new output directory: " + runs_path.string());
        const auto warmup_path = dir / (requested_method+"_warmup_runs.csv");
        if (warmups && std::filesystem::exists(warmup_path))
            throw std::runtime_error("Warmup result already exists; choose a new output directory.");
        std::ofstream measurement_output(runs_path), warmup_output;
        if (warmups) warmup_output.open(warmup_path);
        if (!measurement_output || (warmups && !warmup_output))
            throw std::runtime_error("Cannot create end-to-end results.");
        const std::string header = "Dataset,Method,Run,Rows,OriginalDim,SelectedDim,PageSize,Pages,Seed,RawBytes,TotalBytes,"
                  "ClusterBytes,ALPBytes,CorrectionBytes,MetadataBytes,CorrectionValues,References,CR,"
                  "SelectionMs,PreprocessingMs,ClusterMs,ALPMs,OtherEncodeMs,EncodeMs,DecodeMs,"
                  "EncodeNsPerRow,DecodeNsPerRow,LoadMs,Verified,SelectedColumns,DimensionPath,BaselineBytes,BaselineMs\n";
        measurement_output << header << std::setprecision(12);
        if (warmups) warmup_output << header << std::setprecision(12);
        std::cout << "End-to-end " << name << ": " << expected_rows << " rows x " << dim << " columns.\n";
        for (int run = -warmups; run < runs; ++run) {
            const bool is_warmup = run < 0;
            const int phase_run = is_warmup ? run+warmups : run;
            auto& output = is_warmup ? warmup_output : measurement_output;

            for (size_t index = 0; index < methods.size(); ++index) {
                const auto& method = methods[(index+phase_run)%methods.size()];
                std::cout << "  " << method << (is_warmup ? " warmup " : " run ")
                          << phase_run+1 << '/' << (is_warmup ? warmups : runs) << " encoding...\n";
                auto encoded = encode_end_to_end_table(rows, method);
                if ((method == "ACluster" || method == "KCluster") && encoded.selection.columns != expected_columns)
                    throw std::runtime_error("C++ selector differs from retained Python columns for " + name);
                auto decode_start = std::chrono::steady_clock::now();
                const auto decoded = decode_end_to_end_table(encoded.bytes);
                encoded.metrics.decode_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                    std::chrono::steady_clock::now()-decode_start).count();
                verify_end_to_end_table(rows, decoded);
                const auto& m = encoded.metrics;
                const uint64_t raw_bytes = uint64_t(expected_rows)*dim*8;
                if (m.cluster_bytes+m.alp_bytes+m.baseline_bytes+m.correction_bytes+m.metadata_bytes != encoded.bytes.size())
                    throw std::runtime_error("Serialized byte accounting mismatch.");
                if (run == 0) {
                    const auto binary_path = dir / (method+".etb");
                    if (std::filesystem::exists(binary_path)) throw std::runtime_error("Refusing to replace saved stream.");
                    write_binary_file(binary_path.string(), encoded.bytes);
                    if (std::filesystem::file_size(binary_path) != encoded.bytes.size())
                        throw std::runtime_error("Saved stream size mismatch.");
                }
                const auto join = [](const std::vector<int>& values) {
                    std::ostringstream s;
                    for (size_t i = 0; i < values.size(); ++i) { if (i) s << ';'; s << values[i]; }
                    return s.str();
                };
                const auto other_ns = m.encode_ns-m.selection_ns-m.preprocessing_ns-m.cluster_ns-m.alp_ns-m.baseline_ns;
                output << name << ',' << method << ',' << phase_run+1 << ',' << expected_rows << ',' << dim << ','
                    << encoded.selection.columns.size() << ",10000," << m.pages << ",20260929," << raw_bytes << ','
                    << encoded.bytes.size() << ',' << m.cluster_bytes << ',' << m.alp_bytes << ',' << m.correction_bytes << ','
                    << m.metadata_bytes << ',' << m.correction_values << ',' << m.references << ','
                    << double(raw_bytes)/encoded.bytes.size() << ',' << m.selection_ns/1e6 << ','
                    << m.preprocessing_ns/1e6 << ',' << m.cluster_ns/1e6 << ',' << m.alp_ns/1e6 << ','
                    << other_ns/1e6 << ',' << m.encode_ns/1e6 << ',' << m.decode_ns/1e6 << ','
                    << double(m.encode_ns)/expected_rows << ',' << double(m.decode_ns)/expected_rows << ','
                    << load_ms << ",true," << join(encoded.selection.columns) << ',' << join(encoded.selection.dimension_path)
                    << ',' << m.baseline_bytes << ',' << m.baseline_ns/1e6 << '\n';
                output.flush();
                if (!output) throw std::runtime_error("Cannot persist end-to-end measurement.");
                std::cout << "    verified complete tuples; CR=" << double(raw_bytes)/encoded.bytes.size()
                          << ", bytes=" << encoded.bytes.size() << ", encode=" << double(m.encode_ns)/expected_rows
                          << " ns/row, decode=" << double(m.decode_ns)/expected_rows << " ns/row, corrections=" << m.correction_values << '\n';
            }
        }
    }
    if (!found) throw std::invalid_argument("No matching end-to-end dataset.");
}
