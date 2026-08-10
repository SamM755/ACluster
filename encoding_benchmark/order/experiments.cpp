#include "experiment_hook.h"
#include <algorithm>
#include <array>
#include <cstring>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
namespace fs = std::filesystem;
using Bytes = std::vector<uint8_t>;
using Tuples = std::vector<std::vector<uint64_t>>;
Tuples tuples(const std::vector<std::vector<double>>& rows) {
    Tuples result;
    for (const auto& row : rows) {
        std::vector<uint64_t> tuple;
        for (double x : row) { uint64_t bits; std::memcpy(&bits, &x, 8); tuple.push_back(bits); }
        result.push_back(std::move(tuple));
    }
    std::sort(result.begin(), result.end());
    return result;
}
Tuples expectedTuples(const std::vector<std::string>& lines, int dim) {
    std::vector<std::vector<double>> rows;
    for (const auto& line : lines) {
        std::stringstream ss(line); std::string field; std::vector<double> row;
        while (std::getline(ss, field, ',')) row.push_back(std::stod(field));
        if (int(row.size()) != dim) throw std::runtime_error("Inconsistent dimensions");
        rows.push_back(std::move(row));
    }
    return tuples(rows);
}
void save(const fs::path& path, const Bytes& bytes) {
    std::ofstream out(path, std::ios::binary);
    out.write(reinterpret_cast<const char*>(bytes.data()), bytes.size());
    out.close(); if (!out) throw std::runtime_error("Cannot write output");
}
struct Packed {
    Bytes bytes; size_t bits = 0;
    void append(const Bytes& from, size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            if (bits / 8 == bytes.size()) bytes.push_back(0);
            if ((from[i/8] >> (7-i%8)) & 1) bytes[bits/8] |= uint8_t(1 << (7-bits%8));
            ++bits;
        }
    }
    void pageCount(uint32_t n) { for(int i=0;i<4;++i) bytes[5+i]=uint8_t(n>>(24-8*i)); }
};
struct Totals { long long rows=0, pageBytes=0, initialK=0, k=0, before=0, after=0, deleted=0; Packed shared; long long encode_ns=0,decode_ns=0; };
int run_order_experiment(int argc, char** argv) { try {
    std::cout << std::unitbuf;
    if (argc==2 && std::string(argv[1])=="self-test") { experimentSelfTest(); return 0; }
    if (argc < 3 || argc > 6) throw std::runtime_error("INPUT NEW_OUTPUT [original|page-sort] [both|baseline] [RUNS]");
    fs::path input=argv[1], output=argv[2];
    const std::string mode=argc>3?argv[3]:"original", methods=argc>4?argv[4]:"both";
    if(mode!="original" && mode!="page-sort") throw std::runtime_error("Unknown layout");
    if(methods!="both" && methods!="baseline") throw std::runtime_error("Unknown methods");
    const int runs=argc>5?std::stoi(argv[5]):3;
    if(runs<1) throw std::runtime_error("Invalid repeat count");
    const int treatments=methods=="both"?2:1;
    if(fs::exists(output)) throw std::runtime_error("Use a fresh output folder");
    fs::create_directories(output/"baseline");
    if(treatments==2) fs::create_directories(output/"heuristic");
    std::ifstream source(input); if(!source) throw std::runtime_error("Cannot read input");
    std::ofstream pages(output/"pages.csv");
    pages << "Page,Treatment,Rows,Dim,IndependentBytes,ValidBits,InitialK,K,BeforeBits,AfterBits,Reassigned,Deleted,Checked,Skipped,DeletionQueries,Verified\n";
    std::ofstream timings(output/"timings.csv");
    timings << "Page,Treatment,Run,Rows,EncodeNs,DecodeNs,Verified\n";
    std::array<Totals,2> totals; int dim=0, pageIndex=0; std::string line;
    while(true) {
        std::vector<std::string> page;
        while(page.size()<10000 && std::getline(source,line)) {
            if(!line.empty() && line.back()=='\r') line.pop_back();
            if(!line.empty()) page.push_back(line);
        }
        if(page.empty()) break;
        if(!dim) dim=1+std::count(page[0].begin(),page[0].end(),',');
        const auto expected=expectedTuples(page,dim);
        if(mode=="page-sort") {
            std::vector<double> keys;
            for(const auto& row:page) keys.push_back(std::stod(row));
            std::vector<size_t> order(page.size()); std::iota(order.begin(),order.end(),0);
            std::stable_sort(order.begin(),order.end(),[&](size_t a,size_t b){return keys[a]<keys[b];});
            std::vector<std::string> sorted; sorted.reserve(page.size());
            for(size_t i:order) sorted.push_back(std::move(page[i]));
            page=std::move(sorted);
            if(!std::is_sorted(page.begin(),page.end(),[](const auto& a,const auto& b){return std::stod(a)<std::stod(b);}))
                throw std::runtime_error("Within-page sort failed");
        }
        long long baselineK=0, baselineBits=0;
        for(int treatment=0;treatment<treatments;++treatment) {
            const std::string label=treatment?"heuristic":"baseline";
            experimentBudget=treatment?4:-2;
            ClusterLogic::ExperimentEncodedPage encoded;
            ExperimentStats stats;
            for(int run=-1;run<runs;++run) {
                const auto start=std::chrono::steady_clock::now();
                auto current=ClusterLogic::encodeExperimentPage(page,dim);
                const auto encoded_at=std::chrono::steady_clock::now();
                auto decoded=ClusterLogic::decode_multidim_impl(current.bytes);
                const auto decoded_at=std::chrono::steady_clock::now();
                if(tuples(decoded)!=expected) throw std::runtime_error("Exact binary64 tuple multiset mismatch");
                if(run==-1) { encoded=current; stats=experimentStats; }
                else {
                    if(current.bytes!=encoded.bytes) throw std::runtime_error("Nondeterministic encoding");
                    const auto enc=std::chrono::duration_cast<std::chrono::nanoseconds>(encoded_at-start).count();
                    const auto dec=std::chrono::duration_cast<std::chrono::nanoseconds>(decoded_at-encoded_at).count();
                    totals[treatment].encode_ns+=enc; totals[treatment].decode_ns+=dec;
                    timings << pageIndex << ',' << label << ',' << run+1 << ',' << page.size() << ',' << enc << ',' << dec << ",1\n";
                }
            }
            if(!treatment) { baselineK=stats.k; baselineBits=stats.after_bits; }
            else if(stats.initial_k!=baselineK || stats.before_bits!=baselineBits)
                throw std::runtime_error("Paired baseline mismatch");
            if(encoded.bytes.size()!=(encoded.valid_bits+7)/8 || encoded.valid_bits<88)
                throw std::runtime_error("Page bit accounting");
            auto& t=totals[treatment];
            t.rows+=page.size(); t.pageBytes+=encoded.bytes.size(); t.initialK+=stats.initial_k;
            t.k+=stats.k; t.before+=stats.before_bits; t.after+=stats.after_bits; t.deleted+=stats.deleted;
            t.shared.append(encoded.bytes, pageIndex?88:0, encoded.valid_bits);
            save(output/label/("page_"+std::to_string(pageIndex)+".bin"),encoded.bytes);
            pages << pageIndex << ',' << label << ',' << page.size() << ',' << dim << ',' << encoded.bytes.size()
                << ',' << encoded.valid_bits << ',' << stats.initial_k << ',' << stats.k << ',' << stats.before_bits
                << ',' << stats.after_bits << ',' << stats.reassigned << ',' << stats.deleted << ',' << stats.checked
                << ',' << stats.skipped << ',' << stats.deletion_queries << ",1\n";
            pages.flush();
        }
        ++pageIndex;
        std::cout << input.stem().string() << " page=" << pageIndex << " rows=" << totals[0].rows
            << " baselineK=" << totals[0].k;
        if(treatments==2) std::cout << " heuristicK=" << totals[1].k;
        std::cout << " mode=" << mode << " verified\n";
    }
    if(pageIndex==0) throw std::runtime_error("Empty input");
    std::ofstream summary(output/"summary.csv");
    summary << "Dataset,Treatment,Runs,Rows,Dim,Pages,SharedBytes,SharedCR,IndependentBytes,IndependentCR,InitialK,K,BeforeBits,AfterBits,Deleted,Verified,Mode,EncodeNsPerRow,DecodeNsPerRow\n";
    for(int treatment=0;treatment<treatments;++treatment) {
        auto& t=totals[treatment]; const std::string label=treatment?"heuristic":"baseline";
        t.shared.pageCount(pageIndex); save(output/label/"shared.bin",t.shared.bytes);
        const double raw=double(t.rows)*dim*8;
        summary << input.stem().string() << ',' << label << ',' << runs << ',' << t.rows << ',' << dim << ',' << pageIndex
            << ',' << t.shared.bytes.size() << ',' << std::setprecision(15) << raw/t.shared.bytes.size()
            << ',' << t.pageBytes << ',' << raw/t.pageBytes << ',' << t.initialK << ',' << t.k
            << ',' << t.before << ',' << t.after << ',' << t.deleted << ",1," << mode << ',' << double(t.encode_ns)/runs/t.rows << ',' << double(t.decode_ns)/runs/t.rows << '\n';
        std::cout << "COMPLETE " << input.stem().string() << ' ' << label << " CR=" << std::setprecision(15)
            << raw/t.shared.bytes.size() << " K=" << t.k << " pageCR=" << raw/t.pageBytes << '\n';
    }
    return 0;
} catch(const std::exception& e) { std::cerr << "ERROR " << e.what() << '\n'; return 1; } }
