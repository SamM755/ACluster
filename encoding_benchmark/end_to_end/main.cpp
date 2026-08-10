#include "experiments.h"
int main(int argc, char** argv) {
    try {
        std::cout << std::unitbuf;
        std::map<std::string, std::string> args;
        for (int i=1; i<argc; i+=2) {
            if (i+1==argc) throw std::runtime_error("Expected an option value");
            args[argv[i]]=argv[i+1];
        }
        const std::set<std::string> allowed{"--experiment","--manifest","--data-root","--output","--runs","--warmups","--dataset","--method"};
        for (const auto& arg:args) if (!allowed.count(arg.first)) throw std::runtime_error("Unknown option: "+arg.first);
        auto get=[&](const std::string& key,const std::string& fallback){return args.count(key)?args.at(key):fallback;};
        const auto experiment=get("--experiment","test");
        if (experiment=="test") { test_end_to_end_pipeline(); test_optimal_1d(); return 0; }
        const auto manifest=get("--manifest","data/selected_columns.csv");
        const auto output=get("--output","output/"+experiment);
        const auto dataset=get("--dataset","all");
        const int runs=std::stoi(get("--runs","10"));
        if (experiment=="end-to-end")
            run_end_to_end_experiment(manifest,get("--data-root","data/raw"),output,runs,dataset,get("--method","all"),std::stoi(get("--warmups","3")));
        else if (experiment=="optimal")
            run_optimal_1d_experiment(manifest,get("--data-root","data/final"),output,runs,100,dataset,10000,10,10);
        else throw std::runtime_error("Unknown experiment");
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
