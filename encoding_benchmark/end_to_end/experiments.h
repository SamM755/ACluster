#pragma once
#include "pipeline.h"
struct DatasetEntry { std::string name; int dim, rows; };
std::vector<DatasetEntry> load_selected_dataset_manifest(const std::string&, const std::string&);
void run_optimal_1d_experiment(const std::string&, const std::string&, const std::string&, int, int, const std::string&, int, int, int);
void run_end_to_end_experiment(const std::string&, const std::string&, const std::string&, int, const std::string&, const std::string&, int);
