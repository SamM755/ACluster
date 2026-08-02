#include "benchmark_adapter.h"

#include "algorithms/acluster.h"
#include "algorithms/alp.h"
#include "algorithms/chimp.h"
#include "algorithms/elf.h"
#include "algorithms/gorilla.h"
#include "algorithms/huffman.h"
#include "algorithms/kcluster.h"
#include "algorithms/kcluster_kmeans.h"
#include "algorithms/rle.h"
#include "algorithms/secondary_compressor.h"

#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

extern int NUM_RUNS;
extern std::ofstream resultsFile;
extern std::ofstream memoryResultsFile;
extern std::map<std::string, std::vector<double>> methodMemoryBytes;
extern std::string currentScenario;

void init_results_csv(const std::string& filepath);
void init_memory_csv(const std::string& filepath);
void write_memory_method_average_csv(const std::string& filepath);

template<typename Compressor>
void run_columnar_encoder_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir,
    SecondaryCompressor* sec_comp
);

void run_multidim_reger_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir,
    SecondaryCompressor* sec_comp
);

template<typename ClusterCompressor>
void run_cluster_encoder_test(
    const std::string& dataset_name,
    int dim,
    int num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir,
    SecondaryCompressor* sec_comp,
    int dimension_group_size
);

void run_vortex_test(
    const std::string& dataset_name,
    int expected_dim,
    int expected_num_rows,
    const std::string& csv_file_path,
    const std::string& output_dir
);

void run_cluster_reference_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
);

void run_dimensionality_impact_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
);

void run_subspace_grouping_experiment(
    const std::string& file_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims,
    int dimension_group_size
);

void set_benchmark_runs(int runs) {
    NUM_RUNS = runs;
}

void begin_benchmark_output(const std::string& output_dir, const std::string& scenario, bool memory) {
    std::filesystem::create_directories(output_dir);
    if (resultsFile.is_open()) resultsFile.close();
    if (memoryResultsFile.is_open()) memoryResultsFile.close();
    methodMemoryBytes.clear();
    currentScenario = scenario;
    init_results_csv(output_dir + "/results.csv");
    if (memory) init_memory_csv(output_dir + "/memory.csv");
}

void finish_benchmark_output(const std::string& output_dir, bool memory) {
    if (memory) write_memory_method_average_csv(output_dir + "/memory_average.csv");
    if (resultsFile.is_open()) resultsFile.close();
    if (memoryResultsFile.is_open()) memoryResultsFile.close();
}

void run_baseline_dataset(
    const DatasetConfig& dataset,
    const std::string& input_path,
    const std::string& output_dir,
    bool gzip,
    bool include_vortex
) {
    NoneCompressor none;
    GzipCompressor gzip_compressor;
    SecondaryCompressor* secondary = gzip
        ? static_cast<SecondaryCompressor*>(&gzip_compressor)
        : static_cast<SecondaryCompressor*>(&none);
    run_columnar_encoder_test<Chimp>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_columnar_encoder_test<Elf>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_columnar_encoder_test<Gorilla>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_columnar_encoder_test<Rle>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_columnar_encoder_test<Huffman>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_columnar_encoder_test<Alp>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    run_multidim_reger_test(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary);
    if (include_vortex) run_vortex_test(dataset.name, dataset.dim, dataset.rows, input_path, output_dir);
    run_cluster_encoder_test<KCluster>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary, 0);
    run_cluster_encoder_test<ACluster>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, secondary, 0);
}

void run_kmeans_dataset(
    const DatasetConfig& dataset,
    const std::string& input_path,
    const std::string& output_dir
) {
    NoneCompressor none;
    run_cluster_encoder_test<KCluster>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, &none, 0);
    run_cluster_encoder_test<ACluster>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, &none, 0);
    run_cluster_encoder_test<KClusterKMeans>(dataset.name, dataset.dim, dataset.rows, input_path, output_dir, &none, 0);
}

void run_reference_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    run_cluster_reference_experiment(data_root, output_dir, datasets, dims);
}

void run_dimensionality_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims
) {
    run_dimensionality_impact_experiment(data_root, output_dir, datasets, dims);
}

void run_subspace_core(
    const std::string& data_root,
    const std::string& output_dir,
    const std::vector<std::string>& datasets,
    const std::vector<int>& dims,
    int group_size
) {
    run_subspace_grouping_experiment(data_root, output_dir, datasets, dims, group_size);
}
