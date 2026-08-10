#ifndef CLUSTER_ENCODER_LOGIC_H
#define CLUSTER_ENCODER_LOGIC_H

#include "cluster_common.h"
#include <functional>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace ClusterLogic {

PreprocessorResult preprocessPageFromStrings(const std::vector<std::string>& csv_page, int& dim_ref);

ClusteringResult sortClusteringBySize(std::vector<std::vector<long long>>& medoids,
    std::vector<int>& assignment, std::vector<long long>& sizes);

class BitStreamReader {
public:
    explicit BitStreamReader(const std::vector<uint8_t>& data);
    long long readBits(int numBits);
    size_t position() const { return currentBitPos; }

private:
    const std::vector<uint8_t>& buffer;
    const size_t bufferSize;
    size_t currentBitPos;
};

inline long long zigzagDecode(long long n) { return (n >> 1) ^ (-(n & 1)); }

std::vector<uint8_t> encode_multidim_impl(
    const std::string& csv_file_path,
    int& dim,
    const std::function<ClusteringResult(const SoAData<long long>&, int, int)>& cluster_strategy,
    int k,
    int pack_size,
    int block_size,
    int page_size
);

enum class KMeansReferenceQuantization {
    RoundTwoDecimals,
    TruncateInteger
};

std::vector<uint8_t> encode_multidim_kmeans_quantized_impl(
    const std::string& csv_file_path,
    int& dim,
    int k,
    int pack_size,
    int block_size,
    int page_size,
    KMeansReferenceQuantization quantization
);

std::pair<std::vector<uint8_t>, GlobalDiagnostics> encode_multidim_impl_with_diag(
    const std::string& csv_file_path, int& dim,
    const std::string& algorithm_name,
    int k, int pack_size, int block_size, int page_size,
    int adjacent_reference = 1,
    int kcluster_max_iterations = 10
);

std::vector<std::vector<double>> decode_multidim_impl(
    const std::vector<uint8_t>& compressed_data
);

struct KMedoidResult {
    ClusteringResult result;
    long long final_cost;
};

struct AClusterResult {
    ClusteringResult result;
    long long final_cost;
};

struct ReferenceRedundancyDiagnostics {
    long long reference_count = 0;
    long long individually_removable_reference_count = 0;
    long long reference_cost_bits = 0;
    long long residual_cost_bits = 0;
    long long removable_reference_cost_bits = 0;
};

struct AOptimalResult {
    ClusteringResult result;
    long long final_cost;
};

struct ReferenceCandidateEncoding {
    std::vector<uint8_t> bytes;
    long long residual_objective_bits = 0;
    long long reference_objective_bits = 0;
    int reference_count = 0;
};
ReferenceCandidateEncoding encodePageWithReferences(
    const PreprocessorResult& preprocessed,
    std::vector<std::vector<long long>> references,
    int page_size, int pack_size, int block_size
);

KMedoidResult kMedoidLogCost_with_diag(const SoAData<long long>& data, int k, int max_iter, double tol, int dim,
                                    const uint32_t* seed = nullptr);

struct ClusterPageEncoding {
    std::vector<uint8_t> bytes;
    std::vector<size_t> row_order;
    size_t reference_count = 0;
    ComponentSizes sizes;
    long long residual_objective_bits = -1;
    long long reference_objective_bits = -1;
    long long selection_ns = 0, diagnostics_ns = 0, serialization_ns = 0;
    long long optimal_sort_ns = 0, optimal_prefix_ns = 0, optimal_search_ns = 0, optimal_dp_ns = 0;
    int optimal_dp_states = 0;
    uint64_t optimal_reference_candidates = 0;
};
ClusterPageEncoding encodePreprocessedPage(
    PreprocessorResult preprocessed, const std::string& method,
    int k, int page_size, int pack_size, int block_size, uint32_t seed,
    bool measure_objective = false);
AClusterResult adaptiveGreedyBasePointSelection_with_diag(const SoAData<long long>& page_data, int dim, int adjacent_reference = 1);
ReferenceRedundancyDiagnostics analyzeReferenceRedundancy(
    const SoAData<long long>& page_data,
    const ClusteringResult& clustering_result,
    int dim
);
AOptimalResult adaptiveOptimal1D_with_diag(const SoAData<long long>& page_data, int dim);

ClusteringResult kMedoidLogCost(
    const SoAData<long long>& data,
    int k,
    int max_iter,
    double tol,
    int dim
);

ClusteringResult adaptiveGreedyBasePointSelection(
    const SoAData<long long>& page_data,
    int dim,
    int adjacent_reference = 1
);

ClusteringResult adaptiveOptimal1D(
    const SoAData<long long>& page_data,
    int dim
);

ClusteringResult kMeansEuclidean(
    const SoAData<long long>& data,
    int k,
    int max_iter,
    double tol,
    int dim
);

}

#endif
