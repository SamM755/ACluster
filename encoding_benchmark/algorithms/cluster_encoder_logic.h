// in algorithms/cluster_encoder_logic.h

#ifndef CLUSTER_ENCODER_LOGIC_H
#define CLUSTER_ENCODER_LOGIC_H

#include "cluster_common.h"
#include <functional>
#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

namespace ClusterLogic {

// --- Declare BitStreamReader class ---
class BitStreamReader {
public:
    explicit BitStreamReader(const std::vector<uint8_t>& data);
    long long readBits(int numBits);

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

std::pair<std::vector<uint8_t>, GlobalDiagnostics> encode_multidim_impl_with_diag(
    const std::string& csv_file_path, int& dim,
    const std::string& algorithm_name,
    int k, int pack_size, int block_size, int page_size
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

KMedoidResult kMedoidLogCost_with_diag(const SoAData<long long>& data, int k, int max_iter, double tol, int dim);
AClusterResult adaptiveGreedyBasePointSelection_with_diag(const SoAData<long long>& page_data, int dim);


ClusteringResult kMedoidLogCost(
    const SoAData<long long>& data, 
    int k, 
    int max_iter, 
    double tol, 
    int dim
);

ClusteringResult adaptiveGreedyBasePointSelection(
    const SoAData<long long>& page_data, 
    int dim
);

} // namespace ClusterLogic

#endif // CLUSTER_ENCODER_LOGIC_H