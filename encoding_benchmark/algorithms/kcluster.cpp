// in algorithms/kcluster.cpp
#include "kcluster.h"
#include "cluster_common.h"
#include "cluster_encoder_logic.h"
#include <functional>

namespace ClusterLogic {

    std::vector<uint8_t> encode_multidim_impl(
        const std::string& csv_file_path,
        int& dim,
        const std::function<ClusteringResult(const SoAData<long long>&, int, int)>& cluster_strategy,
        int k,
        int pack_size,
        int block_size,
        int page_size
    );

    std::vector<std::vector<double>> decode_multidim_impl(
        const std::vector<uint8_t>& compressed_data
    );

    ClusteringResult kMedoidLogCost(const SoAData<long long>& data, int k, int max_iter, double tol, int dim);

} // namespace ClusterLogic

std::pair<std::vector<uint8_t>, GlobalDiagnostics> KCluster::encode_multidim_with_diag(
    const std::string& csv_file_path, int& dim,
    int k, int page_size, int pack_size, int block_size) 
{
    // The strategy needs to be wrapped to match the old signature for the old impl function
    // For the new one, we can call the _with_diag versions directly inside.
    return ClusterLogic::encode_multidim_impl_with_diag(
        csv_file_path, dim, "KCluster", k, pack_size, block_size, page_size
    );
}


std::vector<uint8_t> KCluster::encode_multidim(
    const std::string& csv_file_path, int& dim,
    int k, int page_size, int pack_size, int block_size
) {
    auto kcluster_strategy = [](const SoAData<long long>& data, int k_val, int dim_val) {
        return ClusterLogic::kMedoidLogCost(data, k_val, 2, 1.0e-3, dim_val);
    };

    // original version
    // return ClusterLogic::encode_multidim_impl(
    //     csv_file_path, dim, kcluster_strategy, k, pack_size, block_size, page_size
    // );

    //diagonistic version
    return encode_multidim_with_diag(csv_file_path, dim, k, page_size, pack_size, block_size).first;
}

std::vector<std::vector<double>> KCluster::decode_multidim(const std::vector<uint8_t>& compressed_data) {
    return ClusterLogic::decode_multidim_impl(compressed_data);
}