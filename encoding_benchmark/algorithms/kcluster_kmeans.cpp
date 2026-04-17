#include "kcluster_kmeans.h"
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

    ClusteringResult kMeansEuclidean(const SoAData<long long>& data, int k, int max_iter, double tol, int dim);

}

std::vector<uint8_t> KClusterKMeans::encode_multidim(
    const std::string& csv_file_path, int& dim,
    int k, int page_size, int pack_size, int block_size
) {
    auto kmeans_strategy = [](const SoAData<long long>& data, int k_val, int dim_val) {
        return ClusterLogic::kMeansEuclidean(data, k_val, 20, 1.0e-3, dim_val);
    };
    return ClusterLogic::encode_multidim_impl(
        csv_file_path, dim, kmeans_strategy, k, pack_size, block_size, page_size
    );
}

std::vector<std::vector<double>> KClusterKMeans::decode_multidim(const std::vector<uint8_t>& compressed_data) {
    return ClusterLogic::decode_multidim_impl(compressed_data);
}
