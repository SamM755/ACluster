#include "acluster.h"
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

    ClusteringResult adaptiveGreedyBasePointSelection(const SoAData<long long>& page_data, int dim);

}

std::pair<std::vector<uint8_t>, GlobalDiagnostics> ACluster::encode_multidim_with_diag(
    const std::string& csv_file_path, int& dim,
    int k,
    int page_size, int pack_size, int block_size)
{
    return ClusterLogic::encode_multidim_impl_with_diag(
        csv_file_path, dim, "ACluster", k, pack_size, block_size, page_size
    );
}

std::vector<uint8_t> ACluster::encode_multidim(
    const std::string& csv_file_path, int& dim,
    int k,
    int page_size, int pack_size, int block_size
) {

    return encode_multidim_with_diag(csv_file_path, dim, k, page_size, pack_size, block_size).first;
}

std::vector<std::vector<double>> ACluster::decode_multidim(const std::vector<uint8_t>& compressed_data) {
    return ClusterLogic::decode_multidim_impl(compressed_data);
}
