// in algorithms/acluster.h
#ifndef ACLUSTER_H
#define ACLUSTER_H

#include <string>
#include <vector>
#include <cstdint>
#include "cluster_common.h"

class ACluster {
public:
    std::string get_name() const { return "ACluster"; }

    std::vector<uint8_t> encode_multidim(
        const std::string& csv_file_path,
        int& dim,
        int k, // k is ignored by ACluster
        int page_size,
        int pack_size = 10,
        int block_size = 10
    );

    std::pair<std::vector<uint8_t>, GlobalDiagnostics> encode_multidim_with_diag(
        const std::string& csv_file_path, int& dim,
        int k, int page_size, int pack_size, int block_size
    );

    std::vector<std::vector<double>> decode_multidim(const std::vector<uint8_t>& compressed_data);
};

#endif // ACLUSTER_H