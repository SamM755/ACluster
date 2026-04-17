#ifndef KCLUSTER_KMEANS_H
#define KCLUSTER_KMEANS_H

#include <string>
#include <vector>
#include <cstdint>

class KClusterKMeans {
public:
    std::string get_name() const { return "KClusterKMeans"; }

    std::vector<uint8_t> encode_multidim(
        const std::string& csv_file_path,
        int& dim,
        int k,
        int page_size,
        int pack_size = 10,
        int block_size = 10
    );

    std::vector<std::vector<double>> decode_multidim(const std::vector<uint8_t>& compressed_data);
};

#endif
