#ifndef KCLUSTER_SUB_H
#define KCLUSTER_SUB_H

#include <string>
#include <vector>
#include <cstdint>

class KClusterSub {
public:
    explicit KClusterSub(int dimension_group_size = 4) : dimension_group_size_(dimension_group_size) {}
    void set_dimension_group_size(int v) { dimension_group_size_ = (v > 0 ? v : 1); }

    std::string get_name() const { return "KClusterSub"; }

    std::vector<uint8_t> encode_multidim(
        const std::string& csv_file_path,
        int& dim,
        int k,
        int page_size,
        int pack_size = 10,
        int block_size = 10
    );

    std::vector<std::vector<double>> decode_multidim(const std::vector<uint8_t>& compressed_data);

private:
    int dimension_group_size_;
};

#endif
