#pragma once
#include "cluster_codec.hpp"

namespace cluster_codec {
Bytes CompressPage(const Bytes &page, uint32_t level);
Bytes DecompressPage(const Bytes &stored, uint64_t expected_size, uint64_t maximum_size);
uint64_t CompressedPageBound(uint64_t size);
}
