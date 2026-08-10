#include "cluster_compression.hpp"
#include "zstd.h"
#include <stdexcept>

namespace cluster_codec {
uint64_t CompressedPageBound(uint64_t size) {
	return duckdb_zstd::ZSTD_compressBound(size) + 16;
}

Bytes CompressPage(const Bytes &page, uint32_t level) {
	Bytes stored;
	Put(stored, page.size(), 8);
	stored.resize(8 + duckdb_zstd::ZSTD_compressBound(page.size()));
	const auto size = duckdb_zstd::ZSTD_compress(stored.data() + 8, stored.size() - 8,
	                                          page.data(), page.size(), static_cast<int>(level));
	if (duckdb_zstd::ZSTD_isError(size)) {
		throw std::runtime_error(duckdb_zstd::ZSTD_getErrorName(size));
	}
	stored.resize(8 + size);
	Seal(stored);
	return stored;
}

Bytes DecompressPage(const Bytes &stored, uint64_t expected_size, uint64_t maximum_size) {
	Verify(stored);
	Reader reader(stored);
	const auto size = reader.Get(8);
	if (size < 36 || size != expected_size || size > maximum_size || stored.size() <= 16) {
		throw std::runtime_error("Invalid Zstd cluster page size");
	}
	const auto source = stored.data() + 8;
	const auto source_size = stored.size() - 16;
	if (duckdb_zstd::ZSTD_getFrameContentSize(source, source_size) != size ||
	    duckdb_zstd::ZSTD_findFrameCompressedSize(source, source_size) != source_size) {
		throw std::runtime_error("Invalid Zstd cluster frame");
	}
	Bytes page(size);
	const auto decoded = duckdb_zstd::ZSTD_decompress(page.data(), page.size(), source, source_size);
	if (duckdb_zstd::ZSTD_isError(decoded) || decoded != size) {
		throw std::runtime_error("Cannot decompress Zstd cluster page");
	}
	return page;
}
}
