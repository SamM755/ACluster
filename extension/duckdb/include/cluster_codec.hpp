#pragma once
#include <cstdint>
#include <string>
#include <vector>

namespace cluster_codec {
using Bytes = std::vector<uint8_t>;
struct Cell {
	uint64_t bits = 0;
	bool valid = true;
};
struct Options {
	uint32_t method = 1;
	std::vector<uint32_t> selected_columns;
	uint32_t page_size = 10000;
	uint32_t k = 100;
	uint32_t seed = 42;
	uint32_t compression = 0;
	uint32_t compression_level = 3;
};
struct Page {
	uint64_t offset = 0;
	uint64_t length = 0;
	uint64_t uncompressed_length = 0;
	uint32_t rows = 0;
	uint32_t mode = 0;
	uint32_t reason = 0;
	uint32_t references = 0;
	uint32_t corrections = 0;
};
uint64_t Bits(double value);
double Value(uint64_t bits);
void Put(Bytes &out, uint64_t value, unsigned width);
uint64_t Checksum(const Bytes &bytes, size_t count);
void Seal(Bytes &bytes);
void Verify(const Bytes &bytes);
class Reader {
public:
	explicit Reader(const Bytes &bytes) : bytes(bytes) {}
	uint64_t Get(unsigned width);
	Bytes Block(size_t length);
	std::string String();
	void Finish() const;
	size_t Remaining() const { return bytes.size() - pos; }
private:
	const Bytes &bytes;
	size_t pos = 0;
};
Bytes Encode(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number, Page &page);
Bytes EncodeCluster(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number,
                    Page &page, std::vector<uint32_t> &row_order);
Bytes Frame(const Bytes &payload, const Page &page);
Bytes EncodeHybrid(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number, Page &page);
std::vector<Cell> DecodeHybrid(const Bytes &bytes, uint32_t dim, const Page &page);
std::vector<Cell> Decode(const Bytes &bytes, uint32_t dim, const Page &page);
std::string MethodName(uint32_t method);
std::string ReasonName(uint32_t reason);
}
