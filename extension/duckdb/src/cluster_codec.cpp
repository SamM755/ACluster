#include "cluster_codec.hpp"
#include "cluster_encoder_logic.h"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace cluster_codec {
uint64_t Bits(double value) {
	uint64_t result;
	std::memcpy(&result, &value, sizeof(result));
	return result;
}
double Value(uint64_t bits) {
	double result;
	std::memcpy(&result, &bits, sizeof(result));
	return result;
}
void Put(Bytes &out, uint64_t value, unsigned width) {
	for (unsigned i = 0; i < width; i++) {
		out.push_back(static_cast<uint8_t>(value >> (8 * i)));
	}
}
uint64_t Checksum(const Bytes &bytes, size_t count) {
	uint64_t result = UINT64_C(14695981039346656037);
	for (size_t i = 0; i < count; i++) {
		result = (result ^ bytes[i]) * UINT64_C(1099511628211);
	}
	return result;
}
void Seal(Bytes &bytes) { Put(bytes, Checksum(bytes, bytes.size()), 8); }
void Verify(const Bytes &bytes) {
	if (bytes.size() < 8) {
		throw std::runtime_error("Truncated cluster checksum");
	}
	uint64_t expected = 0;
	for (unsigned i = 0; i < 8; i++) {
		expected |= uint64_t(bytes[bytes.size() - 8 + i]) << (8 * i);
	}
	if (expected != Checksum(bytes, bytes.size() - 8)) {
		throw std::runtime_error("Cluster checksum mismatch");
	}
}
uint64_t Reader::Get(unsigned width) {
	if (width < 1 || width > 8 || width > Remaining()) {
		throw std::runtime_error("Truncated cluster field");
	}
	uint64_t result = 0;
	for (unsigned i = 0; i < width; i++) {
		result |= uint64_t(bytes[pos++]) << (8 * i);
	}
	return result;
}
Bytes Reader::Block(size_t length) {
	if (length > Remaining()) {
		throw std::runtime_error("Invalid cluster block length");
	}
	Bytes result(bytes.begin() + pos, bytes.begin() + pos + length);
	pos += length;
	return result;
}
std::string Reader::String() {
	const auto length = Get(4);
	if (length > 4096) {
		throw std::runtime_error("Cluster column name too long");
	}
	auto data = Block(length);
	return std::string(data.begin(), data.end());
}
void Reader::Finish() const {
	if (Remaining()) {
		throw std::runtime_error("Unexpected trailing cluster bytes");
	}
}
std::string MethodName(uint32_t method) {
	if (method == 0) { return "raw"; }
	if (method == 1) { return "acluster"; }
	if (method == 2) { return "kcluster"; }
	if (method == 3) { return "acluster_alp"; }
	if (method == 4) { return "kcluster_alp"; }
	if (method == 5) { return "alp"; }
	throw std::runtime_error("Unknown cluster method");
}
std::string ReasonName(uint32_t reason) {
	static const char *names[] = {"none", "requested_raw", "null", "non_finite", "numeric_range", "size"};
	if (reason > 5) { throw std::runtime_error("Unknown cluster fallback reason"); }
	return names[reason];
}
Bytes Frame(const Bytes &payload, const Page &page) {
	Bytes out;
	Put(out, page.rows, 4);
	Put(out, page.mode, 4);
	Put(out, page.reason, 4);
	Put(out, page.references, 4);
	Put(out, page.corrections, 4);
	Put(out, payload.size(), 8);
	out.insert(out.end(), payload.begin(), payload.end());
	Seal(out);
	return out;
}
Bytes Encode(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number, Page &page) {
	if (options.method >= 3) { return EncodeHybrid(cells, dim, options, page_number, page); }
	std::vector<uint32_t> row_order;
	return EncodeCluster(cells, dim, options, page_number, page, row_order);
}
Bytes EncodeCluster(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number,
                    Page &page, std::vector<uint32_t> &row_order) {
	if (!dim || dim > 255 || cells.empty() || cells.size() % dim || !options.page_size ||
	    options.page_size > 10000 || cells.size() / dim > options.page_size || options.method > 2 ||
	    !options.k || options.k > 10000) {
		throw std::invalid_argument("Invalid cluster page/options");
	}
	page = Page();
	page.rows = static_cast<uint32_t>(cells.size() / dim);
	row_order.resize(page.rows);
	std::iota(row_order.begin(), row_order.end(), 0);
	Bytes raw;
	raw.reserve(cells.size() * 9);
	page.reason = options.method == 0 ? 1 : 0;
	for (const auto &cell : cells) {
		Put(raw, cell.valid ? 1 : 0, 1);
		Put(raw, cell.valid ? cell.bits : 0, 8);
		if (!page.reason && !cell.valid) { page.reason = 2; }
		if (!page.reason && !std::isfinite(Value(cell.bits))) { page.reason = 3; }
		if (!page.reason && std::abs(Value(cell.bits)) >= 1e17) { page.reason = 4; }
	}
	if (page.reason) { return Frame(raw, page); }
	std::vector<std::string> lines(page.rows);
	char buffer[128];
	for (size_t r = 0; r < page.rows; r++) {
		for (size_t c = 0; c < dim; c++) {
			auto converted = std::to_chars(buffer, buffer + sizeof(buffer), Value(cells[r * dim + c].bits));
			if (converted.ec != std::errc()) { throw std::runtime_error("Cluster numeric conversion failed"); }
			if (c) { lines[r].push_back(','); }
			lines[r].append(buffer, converted.ptr);
		}
	}
	int actual_dim = static_cast<int>(dim);
	PreprocessorResult prepared;
	try {
		prepared = ClusterLogic::preprocessPageFromStrings(lines, actual_dim);
	} catch (const std::out_of_range &) {
		page.reason = 4;
		return Frame(raw, page);
	}
	if (!prepared.data || actual_dim != static_cast<int>(dim) || prepared.data->num_rows != int(page.rows)) {
		throw std::runtime_error("Cluster preprocessing changed the page shape");
	}
	std::vector<std::vector<std::pair<uint32_t, uint64_t>>> corrections(page.rows);
	for (uint32_t c = 0; c < dim; c++) {
		const auto scale = std::pow(10.0, prepared.max_decimal_places[c]);
		for (uint32_t r = 0; r < page.rows; r++) {
			const auto normalized = prepared.data->columns[c][r];
			if (normalized < 0 || normalized > INT64_C(2000000000000000000)) {
				page.reason = 4;
				page.corrections = 0;
				return Frame(raw, page);
			}
			const auto restored = static_cast<double>(normalized + prepared.min_values[c]) / scale;
			if (Bits(restored) != cells[size_t(r) * dim + c].bits) {
				corrections[r].emplace_back(c, cells[size_t(r) * dim + c].bits);
				page.corrections++;
			}
		}
	}
	auto encoded = ClusterLogic::encodePreprocessedPage(std::move(prepared),
	    options.method == 1 ? "ACluster" : "KCluster", std::min(options.k, page.rows), options.page_size,
	    10, 10, options.seed + static_cast<uint32_t>(page_number));
	Bytes payload;
	Put(payload, encoded.bytes.size(), 8);
	payload.insert(payload.end(), encoded.bytes.begin(), encoded.bytes.end());
	for (uint32_t r = 0; r < page.rows; r++) {
		for (const auto &correction : corrections[encoded.row_order[r]]) {
			Put(payload, r * dim + correction.first, 4);
			Put(payload, correction.second, 8);
		}
	}
	if (payload.size() >= raw.size()) {
		page.reason = 5;
		page.corrections = 0;
		return Frame(raw, page);
	}
	page.mode = options.method;
	page.references = static_cast<uint32_t>(encoded.reference_count);
	row_order.assign(encoded.row_order.begin(), encoded.row_order.end());
	return Frame(payload, page);
}
std::vector<Cell> Decode(const Bytes &bytes, uint32_t dim, const Page &page) {
	if (page.mode >= 3) { return DecodeHybrid(bytes, dim, page); }
	if (!dim || dim > 255 || !page.rows || page.rows > 10000 || bytes.size() > uint64_t(page.rows) * dim * 9 + 36) {
		throw std::runtime_error("Invalid cluster page shape/size");
	}
	Verify(bytes);
	Reader reader(bytes);
	if (reader.Get(4) != page.rows || reader.Get(4) != page.mode || reader.Get(4) != page.reason ||
	    reader.Get(4) != page.references || reader.Get(4) != page.corrections) {
		throw std::runtime_error("Cluster page directory mismatch");
	}
	const auto payload_length = reader.Get(8);
	if (reader.Remaining() < 8 || payload_length != reader.Remaining() - 8) {
		throw std::runtime_error("Invalid cluster page payload length");
	}
	std::vector<Cell> result(size_t(page.rows) * dim);
	if (page.mode == 0) {
		if (payload_length != result.size() * 9 || page.references || page.corrections || !page.reason || page.reason > 5) {
			throw std::runtime_error("Invalid RAW cluster page");
		}
		for (auto &cell : result) {
			const auto valid = reader.Get(1);
			if (valid > 1) { throw std::runtime_error("Invalid cluster validity byte"); }
			cell.valid = valid == 1;
			cell.bits = reader.Get(8);
		}
	} else {
		if (page.mode > 2 || page.reason || !page.references || page.references > page.rows || page.corrections > result.size()) {
			throw std::runtime_error("Invalid encoded cluster page");
		}
		auto inner = reader.Block(reader.Get(8));
		ClusterLogic::BitStreamReader probe(inner);
		if (probe.readBits(8) != dim || probe.readBits(16) != 10 || probe.readBits(16) != 10 ||
		    probe.readBits(32) != 1 || probe.readBits(16) < page.rows || probe.readBits(16) != page.references ||
		    probe.readBits(16) != page.rows) {
			throw std::runtime_error("Invalid embedded cluster header");
		}
		auto rows = ClusterLogic::decode_multidim_impl(inner);
		if (rows.size() != page.rows) { throw std::runtime_error("Cluster decoded row count mismatch"); }
		for (size_t r = 0; r < rows.size(); r++) {
			if (rows[r].size() != dim) { throw std::runtime_error("Cluster decoded dimension mismatch"); }
			for (size_t c = 0; c < dim; c++) { result[r * dim + c].bits = Bits(rows[r][c]); }
		}
		uint64_t previous = 0;
		for (uint32_t i = 0; i < page.corrections; i++) {
			const auto position = reader.Get(4);
			if (position >= result.size() || (i && position <= previous)) {
				throw std::runtime_error("Invalid cluster correction position");
			}
			result[position].bits = reader.Get(8);
			previous = position;
		}
	}
	reader.Get(8);
	reader.Finish();
	return result;
}
}
