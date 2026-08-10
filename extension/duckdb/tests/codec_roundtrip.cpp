#include "cluster_codec.hpp"
#include "cluster_compression.hpp"
#include <algorithm>
#include <iostream>
#include <random>
#include <stdexcept>

using namespace cluster_codec;
static void Require(bool condition, const char *message) {
	if (!condition) { throw std::runtime_error(message); }
}
using Tuples = std::vector<std::vector<uint64_t>>;
static Tuples Canonical(const std::vector<Cell> &cells, uint32_t dim) {
	Tuples rows;
	for (size_t i = 0; i < cells.size(); i += dim) {
		std::vector<uint64_t> row;
		for (uint32_t c = 0; c < dim; c++) {
			row.push_back(cells[i+c].valid);
			row.push_back(cells[i+c].valid ? cells[i+c].bits : 0);
		}
		rows.push_back(row);
	}
	std::sort(rows.begin(), rows.end());
	return rows;
}
static uint32_t checks = 0;
static uint32_t encoded_pages[3] = {};
static uint32_t patched_pages = 0;
static void CheckHybrid() {
	uint32_t cases = 0;
	for (uint32_t method : {3U, 4U, 5U}) {
		for (uint32_t count : {1U, 1031U, 10000U}) {
			for (uint32_t special = 0; special < 4; special++) {
				std::vector<Cell> cells;
				for (uint32_t r = 0; r < count; r++) {
					const auto group = double((r * 37) % 7);
					cells.push_back({Bits(group), true});
					cells.push_back({Bits(double(r) / 100.0), true});
					cells.push_back({Bits(group * 2.0), true});
					cells.push_back({Bits(-double(r % 13) / 10.0), true});
				}
				if (special == 1) {
					cells[1].bits = Bits(-0.0);
					cells[3].bits = Bits(0.12345678901234567);
				}
				if (special == 2) {
					cells[0].valid = false;
					cells[1].bits = UINT64_C(0x7ff0000000001234);
				}
				if (special == 3) {
					cells[0].bits = Bits(1e200);
					cells[3].bits = Bits(1e-300);
				}
				Options options;
				options.method = method;
				options.k = 8;
				if (method != 5) { options.selected_columns = {2, 0}; }
				Page page;
				auto bytes = Encode(cells, 4, options, 0, page);
				Require(Canonical(cells, 4) == Canonical(Decode(bytes, 4, page), 4), "Hybrid tuple mismatch");
				auto compressed = CompressPage(bytes, 3);
				Require(DecompressPage(compressed, bytes.size(), uint64_t(count) * 4 * 9 + 36) == bytes,
				        "Hybrid Zstd mismatch");
				if (count == 1031 && special == 0) { Require(page.mode == method, "Hybrid path unused"); }
				auto bad = bytes;
				bad.resize(bad.size() - 1);
				bool rejected = false;
				try { Decode(bad, 4, page); } catch (const std::exception &) { rejected = true; }
				Require(rejected, "Hybrid truncation accepted");
				cases++;
			}
		}
	}
	std::cout << "PASS " << cases << " hybrid/ALP page cases including NULL, NaN payloads, signed zero and RAW substreams\n";
}
static void Check(const std::vector<Cell> &cells, uint32_t dim, uint32_t method, int expected_reason = -1) {
	Options options;
	options.method = method;
	options.k = 8;
	Page page;
	auto bytes = Encode(cells, dim, options, 0, page);
	auto decoded = Decode(bytes, dim, page);
	Require(Canonical(cells, dim) == Canonical(decoded, dim), "Tuple bitwise mismatch");
	for (uint32_t level : {1U, 3U}) {
		auto compressed = CompressPage(bytes, level);
		Require(compressed.size() <= CompressedPageBound(bytes.size()), "Zstd size bound violated");
		Require(DecompressPage(compressed, bytes.size(), uint64_t(page.rows) * dim * 9 + 36) == bytes,
		        "Zstd round trip mismatch");
		compressed[compressed.size() / 2] ^= 1;
		bool rejected = false;
		try { DecompressPage(compressed, bytes.size(), bytes.size()); } catch (const std::exception &) { rejected = true; }
		Require(rejected, "Corrupt Zstd page accepted");
	}
	if (expected_reason >= 0) { Require(page.reason == uint32_t(expected_reason), "Unexpected fallback"); }
	encoded_pages[page.mode]++;
	patched_pages += page.corrections > 0;
	Page again;
	Require(bytes == Encode(cells, dim, options, 0, again), "Non-deterministic page encoding");
	for (size_t cut : {size_t(0), size_t(15), bytes.size() - 1}) {
		auto bad = bytes;
		bad.resize(std::min(cut, bad.size()));
		bool rejected = false;
		try { Decode(bad, dim, page); } catch (const std::exception &) { rejected = true; }
		Require(rejected, "Truncation accepted");
	}
	auto bad = bytes;
	bad[bad.size() / 2] ^= 1;
	bool rejected = false;
	try { Decode(bad, dim, page); } catch (const std::exception &) { rejected = true; }
	Require(rejected, "Corrupt page accepted");
	checks++;
}
int main() {
	try {
		CheckHybrid();
		std::mt19937_64 random(42);
		for (uint32_t method = 0; method < 3; method++) {
			for (uint32_t dim : {1U, 3U, 9U, 255U}) {
				std::vector<Cell> cells;
				for (uint32_t r = 0; r < 321; r++) {
					for (uint32_t c = 0; c < dim; c++) {
						cells.push_back({Bits((int(r % 7) - 3) * 0.25 + c), true});
					}
				}
				Check(cells, dim, method);
			}
			std::vector<Cell> page(10000, {Bits(12.5), true});
			Check(page, 1, method);
			Check({{Bits(-3.75), true}}, 1, method);
			page.resize(1000);
			page[1].bits = Bits(-0.0);
			page[2].bits = Bits(1e-200);
			page[3].bits = Bits(0.12345678901234567);
			Check(page, 1, method);
			for (uint64_t bits : {UINT64_C(0x7ff0000000000000), UINT64_C(0xfff0000000000000),
			                     UINT64_C(0x7ff8000000001234), UINT64_C(0x7ff0000000000001)}) {
				page[0] = {bits, true};
				Check(page, 1, method, method ? 3 : 1);
			}
			page[0] = {Bits(1e200), true};
			Check(page, 1, method, method ? 4 : 1);
			page[0] = {0, false};
			Check(page, 1, method, method ? 2 : 1);
			for (int trial = 0; trial < 12; trial++) {
				page.clear();
				for (int r = 0; r < 173; r++) {
					const double x = int64_t(random() % 20001) / 100.0 - 100.0;
					page.push_back({Bits(x), true});
					page.push_back({Bits(x * 2.0 + 0.01), true});
				}
				Check(page, 2, method);
			}
		}
		Require(encoded_pages[1] > 0 && encoded_pages[2] > 0, "Algorithms never used");
		Require(patched_pages > 0, "Precision correction never used");
		std::cout << "PASS " << checks << " page cases; raw=" << encoded_pages[0]
		          << " acluster=" << encoded_pages[1] << " kcluster=" << encoded_pages[2]
		          << " corrected_pages=" << patched_pages << "\n";
		std::cout << "PASS " << checks * 2 << " Zstd page round trips and corruption checks\n";
	} catch (const std::exception &error) {
		std::cerr << error.what() << "\n";
		return 1;
	}
}
