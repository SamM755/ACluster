#include "cluster_codec.hpp"
#include "alp.h"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace cluster_codec {
static std::vector<uint32_t> RemainingColumns(uint32_t dim, const std::vector<uint32_t> &selected) {
	std::vector<bool> seen(dim, false);
	for (auto c : selected) {
		if (c >= dim || seen[c]) { throw std::invalid_argument("Invalid or duplicate cluster column"); }
		seen[c] = true;
	}
	std::vector<uint32_t> remaining;
	for (uint32_t c = 0; c < dim; c++) { if (!seen[c]) { remaining.push_back(c); } }
	return remaining;
}

Bytes EncodeHybrid(const std::vector<Cell> &cells, uint32_t dim, const Options &options, uint64_t page_number, Page &page) {
	if (!dim || dim > 255 || cells.empty() || cells.size() % dim || !options.page_size ||
	    options.page_size > 10000 || cells.size() / dim > options.page_size || options.method < 3 ||
	    options.method > 5 || !options.k || options.k > 10000 ||
	    (options.method == 5 ? !options.selected_columns.empty() : options.selected_columns.empty())) {
		throw std::invalid_argument("Invalid hybrid page/options");
	}
	const auto &selected = options.selected_columns;
	const auto remaining = RemainingColumns(dim, selected);
	page = Page();
	page.rows = static_cast<uint32_t>(cells.size() / dim);
	page.mode = options.method;
	std::vector<uint32_t> order(page.rows);
	std::iota(order.begin(), order.end(), 0);
	Bytes payload;
	Put(payload, selected.size(), 4);
	for (auto c : selected) { Put(payload, c, 4); }
	Bytes cluster;
	if (!selected.empty()) {
		std::vector<Cell> projection;
		projection.reserve(size_t(page.rows) * selected.size());
		for (uint32_t r = 0; r < page.rows; r++) {
			for (auto c : selected) { projection.push_back(cells[size_t(r) * dim + c]); }
		}
		auto cluster_options = options;
		cluster_options.method -= 2;
		cluster_options.selected_columns.clear();
		Page cluster_page;
		cluster = EncodeCluster(projection, static_cast<uint32_t>(selected.size()), cluster_options,
		                        page_number, cluster_page, order);
		page.references = cluster_page.references;
		page.corrections = cluster_page.corrections;
	}
	Put(payload, cluster.size(), 8);
	payload.insert(payload.end(), cluster.begin(), cluster.end());
	for (auto c : remaining) {
		std::vector<double> column;
		column.reserve(page.rows);
		Bytes raw;
		bool supported = true;
		for (auto r : order) {
			const auto &cell = cells[size_t(r) * dim + c];
			Put(raw, cell.valid, 1);
			Put(raw, cell.valid ? cell.bits : 0, 8);
			const auto value = Value(cell.bits);
			supported = supported && cell.valid && std::isfinite(value) && std::abs(value) < 1e17;
			column.push_back(value);
		}
		Bytes encoded;
		if (supported) { encoded = Alp().encode(column); }
		const bool use_alp = supported && encoded.size() < raw.size();
		const auto &stream = use_alp ? encoded : raw;
		Put(payload, use_alp ? 1 : 0, 1);
		Put(payload, stream.size(), 8);
		payload.insert(payload.end(), stream.begin(), stream.end());
	}
	if (payload.size() >= cells.size() * 9) {
		Bytes raw;
		for (const auto &cell : cells) {
			Put(raw, cell.valid, 1);
			Put(raw, cell.valid ? cell.bits : 0, 8);
		}
		page.mode = 0;
		page.reason = 5;
		page.references = page.corrections = 0;
		return Frame(raw, page);
	}
	return Frame(payload, page);
}

std::vector<Cell> DecodeHybrid(const Bytes &bytes, uint32_t dim, const Page &page) {
	if (!dim || dim > 255 || !page.rows || page.rows > 10000 || page.mode < 3 || page.mode > 5 ||
	    page.reason || page.references > page.rows || page.corrections > uint64_t(page.rows) * dim ||
	    bytes.size() > uint64_t(page.rows) * dim * 9 + 36) {
		throw std::runtime_error("Invalid hybrid page shape");
	}
	Verify(bytes);
	Reader frame(bytes);
	if (frame.Get(4) != page.rows || frame.Get(4) != page.mode || frame.Get(4) != page.reason ||
	    frame.Get(4) != page.references || frame.Get(4) != page.corrections) {
		throw std::runtime_error("Hybrid page directory mismatch");
	}
	auto payload = frame.Block(frame.Get(8));
	frame.Get(8);
	frame.Finish();
	Reader reader(payload);
	const auto selected_count = reader.Get(4);
	if (selected_count > dim || (page.mode == 5 ? selected_count != 0 : selected_count == 0)) {
		throw std::runtime_error("Invalid hybrid selection count");
	}
	std::vector<uint32_t> selected;
	for (uint64_t c = 0; c < selected_count; c++) { selected.push_back(static_cast<uint32_t>(reader.Get(4))); }
	const auto remaining = RemainingColumns(dim, selected);
	auto cluster = reader.Block(reader.Get(8));
	std::vector<Cell> result(size_t(page.rows) * dim);
	if (!selected.empty()) {
		Reader probe(cluster);
		Page inner;
		inner.rows = static_cast<uint32_t>(probe.Get(4));
		inner.mode = static_cast<uint32_t>(probe.Get(4));
		inner.reason = static_cast<uint32_t>(probe.Get(4));
		inner.references = static_cast<uint32_t>(probe.Get(4));
		inner.corrections = static_cast<uint32_t>(probe.Get(4));
		if (inner.rows != page.rows || (inner.mode && inner.mode != page.mode - 2) ||
		    inner.references != page.references || inner.corrections != page.corrections) {
			throw std::runtime_error("Invalid hybrid cluster stream");
		}
		auto decoded = Decode(cluster, static_cast<uint32_t>(selected.size()), inner);
		for (uint32_t r = 0; r < page.rows; r++) {
			for (size_t c = 0; c < selected.size(); c++) { result[size_t(r) * dim + selected[c]] = decoded[size_t(r) * selected.size() + c]; }
		}
	} else if (!cluster.empty() || page.references || page.corrections) {
		throw std::runtime_error("Unexpected ALP cluster stream");
	}
	for (auto c : remaining) {
		const auto mode = reader.Get(1);
		auto stream = reader.Block(reader.Get(8));
		if (mode == 0) {
			if (stream.size() != size_t(page.rows) * 9) { throw std::runtime_error("Invalid raw column size"); }
			Reader raw(stream);
			for (uint32_t r = 0; r < page.rows; r++) {
				const auto valid = raw.Get(1);
				if (valid > 1) { throw std::runtime_error("Invalid raw column validity"); }
				result[size_t(r) * dim + c] = {raw.Get(8), valid == 1};
			}
		} else if (mode == 1) {
			Reader header(stream);
			if (header.Get(4) != page.rows || stream.size() > uint64_t(page.rows) * 9) {
				throw std::runtime_error("Invalid ALP stream size");
			}
			auto column = Alp().decode(stream);
			if (column.size() != page.rows) { throw std::runtime_error("Invalid ALP row count"); }
			for (uint32_t r = 0; r < page.rows; r++) { result[size_t(r) * dim + c].bits = Bits(column[r]); }
		} else { throw std::runtime_error("Invalid hybrid column mode"); }
	}
	reader.Finish();
	return result;
}
}
