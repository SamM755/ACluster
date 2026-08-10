#include "cluster_extension.hpp"
#include "cluster_codec.hpp"
#include "cluster_compression.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/vector_iterator.hpp"
#include "duckdb/common/vector/vector_writer.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/client_context.hpp"
#include <algorithm>

namespace duckdb {
namespace {
using cluster_codec::Bytes;
using cluster_codec::Cell;
using cluster_codec::Options;
using cluster_codec::Page;
using cluster_codec::Put;
using cluster_codec::Reader;
constexpr uint64_t HEADER_MAGIC = UINT64_C(0x3130305254534c43);
constexpr uint64_t FOOTER_MAGIC = UINT64_C(0x313030444e454c43);
constexpr uint64_t HEADER_MAGIC_V2 = UINT64_C(0x3230305254534c43);
constexpr uint64_t FOOTER_MAGIC_V2 = UINT64_C(0x323030444e454c43);
constexpr idx_t ENTRY_SIZE = 36;

struct FileInfo {
	vector<Identifier> names;
	Options options;
	std::vector<Page> pages;
	uint64_t total_rows = 0;
	uint64_t file_size = 0;
	uint64_t header_hash = 0;
	uint64_t footer_hash = 0;
};

Bytes ReadBytes(FileHandle &file, uint64_t offset, uint64_t length) {
	const auto size = file.GetFileSize();
	if (offset > size || length > size - offset) { throw IOException("Truncated cluster file"); }
	Bytes bytes(length);
	if (length) { file.Read(bytes.data(), length, offset); }
	return bytes;
}

FileInfo ReadInfoInternal(FileHandle &file) {
	FileInfo info;
	info.file_size = file.GetFileSize();
	if (info.file_size < 72) { throw IOException("Truncated cluster file"); }
	auto prefix = ReadBytes(file, 0, 16);
	Reader first(prefix);
	const auto magic = first.Get(8);
	const bool v2 = magic == HEADER_MAGIC_V2;
	if (magic != HEADER_MAGIC && !v2) { throw IOException("Not a cluster file"); }
	const auto entry_size = ENTRY_SIZE + (v2 ? 8 : 0);
	const auto header_length = first.Get(8);
	if (header_length > 2 * 1024 * 1024) { throw IOException("Invalid cluster header length"); }
	auto header = ReadBytes(file, 16, header_length);
	cluster_codec::Verify(header);
	Reader schema(header);
	const auto dim = schema.Get(4);
	if (!dim || dim > 255) { throw IOException("Invalid cluster column count"); }
	info.options.page_size = static_cast<uint32_t>(schema.Get(4));
	info.options.method = static_cast<uint32_t>(schema.Get(4));
	info.options.k = static_cast<uint32_t>(schema.Get(4));
	info.options.seed = static_cast<uint32_t>(schema.Get(4));
	if (v2) {
		info.options.compression = static_cast<uint32_t>(schema.Get(4));
		info.options.compression_level = static_cast<uint32_t>(schema.Get(4));
		if (info.options.compression != 1 || !info.options.compression_level || info.options.compression_level > 22) {
			throw IOException("Invalid cluster compression options");
		}
	}
	if (!info.options.page_size || info.options.page_size > 10000 || info.options.method > 5 ||
	    !info.options.k || info.options.k > 10000) { throw IOException("Invalid cluster file options"); }
	for (idx_t c = 0; c < dim; c++) { info.names.emplace_back(schema.String()); }
	info.header_hash = schema.Get(8);
	schema.Finish();
	auto trailer = ReadBytes(file, info.file_size - 24, 24);
	Reader tail(trailer);
	const auto footer_offset = tail.Get(8);
	const auto footer_length = tail.Get(8);
	if (tail.Get(8) != (v2 ? FOOTER_MAGIC_V2 : FOOTER_MAGIC) || footer_offset < 16 + header_length ||
	    footer_offset > info.file_size - 24 || footer_length != info.file_size - 24 - footer_offset ||
	    footer_length < 24 || (footer_length - 24) % entry_size) {
		throw IOException("Invalid or incomplete cluster footer");
	}
	auto footer = ReadBytes(file, footer_offset, footer_length);
	cluster_codec::Verify(footer);
	Reader directory(footer);
	info.total_rows = directory.Get(8);
	const auto page_count = directory.Get(8);
	if (page_count != (footer_length - 24) / entry_size) { throw IOException("Invalid cluster page count"); }
	uint64_t next_offset = 16 + header_length;
	uint64_t total = 0;
	for (uint64_t p = 0; p < page_count; p++) {
		Page page;
		page.offset = directory.Get(8);
		page.length = directory.Get(8);
		page.rows = static_cast<uint32_t>(directory.Get(4));
		page.mode = static_cast<uint32_t>(directory.Get(4));
		page.reason = static_cast<uint32_t>(directory.Get(4));
		page.references = static_cast<uint32_t>(directory.Get(4));
		page.corrections = static_cast<uint32_t>(directory.Get(4));
		page.uncompressed_length = v2 ? directory.Get(8) : page.length;
		const auto maximum = uint64_t(page.rows) * dim * 9 + 36;
		if (page.offset != next_offset || page.offset > footer_offset || page.length > footer_offset - page.offset ||
		    page.length < (v2 ? 17 : 36) || !page.rows || page.rows > info.options.page_size ||
		    page.uncompressed_length < 36 || page.uncompressed_length > maximum ||
		    page.length > (v2 ? cluster_codec::CompressedPageBound(maximum) : maximum) || page.mode > 5 || page.reason > 5 ||
		    (page.mode && page.mode != info.options.method) || page.references > page.rows ||
		    page.corrections > uint64_t(page.rows) * dim) {
			throw IOException("Invalid cluster page directory entry");
		}
		next_offset += page.length;
		total += page.rows;
		info.pages.push_back(page);
	}
	info.footer_hash = directory.Get(8);
	directory.Finish();
	if (next_offset != footer_offset || total != info.total_rows) { throw IOException("Cluster file totals mismatch"); }
	return info;
}

FileInfo ReadInfo(FileHandle &file) {
	try {
		return ReadInfoInternal(file);
	} catch (const Exception &) {
		throw;
	} catch (const std::runtime_error &error) {
		throw IOException("Invalid cluster file: %s", error.what());
	}
}

struct WriteBind : public TableFunctionData {
	vector<Identifier> names;
	Options options;
};
struct WriteLocal : public LocalFunctionData {};
struct WriteGlobal : public GlobalFunctionData {
	unique_ptr<FileHandle> file;
	mutex lock;
	std::vector<Cell> cells;
	std::vector<Page> pages;
	uint64_t position = 0;
	uint64_t total_rows = 0;
	bool finished = false;
	~WriteGlobal() override {
		if (file && !finished) {
			try { file->AbortWrite(); } catch (...) { }
		}
	}
	void Write(Bytes &bytes) {
		if (file->Write(bytes.data(), bytes.size()) != int64_t(bytes.size())) { throw IOException("Short cluster write"); }
		position += bytes.size();
	}
	void Flush(const WriteBind &bind) {
		if (cells.empty()) { return; }
		Page page;
		auto bytes = cluster_codec::Encode(cells, static_cast<uint32_t>(bind.names.size()), bind.options, pages.size(), page);
		page.uncompressed_length = bytes.size();
		if (bind.options.compression) { bytes = cluster_codec::CompressPage(bytes, bind.options.compression_level); }
		page.offset = position;
		page.length = bytes.size();
		Write(bytes);
		total_rows += page.rows;
		pages.push_back(page);
		cells.clear();
	}
};

unique_ptr<FunctionData> WriteBindFunction(ClientContext &, CopyFunctionBindInput &input,
                                          const vector<Identifier> &names, const vector<LogicalType> &types) {
	if (names.empty() || names.size() > 255) { throw BinderException("cluster supports 1 to 255 DOUBLE columns"); }
	auto result = make_uniq<WriteBind>();
	result->names = names;
	for (idx_t i = 0; i < types.size(); i++) {
		if (types[i] != LogicalType::DOUBLE) { throw BinderException("cluster v1 requires DOUBLE columns; cast explicitly"); }
		if (names[i].size() > 4096) { throw BinderException("cluster column name too long"); }
	}
	for (const auto &option : input.info.options) {
		if (option.second.size() != 1 || option.second[0].IsNull()) { throw BinderException("Invalid cluster option"); }
		const auto key = StringUtil::Lower(option.first.GetIdentifierName());
		const auto &value = option.second[0];
		if (key == "method") {
			const auto method = StringUtil::Lower(value.ToString());
			if (method == "raw") { result->options.method = 0; }
			else if (method == "acluster") { result->options.method = 1; }
			else if (method == "kcluster") { result->options.method = 2; }
			else if (method == "acluster_alp") { result->options.method = 3; }
			else if (method == "kcluster_alp") { result->options.method = 4; }
			else if (method == "alp") { result->options.method = 5; }
			else { throw BinderException("Unknown cluster METHOD"); }
		} else if (key == "cluster_columns") {
			for (const auto &token : StringUtil::Split(value.ToString(), ',')) {
				if (token.empty() || token.size() > 3 || token.find_first_not_of("0123456789") != string::npos) {
					throw BinderException("CLUSTER_COLUMNS requires comma-separated zero-based column indices");
				}
				const auto c = static_cast<uint32_t>(std::stoul(token));
				if (c >= names.size() || std::find(result->options.selected_columns.begin(), result->options.selected_columns.end(), c) != result->options.selected_columns.end()) {
					throw BinderException("Invalid or duplicate CLUSTER_COLUMNS index");
				}
				result->options.selected_columns.push_back(c);
			}
		} else if (key == "compression") {
			const auto compression = StringUtil::Lower(value.ToString());
			if (compression == "none" || compression == "uncompressed") { result->options.compression = 0; }
			else if (compression == "zstd") { result->options.compression = 1; }
			else { throw BinderException("cluster COMPRESSION must be none or zstd"); }
		} else if (key == "compression_level") {
			const auto level = value.GetValue<int64_t>();
			if (level < 1 || level > 22) { throw BinderException("cluster COMPRESSION_LEVEL must be between 1 and 22"); }
			result->options.compression_level = static_cast<uint32_t>(level);
		} else if (key == "page_size" || key == "k" || key == "seed") {
			const auto number = value.GetValue<int64_t>();
			if (number < 0 || (key != "seed" && (number < 1 || number > 10000)) || uint64_t(number) > UINT32_MAX) {
				throw BinderException("Invalid cluster %s", key);
			}
			if (key == "page_size") { result->options.page_size = static_cast<uint32_t>(number); }
			if (key == "k") { result->options.k = static_cast<uint32_t>(number); }
			if (key == "seed") { result->options.seed = static_cast<uint32_t>(number); }
		} else { throw BinderException("Unknown cluster option: %s", key); }
	}
	const auto hybrid = result->options.method == 3 || result->options.method == 4;
	if (hybrid == result->options.selected_columns.empty()) {
		throw BinderException("CLUSTER_COLUMNS is required only for acluster_alp and kcluster_alp");
	}
	return std::move(result);
}

void CopyOptions(ClientContext &, CopyOptionsInput &input) {
	input.options["method"] = CopyOption(LogicalType::VARCHAR, CopyOptionMode::WRITE_ONLY);
	input.options["cluster_columns"] = CopyOption(LogicalType::VARCHAR, CopyOptionMode::WRITE_ONLY);
	input.options["page_size"] = CopyOption(LogicalType::BIGINT, CopyOptionMode::WRITE_ONLY);
	input.options["k"] = CopyOption(LogicalType::BIGINT, CopyOptionMode::WRITE_ONLY);
	input.options["seed"] = CopyOption(LogicalType::BIGINT, CopyOptionMode::WRITE_ONLY);
	input.options["compression"] = CopyOption(LogicalType::VARCHAR, CopyOptionMode::WRITE_ONLY);
	input.options["compression_level"] = CopyOption(LogicalType::BIGINT, CopyOptionMode::WRITE_ONLY);
}

unique_ptr<GlobalFunctionData> WriteInit(ClientContext &context, FunctionData &data, const string &path) {
	const auto &bind = data.Cast<WriteBind>();
	auto state = make_uniq<WriteGlobal>();
	state->file = FileSystem::GetFileSystem(context).OpenFile(path,
	    FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE | FileFlags::FILE_FLAGS_EXCLUSIVE_CREATE);
	Bytes header;
	Put(header, bind.names.size(), 4);
	Put(header, bind.options.page_size, 4);
	Put(header, bind.options.method, 4);
	Put(header, bind.options.k, 4);
	Put(header, bind.options.seed, 4);
	if (bind.options.compression) {
		Put(header, bind.options.compression, 4);
		Put(header, bind.options.compression_level, 4);
	}
	for (const auto &name : bind.names) {
		Put(header, name.size(), 4);
		const auto &text = name.GetIdentifierName();
		header.insert(header.end(), text.begin(), text.end());
	}
	cluster_codec::Seal(header);
	Bytes prefix;
	Put(prefix, bind.options.compression ? HEADER_MAGIC_V2 : HEADER_MAGIC, 8);
	Put(prefix, header.size(), 8);
	state->Write(prefix);
	state->Write(header);
	return std::move(state);
}

void WriteSink(ExecutionContext &, FunctionData &data, GlobalFunctionData &global, LocalFunctionData &, DataChunk &input) {
	const auto &bind = data.Cast<WriteBind>();
	auto &state = global.Cast<WriteGlobal>();
	lock_guard<mutex> guard(state.lock);
	const auto dim = bind.names.size();
	idx_t consumed = 0;
	while (consumed < input.size()) {
		const auto start = state.cells.size();
		const auto count = MinValue<idx_t>(input.size() - consumed, bind.options.page_size - start / dim);
		state.cells.resize(start + count * dim);
		for (idx_t c = 0; c < dim; c++) {
			auto values = input.data[c].Values<double>();
			for (idx_t r = 0; r < count; r++) {
				auto entry = values[consumed + r];
				auto &cell = state.cells[start + r * dim + c];
				cell.valid = entry.IsValid();
				cell.bits = cell.valid ? cluster_codec::Bits(entry.GetValue()) : 0;
			}
		}
		consumed += count;
		if (state.cells.size() / dim == bind.options.page_size) { state.Flush(bind); }
	}
}

void WriteFinalize(ClientContext &, FunctionData &data, GlobalFunctionData &global) {
	auto &state = global.Cast<WriteGlobal>();
	const auto &bind = data.Cast<WriteBind>();
	state.Flush(bind);
	const auto footer_offset = state.position;
	Bytes footer;
	Put(footer, state.total_rows, 8);
	Put(footer, state.pages.size(), 8);
	for (const auto &page : state.pages) {
		Put(footer, page.offset, 8);
		Put(footer, page.length, 8);
		Put(footer, page.rows, 4);
		Put(footer, page.mode, 4);
		Put(footer, page.reason, 4);
		Put(footer, page.references, 4);
		Put(footer, page.corrections, 4);
		if (bind.options.compression) { Put(footer, page.uncompressed_length, 8); }
	}
	cluster_codec::Seal(footer);
	state.Write(footer);
	Bytes trailer;
	Put(trailer, footer_offset, 8);
	Put(trailer, footer.size(), 8);
	Put(trailer, bind.options.compression ? FOOTER_MAGIC_V2 : FOOTER_MAGIC, 8);
	state.Write(trailer);
	state.file->Sync();
	state.file->Close();
	state.finished = true;
}

struct ReadBind : public TableFunctionData {
	string path;
	FileInfo info;
};
struct ReadGlobal : public GlobalTableFunctionState {
	unique_ptr<FileHandle> file;
	std::vector<Cell> cells;
	idx_t next_page = 0;
	idx_t next_row = 0;
	idx_t MaxThreads() const override { return 1; }
};

unique_ptr<FunctionData> BindReader(ClientContext &context, TableFunctionBindInput &input,
                                    vector<LogicalType> &types, vector<Identifier> &names, bool metadata) {
	if (input.inputs[0].IsNull()) { throw BinderException("cluster path cannot be NULL"); }
	auto bind = make_uniq<ReadBind>();
	bind->path = input.inputs[0].GetValue<string>();
	auto file = FileSystem::GetFileSystem(context).OpenFile(bind->path, FileFlags::FILE_FLAGS_READ);
	bind->info = ReadInfo(*file);
	if (metadata) {
		names = {"page_id", "rows", "mode", "fallback_reason", "references", "corrections", "page_bytes", "file_bytes",
		         "compression", "codec_page_bytes", "compression_level"};
		types = {LogicalType::UBIGINT, LogicalType::UBIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR,
		         LogicalType::UBIGINT, LogicalType::UBIGINT, LogicalType::UBIGINT, LogicalType::UBIGINT,
		         LogicalType::VARCHAR, LogicalType::UBIGINT, LogicalType::UBIGINT};
	} else {
		names = bind->info.names;
		types.assign(names.size(), LogicalType::DOUBLE);
	}
	return std::move(bind);
}

unique_ptr<GlobalTableFunctionState> ReadInit(ClientContext &context, TableFunctionInitInput &input) {
	const auto &bind = input.bind_data->Cast<ReadBind>();
	auto state = make_uniq<ReadGlobal>();
	state->file = FileSystem::GetFileSystem(context).OpenFile(bind.path, FileFlags::FILE_FLAGS_READ);
	auto current = ReadInfo(*state->file);
	if (current.file_size != bind.info.file_size || current.header_hash != bind.info.header_hash ||
	    current.footer_hash != bind.info.footer_hash) { throw IOException("Cluster file changed after binding"); }
	return std::move(state);
}

void Scan(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	const auto &bind = input.bind_data->Cast<ReadBind>();
	auto &state = input.global_state->Cast<ReadGlobal>();
	const auto dim = bind.info.names.size();
	if (state.next_row == state.cells.size() / dim) {
		if (state.next_page == bind.info.pages.size()) { return; }
		context.InterruptCheck();
		const auto &page = bind.info.pages[state.next_page++];
		try {
			auto bytes = ReadBytes(*state.file, page.offset, page.length);
			if (bind.info.options.compression) {
				bytes = cluster_codec::DecompressPage(bytes, page.uncompressed_length, uint64_t(page.rows) * dim * 9 + 36);
			}
			state.cells = cluster_codec::Decode(bytes, static_cast<uint32_t>(dim), page);
		} catch (const Exception &) {
			throw;
		} catch (const std::runtime_error &error) {
			throw IOException("Invalid cluster page: %s", error.what());
		}
		state.next_row = 0;
	}
	const auto count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, state.cells.size() / dim - state.next_row);
	for (idx_t c = 0; c < dim; c++) {
		auto writer = FlatVector::Writer<double>(output.data[c], count);
		for (idx_t r = 0; r < count; r++) {
			const auto &cell = state.cells[(state.next_row + r) * dim + c];
			if (cell.valid) { writer.WriteValue(cluster_codec::Value(cell.bits)); }
			else { writer.WriteNull(); }
		}
	}
	output.CheckCardinality(count);
	state.next_row += count;
}

void MetadataScan(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	const auto &info = input.bind_data->Cast<ReadBind>().info;
	auto &state = input.global_state->Cast<ReadGlobal>();
	const auto count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, info.pages.size() - state.next_page);
	for (idx_t c = 0; c < 11; c++) {
		if (c == 2 || c == 3 || c == 8) {
			auto writer = FlatVector::Writer<string_t>(output.data[c], count);
			for (idx_t r = 0; r < count; r++) {
				const auto &page = info.pages[state.next_page + r];
				const auto text = c == 2 ? cluster_codec::MethodName(page.mode) : c == 3 ? cluster_codec::ReasonName(page.reason) :
				                  info.options.compression ? "zstd" : "none";
				writer.WriteValue(string_t(text));
			}
		} else {
			auto writer = FlatVector::Writer<uint64_t>(output.data[c], count);
			for (idx_t r = 0; r < count; r++) {
				const auto &page = info.pages[state.next_page + r];
				const uint64_t value = c == 0 ? state.next_page + r : c == 1 ? page.rows : c == 4 ? page.references :
				                       c == 5 ? page.corrections : c == 6 ? page.length : c == 7 ? info.file_size :
				                       c == 9 ? page.uncompressed_length : info.options.compression ? info.options.compression_level : 0;
				writer.WriteValue(value);
			}
		}
	}
	output.CheckCardinality(count);
	state.next_page += count;
}

void BitsFunction(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<double, string_t>(args.data[0], result, args.size(), [&result](double value) {
		const char digits[] = "0123456789abcdef";
		auto bits = cluster_codec::Bits(value);
		string text(16, '0');
		for (idx_t i = 0; i < 16; i++) { text[15 - i] = digits[(bits >> (i * 4)) & 15]; }
		return StringVector::AddString(result, text);
	});
}

void LoadInternal(ExtensionLoader &loader) {
	CopyFunction copy("cluster");
	copy.extension = "cluster";
	copy.copy_options = CopyOptions;
	copy.copy_to_bind = WriteBindFunction;
	copy.copy_to_initialize_global = WriteInit;
	copy.copy_to_initialize_local = [](ExecutionContext &, FunctionData &) -> unique_ptr<LocalFunctionData> {
		return make_uniq<WriteLocal>();
	};
	copy.copy_to_sink = WriteSink;

	copy.copy_to_combine = [](ExecutionContext &, FunctionData &, GlobalFunctionData &, LocalFunctionData &) {};
	copy.copy_to_finalize = WriteFinalize;
	copy.execution_mode = [](bool, bool) { return CopyFunctionExecutionMode::REGULAR_COPY_TO_FILE; };
	loader.RegisterFunction(copy);
	TableFunction scan("read_cluster", {LogicalType::VARCHAR}, Scan,
	    [](ClientContext &c, TableFunctionBindInput &i, vector<LogicalType> &t, vector<Identifier> &n) {
		    return BindReader(c, i, t, n, false);
	    }, ReadInit);
	scan.projection_pushdown = false;
	scan.filter_pushdown = false;
	scan.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(scan);
	TableFunction metadata("cluster_metadata", {LogicalType::VARCHAR}, MetadataScan,
	    [](ClientContext &c, TableFunctionBindInput &i, vector<LogicalType> &t, vector<Identifier> &n) {
		    return BindReader(c, i, t, n, true);
	    }, ReadInit);
	metadata.projection_pushdown = false;
	metadata.filter_pushdown = false;
	loader.RegisterFunction(metadata);
	ScalarFunction bits("cluster_bits", {LogicalType::DOUBLE}, LogicalType::VARCHAR, BitsFunction);
	loader.RegisterFunction(bits);
}
}
void ClusterExtension::Load(ExtensionLoader &loader) { LoadInternal(loader); }
}
extern "C" {
DUCKDB_CPP_EXTENSION_ENTRY(cluster, loader) { duckdb::LoadInternal(loader); }
}
