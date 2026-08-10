import ctypes as ct
import os
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
names = ['libduckdb.dll', 'libduckdb.so', 'libduckdb.dylib']
DLL = next((ROOT / '.build/src' / name for name in names if (ROOT / '.build/src' / name).is_file()), None)
if DLL is None:
    raise RuntimeError('Build DuckDB with scripts/build.sh first')
_dll_directories = []
if os.name == 'nt':
    for folder in os.environ.get('PATH', '').split(os.pathsep):
        try:
            if folder and Path(folder).is_dir():
                _dll_directories.append(os.add_dll_directory(str(Path(folder).resolve())))
        except OSError:
            continue
LIB = ct.CDLL(str(DLL))

class Result(ct.Structure):

    _fields_ = [("column_count", ct.c_uint64), ("row_count", ct.c_uint64),
                ("rows_changed", ct.c_uint64), ("columns", ct.c_void_p),
                ("error", ct.c_void_p), ("internal", ct.c_void_p)]

def signature(name, args, result):
    function = getattr(LIB, name)
    function.argtypes, function.restype = args, result
    return function

open_db = signature("duckdb_open", [ct.c_char_p, ct.POINTER(ct.c_void_p)], ct.c_int)
close_db = signature("duckdb_close", [ct.POINTER(ct.c_void_p)], None)
connect = signature("duckdb_connect", [ct.c_void_p, ct.POINTER(ct.c_void_p)], ct.c_int)
disconnect = signature("duckdb_disconnect", [ct.POINTER(ct.c_void_p)], None)
query = signature("duckdb_query", [ct.c_void_p, ct.c_char_p, ct.POINTER(Result)], ct.c_int)
destroy = signature("duckdb_destroy_result", [ct.POINTER(Result)], None)
error = signature("duckdb_result_error", [ct.POINTER(Result)], ct.c_char_p)
rows = signature("duckdb_row_count", [ct.POINTER(Result)], ct.c_uint64)
columns = signature("duckdb_column_count", [ct.POINTER(Result)], ct.c_uint64)
column_name = signature("duckdb_column_name", [ct.POINTER(Result), ct.c_uint64], ct.c_char_p)
value = signature("duckdb_value_varchar", [ct.POINTER(Result), ct.c_uint64, ct.c_uint64], ct.c_void_p)
free = signature("duckdb_free", [ct.c_void_p], None)
version = signature("duckdb_library_version", [], ct.c_char_p)

def literal(value):
    return "'" + str(value).replace("\\", "/").replace("'", "''") + "'"

class Connection:
    def __init__(self):
        self.db, self.connection = ct.c_void_p(), ct.c_void_p()
        if open_db(None, ct.byref(self.db)) or connect(self.db, ct.byref(self.connection)):
            raise RuntimeError("Cannot open locally built DuckDB")
        self.execute("LOAD cluster; SET threads=1; SET memory_limit='16GB'; SET preserve_insertion_order=true;")

    def execute(self, sql):
        result = Result()
        started = time.perf_counter_ns()
        state = query(self.connection, sql.encode("utf-8"), ct.byref(result))
        elapsed = time.perf_counter_ns() - started
        try:
            if state:
                raise RuntimeError(error(ct.byref(result)).decode("utf-8") + "\nSQL: " + sql)
            names = [column_name(ct.byref(result), c).decode("utf-8") for c in range(columns(ct.byref(result)))]
            output = []
            for r in range(rows(ct.byref(result))):
                record = {}
                for c, name in enumerate(names):
                    address = value(ct.byref(result), c, r)
                    try:
                        record[name] = ct.string_at(address).decode("utf-8") if address else None
                    finally:
                        if address:
                            free(address)
                output.append(record)
            return output, elapsed
        finally:
            destroy(ct.byref(result))

    def close(self):
        if self.connection:
            disconnect(ct.byref(self.connection))
        if self.db:
            close_db(ct.byref(self.db))

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

if __name__ == "__main__":
    with Connection() as connection:
        print(version().decode(), connection.execute("SELECT cluster_bits(1::DOUBLE) AS bits")[0])
