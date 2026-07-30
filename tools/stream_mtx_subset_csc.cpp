#include <zlib.h>

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

[[noreturn]] void fail(const std::string& message) {
  std::cerr << "ERROR: " << message << '\n';
  std::exit(EXIT_FAILURE);
}

uint64_t parse_uint(const char*& cursor) {
  while (*cursor == ' ' || *cursor == '\t') ++cursor;
  if (*cursor < '0' || *cursor > '9') fail("expected unsigned integer in Matrix Market input");
  uint64_t value = 0;
  while (*cursor >= '0' && *cursor <= '9') {
    const uint64_t digit = static_cast<uint64_t>(*cursor - '0');
    if (value > (std::numeric_limits<uint64_t>::max() - digit) / 10) fail("integer parse overflow");
    value = value * 10 + digit;
    ++cursor;
  }
  return value;
}

std::vector<std::pair<uint32_t, uint32_t>> read_mapping(const std::string& path) {
  std::ifstream input(path);
  if (!input) fail("could not open subset mapping: " + path);
  std::string header;
  if (!std::getline(input, header) || header != "cell_index\tsubset_index") {
    fail("subset mapping must start with cell_index and subset_index columns");
  }
  std::vector<std::pair<uint32_t, uint32_t>> mapping;
  uint64_t cell = 0;
  uint64_t subset = 0;
  while (input >> cell >> subset) {
    if (cell == 0 || subset == 0 || cell > std::numeric_limits<uint32_t>::max() ||
        subset > std::numeric_limits<uint32_t>::max()) {
      fail("subset mapping index is outside positive 32-bit range");
    }
    mapping.emplace_back(static_cast<uint32_t>(cell), static_cast<uint32_t>(subset));
  }
  if (!input.eof() || mapping.empty()) fail("subset mapping is empty or malformed");
  for (size_t index = 0; index < mapping.size(); ++index) {
    if (mapping[index].second != index + 1 || (index > 0 && mapping[index].first <= mapping[index - 1].first)) {
      fail("subset mapping must have increasing cell_index and consecutive subset_index");
    }
  }
  return mapping;
}

void flush_buffers(
    std::ofstream& i_output,
    std::ofstream& x_output,
    std::vector<uint32_t>& i_buffer,
    std::vector<uint32_t>& x_buffer) {
  if (i_buffer.empty()) return;
  i_output.write(reinterpret_cast<const char*>(i_buffer.data()), static_cast<std::streamsize>(i_buffer.size() * sizeof(uint32_t)));
  x_output.write(reinterpret_cast<const char*>(x_buffer.data()), static_cast<std::streamsize>(x_buffer.size() * sizeof(uint32_t)));
  if (!i_output || !x_output) fail("failed while writing CSC binary buffers");
  i_buffer.clear();
  x_buffer.clear();
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "Usage: stream_mtx_subset_csc MATRIX.mtx.gz SUBSET_MAP.tsv OUTPUT_PREFIX\n";
    return EXIT_FAILURE;
  }
  const std::string matrix_path = argv[1];
  const std::string mapping_path = argv[2];
  const std::string prefix = argv[3];
  const auto mapping = read_mapping(mapping_path);

  gzFile input = gzopen(matrix_path.c_str(), "rb");
  if (input == nullptr) fail("could not open gzip Matrix Market input: " + matrix_path);
  gzbuffer(input, 1 << 20);
  char line[256];
  if (gzgets(input, line, sizeof(line)) == nullptr ||
      std::string(line).find("%%MatrixMarket matrix coordinate integer") != 0) {
    gzclose(input);
    fail("input is not an integer coordinate Matrix Market file");
  }
  do {
    if (gzgets(input, line, sizeof(line)) == nullptr) {
      gzclose(input);
      fail("Matrix Market dimensions line was not found");
    }
  } while (line[0] == '%');
  const char* dimensions = line;
  const uint64_t row_count = parse_uint(dimensions);
  const uint64_t column_count = parse_uint(dimensions);
  const uint64_t expected_nonzero = parse_uint(dimensions);
  if (row_count == 0 || row_count > std::numeric_limits<uint32_t>::max() || column_count == 0) {
    gzclose(input);
    fail("unsupported Matrix Market dimensions");
  }

  std::vector<int32_t> cell_to_subset(static_cast<size_t>(column_count) + 1, -1);
  for (const auto& item : mapping) {
    if (item.first > column_count) {
      gzclose(input);
      fail("subset mapping cell exceeds matrix column count");
    }
    cell_to_subset[item.first] = static_cast<int32_t>(item.second - 1);
  }

  std::ofstream i_output(prefix + "_i.bin", std::ios::binary);
  std::ofstream x_output(prefix + "_x.bin", std::ios::binary);
  if (!i_output || !x_output) {
    gzclose(input);
    fail("could not create CSC binary outputs");
  }
  std::vector<uint64_t> nnz_by_subset(mapping.size(), 0);
  std::vector<uint32_t> i_buffer;
  std::vector<uint32_t> x_buffer;
  i_buffer.reserve(1 << 20);
  x_buffer.reserve(1 << 20);
  uint64_t observed_nonzero = 0;
  uint64_t selected_nonzero = 0;
  uint64_t previous_column = 0;
  uint64_t previous_row = 0;

  while (gzgets(input, line, sizeof(line)) != nullptr) {
    if (line[0] == '%') continue;
    const char* cursor = line;
    const uint64_t row = parse_uint(cursor);
    const uint64_t column = parse_uint(cursor);
    const uint64_t value = parse_uint(cursor);
    if (row == 0 || row > row_count || column == 0 || column > column_count ||
        value > std::numeric_limits<uint32_t>::max()) {
      gzclose(input);
      fail("Matrix Market coordinate or value is outside supported range");
    }
    if (column < previous_column || (column == previous_column && row <= previous_row)) {
      gzclose(input);
      fail("Matrix Market entries must be sorted by column and strictly by row within column");
    }
    previous_row = column == previous_column ? row : row;
    previous_column = column;
    const int32_t subset = cell_to_subset[static_cast<size_t>(column)];
    if (subset >= 0) {
      i_buffer.push_back(static_cast<uint32_t>(row - 1));
      x_buffer.push_back(static_cast<uint32_t>(value));
      ++nnz_by_subset[static_cast<size_t>(subset)];
      ++selected_nonzero;
      if (i_buffer.size() == i_buffer.capacity()) flush_buffers(i_output, x_output, i_buffer, x_buffer);
    }
    ++observed_nonzero;
  }
  flush_buffers(i_output, x_output, i_buffer, x_buffer);
  if (gzclose(input) != Z_OK) fail("gzip stream ended with an error");
  if (observed_nonzero != expected_nonzero) fail("observed entries do not match declared nonzero count");
  if (selected_nonzero > static_cast<uint64_t>(std::numeric_limits<int32_t>::max())) {
    fail("selected matrix exceeds the 32-bit Matrix CSC limit");
  }

  std::ofstream p_output(prefix + "_p.tsv");
  p_output << "column_pointer\n0\n";
  uint64_t cumulative = 0;
  for (const uint64_t count : nnz_by_subset) {
    cumulative += count;
    p_output << cumulative << '\n';
  }
  if (!p_output) fail("failed while writing CSC column pointers");

  std::ofstream manifest(prefix + "_manifest.tsv");
  manifest << "key\tvalue\n"
           << "features\t" << row_count << '\n'
           << "cells\t" << mapping.size() << '\n'
           << "nonzero\t" << selected_nonzero << '\n'
           << "endian\tlittle\n"
           << "index_base\tzero\n"
           << "value_type\tuint32\n";
  if (!manifest) fail("failed while writing CSC manifest");
  return EXIT_SUCCESS;
}
