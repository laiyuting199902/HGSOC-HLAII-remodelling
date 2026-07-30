#include <zlib.h>

#include <cerrno>
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
  while (*cursor == ' ' || *cursor == '\t') {
    ++cursor;
  }
  if (*cursor < '0' || *cursor > '9') {
    fail("expected an unsigned integer in Matrix Market input");
  }
  uint64_t value = 0;
  while (*cursor >= '0' && *cursor <= '9') {
    const uint64_t digit = static_cast<uint64_t>(*cursor - '0');
    if (value > (std::numeric_limits<uint64_t>::max() - digit) / 10) {
      fail("integer overflow while parsing Matrix Market input");
    }
    value = value * 10 + digit;
    ++cursor;
  }
  return value;
}

std::vector<std::pair<uint32_t, uint32_t>> read_mapping(const std::string& path) {
  std::ifstream input(path);
  if (!input) {
    fail("could not open mapping file: " + path);
  }
  std::string line;
  if (!std::getline(input, line) || line != "cell_index\tgroup_index") {
    fail("mapping file must start with cell_index and group_index columns");
  }
  std::vector<std::pair<uint32_t, uint32_t>> mapping;
  uint64_t cell = 0;
  uint64_t group = 0;
  while (input >> cell >> group) {
    if (cell == 0 || cell > std::numeric_limits<uint32_t>::max() ||
        group == 0 || group > std::numeric_limits<uint32_t>::max()) {
      fail("mapping indices must be positive 32-bit integers");
    }
    mapping.emplace_back(static_cast<uint32_t>(cell), static_cast<uint32_t>(group));
  }
  if (!input.eof()) {
    fail("invalid row in mapping file");
  }
  if (mapping.empty()) {
    fail("mapping file does not contain selected cells");
  }
  return mapping;
}

}  // namespace

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "Usage: stream_mtx_pseudobulk MATRIX.mtx.gz CELL_MAP.tsv OUTPUT.tsv QC.tsv GROUP_COUNT\n";
    return EXIT_FAILURE;
  }

  const std::string matrix_path = argv[1];
  const std::string mapping_path = argv[2];
  const std::string output_path = argv[3];
  const std::string qc_path = argv[4];
  const uint64_t parsed_group_count = std::strtoull(argv[5], nullptr, 10);
  if (parsed_group_count == 0 || parsed_group_count > std::numeric_limits<uint32_t>::max()) {
    fail("GROUP_COUNT must be a positive 32-bit integer");
  }
  const uint32_t group_count = static_cast<uint32_t>(parsed_group_count);

  gzFile input = gzopen(matrix_path.c_str(), "rb");
  if (input == nullptr) {
    fail("could not open gzip Matrix Market input: " + matrix_path);
  }
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
  if (row_count == 0 || column_count == 0 ||
      row_count > std::numeric_limits<size_t>::max() / group_count) {
    gzclose(input);
    fail("invalid or unsupported Matrix Market dimensions");
  }

  const auto mapping = read_mapping(mapping_path);
  std::vector<int32_t> cell_to_group(static_cast<size_t>(column_count) + 1, -1);
  std::vector<uint64_t> selected_cells(group_count, 0);
  for (const auto& item : mapping) {
    const uint32_t cell = item.first;
    const uint32_t one_based_group = item.second;
    if (cell > column_count || one_based_group > group_count) {
      gzclose(input);
      fail("mapping index exceeds matrix columns or GROUP_COUNT");
    }
    if (cell_to_group[cell] != -1) {
      gzclose(input);
      fail("mapping contains a duplicate cell index");
    }
    const int32_t zero_based_group = static_cast<int32_t>(one_based_group - 1);
    cell_to_group[cell] = zero_based_group;
    ++selected_cells[zero_based_group];
  }

  std::vector<uint64_t> counts(static_cast<size_t>(row_count) * group_count, 0);
  std::vector<uint64_t> total_counts(group_count, 0);
  uint64_t observed_nonzero = 0;
  while (gzgets(input, line, sizeof(line)) != nullptr) {
    if (line[0] == '%') {
      continue;
    }
    const char* cursor = line;
    const uint64_t row = parse_uint(cursor);
    const uint64_t column = parse_uint(cursor);
    const uint64_t value = parse_uint(cursor);
    if (row == 0 || row > row_count || column == 0 || column > column_count) {
      gzclose(input);
      fail("Matrix Market coordinate is outside declared dimensions");
    }
    const int32_t group = cell_to_group[static_cast<size_t>(column)];
    if (group >= 0) {
      const size_t offset = static_cast<size_t>(row - 1) * group_count + static_cast<size_t>(group);
      if (counts[offset] > std::numeric_limits<uint64_t>::max() - value ||
          total_counts[group] > std::numeric_limits<uint64_t>::max() - value) {
        gzclose(input);
        fail("pseudobulk count overflow");
      }
      counts[offset] += value;
      total_counts[group] += value;
    }
    ++observed_nonzero;
  }
  if (gzclose(input) != Z_OK) {
    fail("gzip stream ended with an error");
  }
  if (observed_nonzero != expected_nonzero) {
    fail("observed Matrix Market entries do not match the declared nonzero count");
  }

  std::ofstream output(output_path);
  if (!output) {
    fail("could not create pseudobulk output: " + output_path);
  }
  output << "feature_index";
  for (uint32_t group = 0; group < group_count; ++group) {
    output << "\tgroup_" << group + 1;
  }
  output << '\n';
  for (uint64_t row = 0; row < row_count; ++row) {
    output << row + 1;
    const size_t offset = static_cast<size_t>(row) * group_count;
    for (uint32_t group = 0; group < group_count; ++group) {
      output << '\t' << counts[offset + group];
    }
    output << '\n';
  }
  if (!output) {
    fail("failed while writing pseudobulk output");
  }

  std::ofstream qc(qc_path);
  if (!qc) {
    fail("could not create pseudobulk QC output: " + qc_path);
  }
  qc << "group_index\tselected_cells\ttotal_counts\n";
  for (uint32_t group = 0; group < group_count; ++group) {
    qc << group + 1 << '\t' << selected_cells[group] << '\t' << total_counts[group] << '\n';
  }
  if (!qc) {
    fail("failed while writing pseudobulk QC output");
  }

  return EXIT_SUCCESS;
}
