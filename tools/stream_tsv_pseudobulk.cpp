#include <zlib.h>

#include <cerrno>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool gz_getline(gzFile handle, std::string& line) {
  line.clear();
  constexpr int buffer_size = 1 << 16;
  char buffer[buffer_size];
  while (true) {
    char* result = gzgets(handle, buffer, buffer_size);
    if (result == nullptr) return !line.empty();
    line.append(buffer);
    if (!line.empty() && line.back() == '\n') {
      line.pop_back();
      if (!line.empty() && line.back() == '\r') line.pop_back();
      return true;
    }
    if (gzeof(handle)) return !line.empty();
  }
}

std::size_t field_count(const std::string& line) {
  std::size_t count = 1;
  for (char value : line) {
    if (value == '\t') ++count;
  }
  return count;
}

bool load_mapping(
    const std::string& path,
    std::size_t cell_count,
    int group_count,
    std::vector<int>& cell_to_group,
    std::size_t& selected_cells) {
  std::ifstream input(path);
  if (!input) return false;
  cell_to_group.assign(cell_count, -1);
  selected_cells = 0;
  std::string line;
  if (!std::getline(input, line)) return false;
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    std::istringstream stream(line);
    long long cell_index = 0;
    long long group_index = 0;
    if (!(stream >> cell_index >> group_index)) return false;
    if (cell_index < 1 || static_cast<std::size_t>(cell_index) > cell_count ||
        group_index < 1 || group_index > group_count) {
      return false;
    }
    std::size_t zero_cell = static_cast<std::size_t>(cell_index - 1);
    if (cell_to_group[zero_cell] != -1) return false;
    cell_to_group[zero_cell] = static_cast<int>(group_index - 1);
    ++selected_cells;
  }
  return selected_cells > 0;
}

long long parse_integer(const char* start) {
  errno = 0;
  char* end = nullptr;
  long long value = std::strtoll(start, &end, 10);
  if (errno != 0 || end == start) {
    throw std::runtime_error("invalid integer count");
  }
  return value;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc != 6) {
    std::cerr << "Usage: stream_tsv_pseudobulk MATRIX.tsv.gz MAPPING.tsv OUTPUT.tsv QC.tsv GROUP_COUNT\n";
    return 2;
  }
  const std::string matrix_path = argv[1];
  const std::string mapping_path = argv[2];
  const std::string output_path = argv[3];
  const std::string qc_path = argv[4];
  const int group_count = std::atoi(argv[5]);
  if (group_count < 1) {
    std::cerr << "GROUP_COUNT must be positive\n";
    return 2;
  }

  gzFile matrix = gzopen(matrix_path.c_str(), "rb");
  if (matrix == nullptr) {
    std::cerr << "Unable to open matrix: " << matrix_path << "\n";
    return 3;
  }
  std::string line;
  if (!gz_getline(matrix, line)) {
    std::cerr << "Matrix header is missing\n";
    gzclose(matrix);
    return 3;
  }
  const std::size_t fields = field_count(line);
  if (fields < 2) {
    std::cerr << "Matrix header has no cell columns\n";
    gzclose(matrix);
    return 3;
  }
  const std::size_t cell_count = fields - 1;
  std::vector<int> cell_to_group;
  std::size_t selected_cells = 0;
  if (!load_mapping(mapping_path, cell_count, group_count, cell_to_group, selected_cells)) {
    std::cerr << "Invalid aggregation mapping\n";
    gzclose(matrix);
    return 4;
  }

  std::ofstream output(output_path);
  if (!output) {
    std::cerr << "Unable to create output: " << output_path << "\n";
    gzclose(matrix);
    return 5;
  }
  output << "gene";
  for (int group = 1; group <= group_count; ++group) output << "\tcount_" << group;
  for (int group = 1; group <= group_count; ++group) output << "\tdetected_" << group;
  output << '\n';

  std::size_t feature_count = 0;
  try {
    while (gz_getline(matrix, line)) {
      if (line.empty()) continue;
      const std::size_t first_tab = line.find('\t');
      if (first_tab == std::string::npos) throw std::runtime_error("feature row has no counts");
      std::vector<long long> counts(group_count, 0);
      std::vector<long long> detected(group_count, 0);
      std::size_t position = first_tab + 1;
      for (std::size_t cell = 0; cell < cell_count; ++cell) {
        if (position > line.size()) throw std::runtime_error("feature row has too few columns");
        const long long value = parse_integer(line.c_str() + position);
        const int group = cell_to_group[cell];
        if (group >= 0) {
          counts[group] += value;
          if (value > 0) ++detected[group];
        }
        const std::size_t next_tab = line.find('\t', position);
        if (cell + 1 < cell_count) {
          if (next_tab == std::string::npos) throw std::runtime_error("feature row has too few columns");
          position = next_tab + 1;
        }
      }
      output << line.substr(0, first_tab);
      for (long long value : counts) output << '\t' << value;
      for (long long value : detected) output << '\t' << value;
      output << '\n';
      ++feature_count;
    }
  } catch (const std::exception& error) {
    std::cerr << "Aggregation error at feature " << (feature_count + 1) << ": " << error.what() << "\n";
    gzclose(matrix);
    return 6;
  }
  gzclose(matrix);
  output.close();

  std::ofstream qc(qc_path);
  if (!qc) {
    std::cerr << "Unable to create QC output: " << qc_path << "\n";
    return 5;
  }
  qc << "key\tvalue\n"
     << "features\t" << feature_count << '\n'
     << "matrix_cells\t" << cell_count << '\n'
     << "selected_cells\t" << selected_cells << '\n'
     << "groups\t" << group_count << '\n';
  return 0;
}
