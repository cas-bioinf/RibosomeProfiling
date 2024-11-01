// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2023-03-29
// Released under Apache License 2.0

#include <iostream>
#include <sstream>
#include <map>

int main(int argc, char* argv[]) {
  if (argc != 1) {
    std::cout << "read_counts\t Read file in SAM format from standard input, group reads by RNAME and POS and\n";
    std::cout << "           \t print read counts to standard output in TAB-separated values file format.\n";
    std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2023-03-29; license: Apache License 2.0" << std::endl;
    return 0;
  }

  std::map<std::string, std::map<size_t, size_t>> counts;
  for (std::string line; std::getline(std::cin, line); ) {
    if (line.size() > 0 && line[0] != '@') {
      std::stringstream parts(line);
      std::string part;
      std::string rname;
      size_t pos;
      for (size_t i = 0; i < 4; i++) {
        if (!std::getline(parts, part, '\t')) {
          std::cerr << "Unexpected line format - not enough columns: " << line << std::endl;
          part = "";
          break;
        }
        if (i == 2) {
          rname = part;
        } else if (i == 3) {
          pos = std::stoull(part);
        }
      }
      if (part.empty()) {
        continue;
      }
      ++counts[rname][pos];
    }
  }
  for (auto counts_it = counts.begin(); counts_it != counts.end(); ++counts_it) {
    for (auto counts_it_it = counts_it->second.begin(); counts_it_it != counts_it->second.end(); ++counts_it_it) {
      std::cout << counts_it->first << '\t' << counts_it_it->first << '\t' << counts_it_it->second << '\n';
    }
  }
}
