// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2023-03-31
// Released under Apache License 2.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>

int main(int argc, char* argv[]) {
  if (argc != 3) {
	std::cout << "region_readcounts <ranges> <counts>\t Reads ranges [from; to) or lengths for each identifier from <ranges> in tab-separated values file format; and\n";
	std::cout << "                                   \t computes an total read count within the region from <counts> file in tab-separated values file format.\n\n";
	std::cout << "                                   \t <ranges> should have lines in format '[identifier]\\t[from]\\t[to]' or '[identifier]\\t[length]'; and\n";
	std::cout << "                                   \t <counts> should have lines in format '[identifier]\\t[position]\\t[count]'.\n";
	std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2023-03-31; license: Apache License 2.0" << std::endl;
	return 0;
  }

  // <id; <from; to>> Boundaries of the examined region [from; to) for each identifier
  std::map<std::string, std::pair<size_t, size_t>> ranges;
  {
	// Tab-separated input file with lines in '<id>\t<length>' or '<id>\t<from>\t<to>' format
	std::ifstream ranges_file(argv[1]);
	for (std::string line; std::getline(ranges_file, line);) {
	  // Position of the first separator
	  size_t tab_first = line.find('\t');
	  if (tab_first == line.npos) {
		std::cerr << "Unexpected line format, three columns expected, but only one occured: " << line << std::endl;
		continue;
	  }
	  // Position of the second separator
	  size_t tab_second = line.find('\t', tab_first + 1);
	  if (tab_second == line.npos) {
		ranges[line.substr(0, tab_first)] = std::pair<size_t, size_t>(1, 1+std::stoull(line.substr(tab_first + 1)));
	  } else {
		if (line.find('\t', tab_second + 1) != line.npos) {
		  std::cerr << "Unexpected line format, three columns expected, but at least four occured: " << line << std::endl;
		  continue;
		}
		ranges[line.substr(0, tab_first)] = std::pair<size_t, size_t>(std::stoull(line.substr(tab_first + 1, tab_second - tab_first - 1)), std::stoull(line.substr(tab_second + 1)));
	  }
	}
  }

  // <id, coef> Total read count within a region for each identifier
  std::map<std::string, double> coefs;
  {
	// Tab-separated input file with lines in '<id>\t<position>\t<count>' format
	std::ifstream counts_file(argv[2]);
	// Set of identifiers missing in ranges to do not repat the error message
	std::set<std::string> missing;
	for (std::string line; std::getline(counts_file, line);) {
	  // Position of the first separator
	  size_t tab_first = line.find('\t');
	  if (tab_first == line.npos) {
		std::cerr << "Unexpected line format, three columns expected, but only one occured: " << line << std::endl;
		continue;
	  }
	  // Position of the second separator
	  size_t tab_second = line.find('\t', tab_first + 1);
	  if (tab_second == line.npos) {
		std::cerr << "Unexpected line format, three columns expected, but only two occured: " << line << std::endl;
		continue;
	  }
	  if (line.find('\t', tab_second + 1) != line.npos) {
		std::cerr << "Unexpected line format, three columns expected, but at least four occured: " << line << std::endl;
		continue;
	  }
	  // The current identifier
	  std::string id = line.substr(0, tab_first);
	  // Pointer to the current identifier in ranges map
	  auto ranges_it = ranges.find(id);
	  if (ranges_it == ranges.end()) {
		if (missing.find(id) == missing.end()) {
		  std::cerr << "Identifier '" << id << "' is missing in the ranges file" << std::endl;
		  missing.emplace(id);
		}
		continue;
	  }
	  // Position from the current line
	  size_t pos = std::stoull(line.substr(tab_first + 1, tab_second - tab_first - 1));
	  if (ranges_it->second.first <= pos && pos < ranges_it->second.second) {
		coefs[id] += std::stod(line.substr(tab_second + 1));
	  }
	}
  }

  if (coefs.empty()) {
	std::cerr << "No coefficent was loaded, it is not possible to normalize" << std::endl;
	return 2;
  }
  std::cout.precision(10);
  for (auto coefs_it = coefs.begin(); coefs_it != coefs.end(); ++coefs_it) {
	std::cout << coefs_it->first << '\t' << coefs_it->second << '\n';
  }

  return 0;
}
