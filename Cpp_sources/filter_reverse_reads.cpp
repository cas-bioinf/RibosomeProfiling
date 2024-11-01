// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2021-08-26
// Released under Apache License 2.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>

/// <summary>
/// Checks, whether a line has set a flag in FLAG field
/// </summary>
/// <param name="line">The examined line.</param>
/// <param name="flag">The examined flag.</param>
/// <returns>TRUE, if the flag is not set; FALSE if flag is set or if an unexpected line format occured.</returns>
inline bool check_flag(const std::string& line, const size_t flag) {
  size_t from = line.find('\t');
  if (from == line.npos) {
    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
    return false;
  }
  ++from;
  size_t to = line.find('\t', from);
  if (to == line.npos) {
    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
    return false;
  }
  return !(std::stoull(line.substr(from, to - from)) & flag);
}

/// <summary>
/// Checks, whether a SEQ is reverse complemented.
/// </summary>
/// <param name="line">The examined line.</param>
/// <returns>TRUE, if a flag 16 is not set in FLAG field; FALSE if the flag is set or if an unexpected line format occured.</returns>
inline bool forward_flag(const std::string &line) {
  return check_flag(line, 16);
}

/// <summary>
/// Checks, whether it is a primary alignment.
/// </summary>
/// <param name="line">The examined line;</param>
/// <returns>TRUE if a flag 256 is not set in FLAG field; FALSE if the flag is set or if an unexpected line format occured.</returns>
inline bool primary_flag(const std::string& line) {
  return check_flag(line, 256);
}

int main(int argc, char* argv[]) {
  if (argc == 1 || (argc % 2) != 1) {
    std::cout << "filter_reverse_reads (<input> <output>)+\t Takes <input> file in SAM format, filter out all reads that are mapped\n";
    std::cout << "                                        \t to reverse strand, and write the rest to <output> file.\n";
    std::cout << "                                        \t It expectes that the input file has grouped QNAMEs and that NH:i:Nmap\n";
    std::cout << "                                        \t is valid.\n";
    std::cout << "                                        \t It updates FLAG with respect by choosing a new primary alignment, MAPQ,\n";
    std::cout << "                                        \t NH:i:Nmap and HI:i:I.\n";
    std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2021-08-26; license: Apache License 2.0" << std::endl;
    return 0;
  }
  
  for (size_t argi = 1; argi < argc; argi += 2) { // Foreach pair of filenames
    std::ifstream input(argv[argi]);
    std::ofstream output(argv[argi + 1]);
    std::string line;
    while (std::getline(input, line)) {
      if (line.empty()) { // Not recognized
        std::cerr << "Unexpected empty line in file '" << argv[argi] << "'" << std::endl;
      } else if (line[0] == '@') { // Header; well it could be filtered to leave out missing transcript ids, but it would result in two-pass read, or store whole file in RAM
        output << line << '\n';
      } else { // Alignments
        // First we need to know, how many alignments there are for the current read
        size_t count = line.rfind("\tNH:i:");
        if (count == line.npos) {
          std::cerr << "Unexpected file format: missing NH:i: tag '" << line << "'" << std::endl;
        } else {
          count += 6;
          size_t count_end = line.find('\t', count);
          count = std::stoull(count_end == line.npos ? line.substr(count) : line.substr(count, count_end - count));

          if (count == 1) { // If there is just a single read, there is nothing to modify
            if (forward_flag(line)) {
              output << line << '\n';
            }
          } else { // Not unique alignment
            // Multiple lines must be processed together to correctly update NH:i tag, MAPQ score etc.
            std::vector<std::string> group;
            // Whether a preserved line is marked as a primary alignment
            size_t primary = -1;
            // Whether an alignment was leaved out
            bool modification = false;
            if (forward_flag(line)) {
              if (primary_flag(line)) {
                primary = group.size();
              }
              group.push_back(line);
            } else {
              modification = true;
            }
            for (size_t i = 1; i < count; i++) {
              if (std::getline(input, line)) {
                if (forward_flag(line)) {
                  if (primary_flag(line)) {
                    primary = group.size();
                  }
                  group.push_back(line);
                } else {
                  modification = true;
                }
              } else {
                std::cerr << "Unexpected end of file file '" << argv[argi] << "'" << std::endl;
                return 17;
              }
            }

            if (modification) { // Some alignment was leaved out, so other lines must be updated
              if (primary == -1) { // No preserved alignment is labeled as a primary, so one of them must be choosen
                // CIGAR strings describing the alignment
                std::vector<std::string> cigars;
                for (size_t i = 0; i < group.size(); i++) {
                  size_t index = 0;
                  for (size_t j = 0; j < 5; j++) {
                    index = group[i].find('\t', index + 1);
                    if (index == group[i].npos) {
                      std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
                      return 12;
                    }
                  }
                  size_t index2 = group[i].find('\t', index + 1);
                  if (index2 == group[i].npos) {
                    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
                    return 11;
                  }
                  cigars.push_back(group[i].substr(index, index2 - index));
                }
                // For now, all alignments have the some CIGAR so it is not necessary to try to identify alignment score
                for (size_t i = 1; i < cigars.size(); i++) {
                  if (cigars[i] != cigars[0]) {
                    std::cerr << "Not implemented yet '" << line << "'" << std::endl;
                  }
                }
                primary = 0;
              } else { // Some preserved alignment was marked as a primary, so no one flag must be changed
                primary = -1;
              }
              // It is constant for all alignments
              std::string mapq = std::to_string((size_t)(group.size() <= 1 ? 255 : (-10 * std::log10(1 - 1.0 / group.size()))));
              for (size_t i = 0; i < group.size(); i++) {
                size_t j = 0;
                std::stringstream parts(group[i]);
                std::string part;
                while (std::getline(parts, part, '\t')) {
                  if (++j == 2 && primary == i) { // New primary alignment must be changed
                    part = std::to_string(std::stoull(part) ^ 256);
                  } else if (j == 5) { // MAPQ score must be recomputed
                    part = mapq;
                  } else if (part.rfind("NH:i:", 0) == 0) { // Number of alignments was changed
                    part = "NH:i:" + std::to_string(group.size());
                  } else if (part.rfind("HI:i:", 0) == 0) { // Index of the alignment could be changed
                    part = "HI:i:" + std::to_string(i+1);
                  }
                  if (j != 1) { // Trimming separator
                    output << '\t';
                  }
                  output << part;
                }
                output << '\n';
              }
            } else { // No alignment leaved out, so no changes in lines
              for (size_t i = 0; i < group.size(); i++) {
                output << group[i] << '\n';
              }
            }
          }
        }
      }
    }
    // Cleaning
    input.close();
    output.flush();
    output.close();
  }
  return 0;
}
