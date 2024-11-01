// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2021-08-26
// Released under Apache License 2.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include "filter_ambiguous_genes.h"

/// <summary>
/// Returns reference sequence name for the given alignment.
/// </summary>
/// <param name="line">The examined line.</param>
/// <returns>RNAME, or empty string if an error occures.</returns>
inline std::string extract_transcript_id(const std::string& line) {
  size_t transcript_from = line.find('\t');
  if (transcript_from == line.npos) {
    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
    return "";
  }
  transcript_from = line.find('\t', transcript_from + 1);
  if (transcript_from == line.npos) {
    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
    return "";
  }
  size_t transcript_to = line.find('\t', ++transcript_from);
  if (transcript_to == line.npos) {
    std::cerr << "Unexpected file format: not enough columns '" << line << "'" << std::endl;
    return "";
  }
  return line.substr(transcript_from, transcript_to - transcript_from);
}

int main(int argc, char* argv[]) {
  if (argc == 1 || (argc % 2) != 0) {
    std::cout << "filter_ambiguous_genes <annotations> (<input> <output>)+\t It takes transcript_id => gene_id mapping from\n";
    std::cout << "                                                        \t <annotations> file in GTF format and then it read\n";
    std::cout << "                                                        \t <input> file in SAM format, filter out all reads that\n";
    std::cout << "                                                        \t are mapped into multiple transcripts from different\n";
    std::cout << "                                                        \t genes (multiple transcripts from the same gene are\n";
    std::cout << "                                                        \t allowed), and write the rest to <output> file.\n";
    std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2021-08-26; license: Apache License 2.0" << std::endl;
    return 0;
  }

  // Mapping saying, what gene_id corresponds to a given transcript_id
  std::map<std::string, std::string> transcript_gene;
  {
    std::string line;
    std::ifstream input(argv[1]);
    while (std::getline(input, line)) {
      if (!line.empty() && line[0] != '#') {
        size_t transcript_from = line.find(" transcript_id \"");
        if (transcript_from != line.npos) {
          transcript_from += 16;
          size_t transcript_to = line.find("\";", transcript_from);
          if (transcript_to == line.npos) {
            std::cerr << "Unexpected line format: incomplete 'transcript_id' tag: '" << line << "'" << std::endl;
            return 2;
          }
          size_t gene_from = line.find("\tgene_id \"");
          if (gene_from == line.npos) {
            std::cerr << "Unexpected line format: missing 'gene_id' tag: '" << line << "'" << std::endl;
            return 2;
          }
          gene_from += 10;
          size_t gene_to = line.find("\";", gene_from);
          if (gene_to == line.npos) {
            std::cerr << "Unexpected line format: incomplete 'gene_id' tag: '" << line << "'" << std::endl;
            return 2;
          }
          transcript_gene[line.substr(transcript_from, transcript_to - transcript_from)] = line.substr(gene_from, gene_to - gene_from);
        }
      }
    }
    input.close();
  }

  for (size_t argi = 2; argi < argc; argi += 2) { // Foreach pair of filenames
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
          if (count == 1) { // If there is just a single read, there is nothing to check
            output << line << '\n';
          } else { // Not unique alignment, it is necessary to check, whether all alignments correspond to the same gene
            std::string transcript_id = extract_transcript_id(line);
            auto transcript_gene_it = transcript_gene.find(transcript_id);
            if (transcript_gene_it == transcript_gene.end()) {
              std::cerr << "Unknown gene_id: a transcript_id '" << transcript_id << "' did not occure in the annotations file '" << argv[1] << "': '" << line << "'" << std::endl;
              return 6;
            }


            // Multiple lines must be processed together to check unambiguity of corresponding gene_ids
            std::vector<std::string> group;
            group.push_back(line);
            for (size_t i = 1; i < count; i++) {
              if (std::getline(input, line)) {
                transcript_id = extract_transcript_id(line);
                auto it = transcript_gene.find(transcript_id);
                if (it == transcript_gene.end()) {
                  std::cerr << "Unknown gene_id: a transcript_id '" << transcript_id << "' did not occure in the annotations file '" << argv[1] << "': '" << line << "'" << std::endl;
                  return 7;
                }
                if (transcript_gene_it->second == it->second) {
                  group.push_back(line);
                } else {
                  group.clear();
                }
              } else {
                std::cerr << "Unexpected end of file file '" << argv[argi] << "'" << std::endl;
                return 17;
              }
            }

            if (group.size() == count) { // All lines was added
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
