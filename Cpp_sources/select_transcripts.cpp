// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2021-09-01
// Released under Apache License 2.0

#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <vector>
#include <cmath>

/// <summary>
/// Extracts template sequence id (transcript) - the third column;
/// </summary>
/// <param name="line">The examined line;</param>
/// <returns>Transcript id</returns>
inline std::string extract_transcript_id(const std::string& line) {
  size_t from = line.find('\t');
  if (from == line.npos) {
	std::cerr << "Unexpected line format: not enough columns: '" << line << "'." << std::endl;
	return "";
  }
  size_t to = line.find('\t', ++from);
  if (to == line.npos) {
	std::cerr << "Unexpected line format: not enough columns: '" << line << "'." << std::endl;
	return "";
  }
  from = ++to;
  to = line.find('\t', from);
  if (to == line.npos) {
	std::cerr << "Unexpected line format: not enough columns: '" << line << "'." << std::endl;
	return "";
  }
  return line.substr(from, to - from);
}

/// <summary>
/// Checks, whether it is a primary alignment.
/// </summary>
/// <param name="line">The examined line;</param>
/// <returns>TRUE if a flag 256 is not set in FLAG field; FALSE if the flag is set or if an unexpected line format occured.</returns>
inline bool primary_flag(const std::string& line) {
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
  return !(std::stoull(line.substr(from, to - from)) & 256);
}


int main(int argc, char* argv[]) {
  if (argc < 4 || argc % 2 != 0) {
	std::cout << "select_transcripts <transcript_ids> <input> <output>\t Filters <input> SAM file only for transcripts from\n";
	std::cout << "                                                    \t <transcript_ids> file (one id per line) and store them in\n";
	std::cout << "                                                    \t <output> SAM file.\n";
	std::cout << "                                                    \t '@SQ', Flags, MAPq, 'NH:i:Nmap' and 'HI:i:id' fileds are\n";
	std::cout << "                                                    \t updated.\n";
	std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2021-09-01; license: Apache License 2.0" << std::endl;
	return 0;
  }

  std::set<std::string> transcript_ids;
  { // Load, what transcript_ids should be preserved
	std::ifstream input(argv[1]);
	for (std::string line; std::getline(input, line); ) {
	  transcript_ids.insert(line);
	}
	input.close();
  }

  // Filter input files
  for (size_t argi = 2; argi < argc; argi += 2) {
	std::ifstream input(argv[argi]);
	std::ofstream output(argv[argi + 1]);
	for (std::string line; std::getline(input, line); ) {
	  if (line.empty()) {
		std::cerr << "Unexpected empty line in file '" << argv[argi] << "'." << std::endl;
		continue;
	  } if (line[0] == '@') {
		if (line.rfind("@SQ\t", 0) == 0) {
		  size_t from = line.find("\tSN:");
		  if (from == line.npos) {
			std::cerr << "Unexpected line format: missing SN field within '@SQ' line: '" << line << "'." << std::endl;
			continue;
		  }
		  from += 4;
		  size_t to = line.find('\t', from);
		  if (transcript_ids.find(to == line.npos ? line.substr(from) : line.substr(from, to - from)) != transcript_ids.end()) {
			output << line << '\n';
		  }
		} else {
		  output << line << '\n';
		}
	  } else {
		// First we need to know, how many alignments there are for the current read
		size_t count = line.rfind("\tNH:i:");
		if (count == line.npos) {
		  std::cerr << "Unexpected file format: missing NH:i: tag '" << line << "'" << std::endl;
		  continue;
		} else {
		  count += 6;
		  size_t count_end = line.find('\t', count);
		  count = std::stoull(count_end == line.npos ? line.substr(count) : line.substr(count, count_end - count));

		  if (count == 1) { // If there is just a single read, there is nothing to modify
			if (transcript_ids.find(extract_transcript_id(line)) != transcript_ids.end()) {
			  output << line << '\n';
			}
		  } else { // Not unique alignment
						   // Multiple lines must be processed together to correctly update NH:i tag, MAPQ score etc.
			std::vector<std::string> group;
			// Whether a preserved line is marked as a primary alignment
			size_t primary = -1;
			if (transcript_ids.find(extract_transcript_id(line)) != transcript_ids.end()) {
			  if (primary_flag(line)) {
				primary = group.size();
			  }
			  group.push_back(line);
			}
			for (size_t i = 1; i < count; i++) {
			  if (std::getline(input, line)) {
				if (transcript_ids.find(extract_transcript_id(line)) != transcript_ids.end()) {
				  if (primary_flag(line)) {
					primary = group.size();
				  }
				  group.push_back(line);
				}
			  } else {
				std::cerr << "Unexpected end of file file '" << argv[argi] << "'" << std::endl;
				return 17;
			  }
			}

			if (group.size() == count) { // No alignment leaved out, so no changes in lines
			  for (size_t i = 0; i < group.size(); i++) {
				output << group[i] << '\n';
			  }
			} else { // Some alignment was leaved out, so other lines must be updated
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
					part = "HI:i:" + std::to_string(i + 1);
				  }
				  if (j != 1) { // Trimming separator
					output << '\t';
				  }
				  output << part;
				}
				output << '\n';
			  }
			}
		  }
		}
	  }
	}
	input.close();
	output.flush();
	output.close();
  }

  return 0;
}
