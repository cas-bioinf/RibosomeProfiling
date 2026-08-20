// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2025-01-05
// Released under Apache License 2.0

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>

struct Settings {
public:
	std::ofstream output;
	char reverse_action;
	long long left_shift;
	bool left_at_least;
	bool left_over;
	long long right_shift;
	bool right_at_least;
	bool right_over;

	Settings() {};

	Settings(const std::string& output, const long long left, const std::string& left_mode, const long long right, const std::string& right_mode, const std::string& reverse_mode)
		: output(output), reverse_action(reverse_mode == "s" ? 0 : (reverse_mode == "r" ? 2 : 1)),
		left_shift(left), left_at_least(left_mode.find('m') == left_mode.npos), left_over(left_mode.find('e') == left_mode.npos),
		right_shift(right), right_at_least(right_mode.find('m') == right_mode.npos), right_over(right_mode.find('e') == right_mode.npos) {
		if (left_mode.size() != 2) {
			throw std::invalid_argument("'" + left_mode + "' is not a valid <left_mode> argument: unexpected length.");
		}
		if (left_at_least && left_mode.find('l') == left_mode.npos) {
			throw std::invalid_argument("'" + left_mode + "' is not a valid <left_mode> argument: at 'm'ost or at 'l'east switcher is missing.");
		}
		if (left_over && left_mode.find('i') == left_mode.npos) {
			throw std::invalid_argument("'" + left_mode + "' is not a valid <left_mode> argument: 'i'nclude or at 'e'xclude switcher is missing.");
		}
		if (right_mode.size() != 2) {
			throw std::invalid_argument("'" + right_mode + "' is not a valid <left_mode> argument: unexpected length.");
		}
		if (right_at_least && right_mode.find('l') == right_mode.npos) {
			throw std::invalid_argument("'" + right_mode + "' is not a valid <left_mode> argument: at 'm'ost or at 'l'east switcher is missing.");
		}
		if (right_over && right_mode.find('i') == right_mode.npos) {
			throw std::invalid_argument("'" + right_mode + "' is not a valid <left_mode> argument: 'i'nclude or at 'e'xclude switcher is missing.");
		}
		if (reverse_action == 1 && reverse_mode != "u") {
			throw std::invalid_argument("'" + reverse_mode + "' is not a valid <left_mode> argument: 's'kip, 'u'se and 'r'evert are the only supported options.");
		}
	}
};

inline bool check(const std::string& transcript, const std::string& chromosome, const bool direction, const std::vector<std::string>& line) {
	if ((!direction || line[3] != "+") && (direction || line[3] != "-")) {
		std::cerr << "Unexpected file format - unknown or inconsistent strand specifier: " << transcript;
		for (size_t i = 0; i < line.size(); ++i) {
			std::cerr << '\t' << line[i];
		}
		std::cerr << std::endl;
		return false;
	}
	if (chromosome != line[0]) {
		std::cerr << "Unexpected file format - inconsistent chromosome specifier: " << transcript;
		for (size_t i = 0; i < line.size(); ++i) {
			std::cerr << '\t' << line[i];
		}
		std::cerr << std::endl;
		return false;
	}
	return true;
}

int main(int argc, char* argv[]) {
	if (argc < 10 || argc % 7 != 3) {
    std::cout << "select_features <input> <annotations> (<feature> <left_shift> <left_mode> <right_shift> <right_mode> <complemented_mode> <output>)+\n";
    std::cout << "\tSelects reads from the input file in BAM format that cover a feature within a transcript and writes them to an output file.\n";
    std::cout << "\tThe program supports defining multiple features, each with its own output file.\n\n";
    std::cout << "Parameters:\n";
    std::cout << "\t<input>             Input file in SAM format.\n";
    std::cout << "\t<annotations>       Annotation file in GTF format.\n";
    std::cout << "\t<feature>           Name of the feature to search for (e.g. 'five_prime_utr', 'start_codon', 'CDS').\n";
    std::cout << "\t                    Note: Values 'gene', 'transcript', and 'exon' are not supported as filtering features.\n";
    std::cout << "\t<left_shift>\n";
    std::cout << "\t<right_shift>       Number of bases to shift from the start/ end of the feature. A positive integer enlarges the region, while a negative value shrinks it.\n";
    std::cout << "\t<left_mode>\n";
    std::cout << "\t<right_mode>        Mode for the left/ right boundary (2 characters):\n";
    std::cout << "\t                    - The first character: 'l' (at least, i.e. overlap the boundary) or 'm' (at most, i.e within the region only).\n";
    std::cout << "\t                    - The second character says what to do if the boundary of the region is shifted out of a transcript and a read starts at the first possible position:\n";
    std::cout << "\t                      'i' (include boundary) or 'e' (exclude boundary).\n";
    std::cout << "\t<complemented_mode> Action for reads on the reverse strand:\n";
    std::cout << "\t                    - 'p': Preserve (preserve original boundaries for reverse reads).\n";
    std::cout << "\t                    - 'r': Revert (swap left/right logic for reverse reads).\n";
    std::cout << "\t                    - 's': Skip (ignore reverse reads).\n";
    std::cout << "\t<output>            Output filename for the specific rule.\n\n";
    std::cout << "\tExample:\n";
    std::cout << "\tselect_features input.sam annotations.gtf five_prime_utr 1000 mi 5 mi s 5utr.sam start_codon 3 le 3 le s start.bam CDS 2 mi 2 mi s cds.sam stop_codon 6 le 0 le s stop.sam three_prime_utr 8 mi 1000 mi s 3utr.sam\n\n";
    std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2023-01-05; license: Apache License 2.0" << std::endl;
		return 0;
	}

	// transcript => [length; feature => [start; end]]
	std::map<std::string, std::pair<long long, std::map<std::string, std::pair<long long, long long>>>> selection;
	// Parse location of selected features
	{
		std::set<std::string> features;
		for (size_t i = 3; i < argc; i += 7) {
			std::string feature(argv[i]);
			if (feature == "gene" || feature == "transcript" || feature == "exon") {
				std::cerr << "'" << feature << "' is corrently not considered to be reasonable feature for this filtering." << std::endl;
				return 1;
			}
			features.emplace(feature);
		}

		std::string line;
		std::map<std::string, std::map<std::string, std::map<long long, std::vector<std::string>>>> lines;
		std::ifstream annotations(argv[2]);
		while (std::getline(annotations, line)) {
			if (line.empty()) {
				std::cerr << "An empty line within fle '" << argv[2] << "'." << std::endl;
			} else if (line[0] != '#') {
				// Parsing columns
				std::vector<std::string> parts;
				parts.reserve(8);

				size_t from = 0;
				for (size_t i = 0; i < 8; i++) {
					size_t to = line.find('\t', from);
					if (to == line.npos) {
						std::cerr << "Unexpected line format - not enough columns: " << line << std::endl;
						break;
					}
					parts.push_back(line.substr(from, to - from));
					from = to + 1;
				}
				if (parts.size() != 8) {
					continue;
				}
				if (line.find('\t', from) != line.npos) {
					std::cerr << "Unexpected line format - too many columns: " << line << std::endl;
					continue;
				}
				
				if (parts[2] == "exon" || features.find(parts[2]) != features.end()) {
					from = line.find("transcript_id \"", from);
					if (from == line.npos) {
						std::cerr << "Unexpected line format - missing transcript_id tag: " << line << std::endl;
					}
					from += 15;
					size_t to = line.find("\";", from);
					if (to == line.npos) {
						std::cerr << "Unexpected line format - unfinished transcript_id tag: " << line << std::endl;
					}
					std::string transcript = line.substr(from, to - from);
					std::vector<std::string>& preserve = lines[transcript][parts[2]][parts[6] == "+" ? std::stoll(parts[3]) : 0 - std::stoll(parts[4])];
					preserve.reserve(4);
					preserve.push_back(parts[0]);
					preserve.push_back(parts[3]);
					preserve.push_back(parts[4]);
					preserve.push_back(parts[6]);
				}
			}
		}

		for (auto transcripts_it = lines.begin(); transcripts_it != lines.end(); ++transcripts_it) {
			// end => [start; shift]
			std::map<long long, std::pair<long long, long long>> intervals;
			std::string chromosome = transcripts_it->second["exon"].begin()->second[0];
			bool direction = transcripts_it->second["exon"].begin()->second[3] == "+";
			long long sum = 0;
			for (auto lines_it = transcripts_it->second["exon"].begin(); lines_it != transcripts_it->second["exon"].end(); ++lines_it) {
				if (!check(transcripts_it->first, chromosome, direction, lines_it->second)) {
					continue;
				}
				long long start = std::stoll(lines_it->second[1]);
				long long end = std::stoll(lines_it->second[2]);
				intervals[end] = std::pair<long long, long long>(start, sum);
				sum += end - start + 1;
			}
			selection[transcripts_it->first] = std::pair<long long, std::map<std::string, std::pair<long long, long long>>>(sum, std::map<std::string, std::pair<long long, long long>>());

			for (auto features_it = features.begin(); features_it != features.end(); ++features_it) {
				if (transcripts_it->second.find(*features_it) == transcripts_it->second.end()) {
					std::cerr << "Transcript '" << transcripts_it->first << "' does not contains any line with the feature '" << *features_it << "'." << std::endl;
					continue;
				}
				long long start = 0;
				long long end = 0;
				long long length = 0;
				for (auto lines_it = transcripts_it->second[*features_it].begin(); lines_it != transcripts_it->second[*features_it].end(); ++lines_it) {
					if (!check(transcripts_it->first, chromosome, direction, lines_it->second)) {
						continue;
					}
					long long start_curr = std::stoll(lines_it->second[1]);
					long long end_curr = std::stoll(lines_it->second[2]);
					long long length_curr = end_curr - start_curr + 1;
					auto intervals_it = intervals.lower_bound(start_curr);
					if (intervals_it == intervals.end()) {
						std::cerr << "Unexpected file format - feature is outside an exon region: " << line << std::endl;
						continue;
					}
					if (direction) {
						start_curr = 1 + intervals_it->second.second + start_curr - intervals_it->second.first;
					} else {
						start_curr = 1 + intervals_it->second.second + intervals_it->first - end_curr;
					}
					end_curr = start_curr + length_curr - 1;
					if (length == 0) {
						start = start_curr;
						end = end_curr;
						length = length_curr;
					} else {
						if (end + 1 != start_curr) {
							std::cerr << "Unexpected file format - feature is discontinued: " << transcripts_it->first;
							for (auto i = 0; i < lines_it->second.size(); ++i) {
								std::cerr << '\t' << lines_it->second[i];
							}
							std::cerr << std::endl;
							continue;
						}
						end = end_curr;
						length += length_curr;
					}
				}
				selection[transcripts_it->first].second[*features_it] = std::pair<long long, long long>(start, end);
			}
		}
	}

	// Filter;
	{
		std::map<std::string, size_t> stats;
		std::map<std::string, Settings> outputs;
		for (size_t i = 3; i < argc; i += 7) {
			outputs.emplace(argv[i], Settings(argv[i + 6], std::stoll(argv[i + 1]), argv[i + 2], std::stoll(argv[i + 3]), argv[i + 4], argv[i + 5]));
		}
		std::ifstream input(argv[1]);
		std::string line;
		while (std::getline(input, line)) {
			if (line.empty()) {
				std::cerr << "An empty line within fle '" << argv[1] << "'." << std::endl;
			} else if (line[0] == '@') {
				for (auto outputs_it = outputs.begin(); outputs_it != outputs.end(); ++outputs_it) {
					outputs_it->second.output << line << '\n';
				}
			} else {
				std::vector<std::string> lines;
				std::set<std::string> matches;
				bool matches_init = true;
				size_t count = -1;
				bool unknown = false;
				bool inconsistent = false;
				for (size_t hi = 0; hi < count; ) {
					lines.push_back(line);
					std::vector<std::string> parts;
					std::string part;
					std::stringstream stream(line);
					while (std::getline(stream, part, '\t')) {
						parts.push_back(part);
					}
					if (hi == 0) {
						for (size_t i = 0; i < parts.size(); i++) {
							if (parts[i].rfind("NH:i:", 0) == 0) {
								count = std::stoull(parts[i].substr(5));
							}
						}
						if (count == -1) {
							std::cerr << "An unexpected line format of the input file, missing 'NH:i:Nmap' field: " << line << std::endl;
							break;
						}
					}
					if (parts.size() < 11) {
						std::cerr << "An unexpected line format of the input file, not enough columns: " << line << std::endl;
						continue;
					}

					auto selection_it = selection.find(parts[2]);
					if (selection_it == selection.end()) {
						unknown = true;
					} else {
						bool reverse = (std::stoull(parts[1]) & 16) > 0;
						long long start = std::stoll(parts[3]);
						long long end = start - 1;
						long long length = 0;
						for (size_t i = 0; i < parts[5].size(); ++i) {
							if (parts[5][i] >= '0' && parts[5][i] <= '9') {
								length *= 10;
								length += parts[5][i] - '0';
							} else {
								switch (parts[5][i]) {
								case 'M':
								case 'D':
								case 'N':
								case '=':
								case 'X':
									end += length;
									length = 0;
									break;
								case 'I':
								case 'S':
								case 'H':
								case 'P':
									length = 0;
									break;
								default:
									std::cerr << "An unsupported CIGAR code '" << parts[5][i] << "': " << line << std::endl;
									break;
								}
							}
						}
						if (length > 0) {
							std::cerr << "An unfinished CIGAR string '" << parts[5] << "': " << line << std::endl;
							continue;
						}
						std::set<std::string> ms;
						for (auto features_it = selection_it->second.second.begin(); features_it != selection_it->second.second.end(); ++features_it) {
							auto& output = outputs[features_it->first];
							// Reverse strands should be filtered out
							if (output.reverse_action == 0 && reverse) {
								continue;
							}
							bool revert = reverse && output.reverse_action == 2;

							// Check the beginning of the read
							if (start <= 1 && features_it->second.first < 1 + (revert ? output.right_shift : output.left_shift)) {
								// Check whether touching the beginning of the transcript is allowed
								if (!(revert ? output.right_over : output.left_over)) {
									continue;
								}
								// Elsewhere it is valid end
							} else if (revert ? output.right_at_least : output.left_at_least) { // At least variant
																																									// Check whether read starts before the boundary
								if (start + (revert ? output.right_shift : output.left_shift) > features_it->second.first) {
									continue;
								}
							} else { // At most variant
											 // Check whether read starts after the boundary
								if (start + (revert ? output.right_shift : output.left_shift) < features_it->second.first) {
									continue;
								}
							}

							// Check the end of the read
							if (end >= selection_it->second.first && features_it->second.second + (revert ? output.left_shift : output.right_shift) > selection_it->second.first) {
								// Check whether touching the beginning of the transcript is allowed
								if (!(revert ? output.left_over : output.right_over)) {
									continue;
								}
								// Elsewhere it is valid end
							} else if (revert ? output.left_at_least : output.right_at_least) { // At least variant
																																									// Check whether read ends after the boundary
								if (end < features_it->second.second + (revert ? output.left_shift : output.right_shift)) {
									continue;
								}
							} else { // At most variant
											 // Check whether read ends before the boundary
								if (end > features_it->second.second + (revert ? output.left_shift : output.right_shift)) {
									continue;
								}
							}

							// Passed
							ms.emplace(features_it->first);
						}
						if (matches_init) {
							matches = ms;
							matches_init = false;
						} else if (matches != ms) {
							inconsistent = true;
							if (matches.size() < ms.size()) {
								matches = ms;
							}
						}
					}
					if (++hi < count && !std::getline(input, line)) {
						std::cerr << "Unexpected end of input file, the last line is: '" << line << "'." << std::endl;
						matches.clear();
						return 45;
					}
				}
				switch (matches.size()) {
				case 0:
					if (inconsistent) {
						std::cerr << "Unexpected state: ambiguous features, but at most 0 features at once. The first line: '" << lines.front() << "'." << std::endl;
					} else {
						if (unknown) {
							++stats[matches_init ? "Unknown transcript" : "No feature-Unknown transcript"];
							std::cerr << (matches_init ? "Unknown transcript:\t" : "No feature-Unknown transcript:\t") << lines.front() << std::endl;
						} else {
							++stats["No feature"];
							std::cerr << "No feature:\t" << lines.front() << std::endl;
						}
					}
					break;
				case 1:
					if (unknown) {
						if (inconsistent) {
							++stats["Inconsistent feature-Unknown transcript"];
							std::cerr << "Inconsistent feature-Unknown transcript:\t" << lines.front() << std::endl;
						} else {
							++stats["Unknown transcript"];
							std::cerr << "Unknown transcript:\t" << lines.front() << std::endl;
						}
					} else {
						if (inconsistent) {
							++stats["Inconsistent feature"];
							std::cerr << "Inconsistent feature:\t" << lines.front() << std::endl;
						} else {
							++stats[*matches.begin()];
							auto output = outputs.find(*matches.begin());
							for (auto lines_it = lines.begin(); lines_it != lines.end(); ++lines_it) {
								output->second.output << *lines_it << '\n';
							}
						}
					}
					break;
				default:
					if (unknown) {
						if (inconsistent) {
							++stats["Ambiguous feature-Inconsistent feature-Unknown transcript"];
							std::cerr << "Ambiguous feature-Inconsistent feature-Unknown transcript:\t" << lines.front() << std::endl;
						} else {
							++stats["Ambiguous feature-Unknown transcript"];
							std::cerr << "Ambiguous feature-Unknown transcript:\t" << lines.front() << std::endl;
						}
					} else {
						if (inconsistent) {
							++stats["Ambiguous feature-Inconsistent feature"];
							std::cerr << "Ambiguous feature-Inconsistent feature:\t" << lines.front() << std::endl;
						} else {
							++stats["Ambiguous feature"];
							std::cerr << "Ambiguous feature:\t" << lines.front() << std::endl;
						}
					}
					for (auto matches_it = matches.begin(); matches_it != matches.end(); ++matches_it) {
						std::cerr << *matches_it << std::endl;
					}
					break;
				}
			}
		}
		input.close();
		for (auto outputs_it = outputs.begin(); outputs_it != outputs.end(); ++outputs_it) {
			outputs_it->second.output.flush();
			outputs_it->second.output.close();
		}
		for (auto stats_it = stats.begin(); stats_it != stats.end(); ++stats_it) {
			std::cout << stats_it->first << '\t' << stats_it->second << '\n';
		}
	}
}
