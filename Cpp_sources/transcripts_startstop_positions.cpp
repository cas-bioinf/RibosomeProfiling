// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2023-03-29
// Released under Apache License 2.0

#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>

class Transcript {
private:
	std::string id;
	bool strand;
	std::map<size_t, size_t> exons;
	size_t start_codon;
	size_t stop_codon;
	bool error;

public:
	static const std::pair<size_t, size_t> UNDEFINED;

	/// <summary>
	/// Initialize the current transcript.
	/// </summary>
	/// <param name="transcript_id">Identifier of the transcript.</param>
	/// <param name="strand">Strand direction of the transcript.</param>
	Transcript(const std::string transcript_id, const bool strand) : id(transcript_id), strand(strand), start_codon(0), stop_codon(0), error(false) {}

	/// <summary>
	/// Provides an access to the transcript_id.
	/// </summary>
  /// <returns>Returns the current transcript_id.</returns>
	inline const std::string& transcript_id() const { return id; }

	/// <summary>
	/// Add an interval to the list of exons.
	/// For the simplicity, only start index of the exon is checked for duplicity; the rest will be checked during coordinates extraction.
	/// </summary>
	/// <param name="from">Start position of the exon.</param>
	/// <param name="to">Stop position of the exon.</param>
	inline void add_exon(const size_t from, const size_t to) {
		size_t from_ = strand ? from : -to;
		size_t to_ = strand ? to : -from;
		if (from > to) {
			std::cerr << "Transcript '" << id << "' contains a line with unordered start-stop positions: " << from << ", " << to << std::endl;
			error = true;
			return;
		}
		auto exons_it = exons.find(from_);
		if (exons_it != exons.end()) {
			std::cerr << "Transcript '" << id << "' contains overlapping exons at position " << (strand ? from : to) << std::endl;
			error = true;
			return;
		}
		exons[from_] = to_;
	}

	/// <summary>
	/// Compares a strand direction with the strand direction of the current transcript
	/// </summary>
	/// <param name="strand">Strand direction.</param>
	/// <returns>Returns whether strains have the same orientation.</returns>
	inline bool check_strand(const bool strand) {
		return this->strand == strand;
	}

	/// <summary>
	/// Update start codon position
	/// </summary>
	/// <param name="from">Start position of the fragment.</param>
	/// <param name="to">Stop position of the fragment.</param>
	inline void update_start_codon(const size_t from, const size_t to) {
		if (this->start_codon == 0) {
			this->start_codon = strand ? from : -to;
		} else {
			this->start_codon = std::min(this->start_codon, strand ? from : -to);
		}
	}

	/// <summary>
	/// Update stop codon position
	/// </summary>
	/// <param name="from">Start position of the fragment.</param>
	/// <param name="to">Stop position of the fragment.</param>
	inline void update_stop_codon(const size_t from, const size_t to) {
		if (this->stop_codon == 0) {
			this->stop_codon = strand ? from : -to;
		} else {
			this->stop_codon = std::min(this->stop_codon, strand ? from : -to);
		}
	}

	/// <summary>
	/// Performs rest of checks on the transcript and identify start and stop codon position in transcript coordinates.
	/// </summary>
	/// <returns>Returns start and stop codon 1-based coordinates.</returns>
	std::pair<size_t, size_t> get_coordinates() {
		if (error || id.empty()) {
			return UNDEFINED;
		}
		if (start_codon == 0) {
			std::cerr << "Transcript '" << id << "' does not have defined start_codon" << std::endl;
			error = true;
			return UNDEFINED;
		}
		if (stop_codon == 0) {
			std::cerr << "Transcript '" << id << "' does not have defined stop_codon" << std::endl;
			error = true;
			return UNDEFINED;
		}
		if (start_codon > stop_codon) {
			std::cerr << "Start and stop codons have the wrong order in transcript '" << id << "'" << std::endl;
			error = true;
			return UNDEFINED;
		}
		if (exons.empty()) {
			std::cerr << "No exon defined for transcript '" << id << "'" << std::endl;
			error = true;
			return UNDEFINED;
		}
		{
			size_t last = 0;
			for (auto exons_it = exons.begin(); exons_it != exons.end(); ++exons_it) {
				if (last >= exons_it->first) {
					std::cerr << "No exon defined for transcript '" << id << "'" << std::endl;
					error = true;
					return UNDEFINED;
				}
				last = exons_it->second;
			}
		}
		size_t start_position = 0;
		for (auto exons_it = exons.begin(); exons_it != exons.end(); ++exons_it) {
			if (start_codon < exons_it->first) {
				std::cerr << "Transcript '" << id << "' has start_codon outside exons" << std::endl;
				error = true;
				return UNDEFINED;
			}
			if (start_codon <= exons_it->second) {
				size_t stop_position = start_position;
				start_position += start_codon - exons_it->first + 1;
				do {
					if (stop_codon < exons_it->first) {
						std::cerr << "Transcript '" << id << "' has stop_codon outside exons" << std::endl;
						error = true;
						return UNDEFINED;
					}
					if (stop_codon <= exons_it->second) {
						stop_position += stop_codon - exons_it->first + 1;
						return std::pair<size_t, size_t>(start_position, stop_position);
					}
					stop_position += exons_it->second - exons_it->first + 1;
				} while (++exons_it != exons.end());
				std::cerr << "Transcript '" << id << "' has stop_codon outside exons" << std::endl;
				error = true;
				return UNDEFINED;
			}
			start_position += exons_it->second - exons_it->first + 1;
		}
		std::cerr << "Transcript '" << id << "' has start_codon outside exons" << std::endl;
		error = true;
		return UNDEFINED;
	}
};
/// <summary>
/// Default pair to be returned if something went wrong.
/// </summary>
const std::pair<size_t, size_t> Transcript::UNDEFINED(-1, -1);

int main(int argc, char* argv[]) {
	if (argc != 2) {
		std::cout << "transcripts_startstop_positions <GTF_file>\t Parses annotations file in GTF format and identifies start and stop codon positions for each transcript\n";
		std::cout << "                                          \t in coordinates relative to the transcript.\n";
		std::cout << "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2023-03-29; license: Apache License 2.0" << std::endl;
		return 0;
	}

	std::map<std::string, std::pair<size_t, size_t>> coordinates;
	{
		std::ifstream file(argv[1]);
		Transcript transcript("", false);
		for (std::string line; std::getline(file, line); ) {
			if (!line.empty() && line[0] != '#') {
				std::stringstream parts(line);
				std::string part("");
				std::string type("");
				size_t start = 0;
				size_t end = 0;
				bool strand = false;
				std::string transcript_id;
				for (size_t i = 0; i < 9; i++) {
					if (!std::getline(parts, part, '\t')) {
						std::cerr << "Unexpected line format - not enough columns: " << line << std::endl;
						part = "";
						break;
					}
					// Well, switch would require an extra variable or goto or needless iterations
					if (i == 2) {
						if (part == "exon" || part == "start_codon" || part == "stop_codon") {
							type = part;
						} else {
							break;
						}
					} else if (i == 3) {
						start = std::stoull(part);
					} else if (i == 4) {
						end = std::stoull(part);
					} else if (i == 6) {
						if (part == "+") {
							strand = true;
						} else if (part != "-") {
							std::cerr << "Unexpected or unsupported strand identifier '" << part << "' within line: " << line << std::endl;
							type = "";
							break;
						}
					} else if (i == 8) {
						size_t from = part.find("transcript_id \"");
						if (from == part.npos) {
							std::cerr << "Missing transcript_id attribute: " << line << std::endl;
							type = "";
							break;
						}
						from += 15;
						size_t to = part.find('"', from);
						if (to == part.npos) {
							std::cerr << "Unfinished transcript_id attribute: " << line << std::endl;
							type = "";
							break;
						}
						transcript_id = part.substr(from, to - from);
						if (transcript_id.empty()) {
							std::cerr << "Empty transcript_id attribute: " << line << std::endl;
							type = "";
							break;
						}
					}
				}
				if (!type.empty()) {
					if (transcript.transcript_id() != transcript_id) {
						coordinates[transcript.transcript_id()] = transcript.get_coordinates();
						transcript = Transcript(transcript_id, strand);
					} else if (!transcript.check_strand(strand)) {
						std::cerr << "Ambiguous strand for transcript '" << transcript.transcript_id() << "'" << std::endl;
						continue;
					}
					if (type == "exon") {
						transcript.add_exon(start, end);
					} else if (type == "start_codon") {
						transcript.update_start_codon(start, end);
					} else if (type == "stop_codon") {
						transcript.update_stop_codon(start, end);
					}
				}
			}
		}
		if (!transcript.transcript_id().empty()) {
			coordinates[transcript.transcript_id()] = transcript.get_coordinates();
		}
	}
	for (auto coordinates_it = coordinates.begin(); coordinates_it != coordinates.end(); ++coordinates_it) {
		if (coordinates_it->second != Transcript::UNDEFINED) {
			std::cout << coordinates_it->first << '\t' << coordinates_it->second.first << '\t' << coordinates_it->second.second << '\n';
		}
	}
	std::cout.flush();
}
