// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2021-08-02
// Released under Apache License 2.0

#include <set>
#include <map>
#include <sstream>
#include <iostream>
#include <fstream>

/// <summary>
/// Add new chromosome to a dictionary
/// </summary>
/// <param name="dictionary">Dictionary of chromosomes</param>
/// <param name="name">Name of the chromosome</param>
/// <param name="sequence">Sequence of the new chromosome</param>
/// <returns>TRUE if no error occured</returns>
bool add_sequence(std::map<std::string, std::string>& dictionary, const std::string& name, std::stringstream& sequence) {
  if (!name.empty() || sequence.tellp() != 0) { // If it is not first (or empty) chromosome, add it to the dictionary
	if (dictionary.emplace(name, sequence.str()).second) {
	  sequence.str("");
	  sequence.clear();
	} else {
	  std::cerr << "Multiple sequences with the same id '" << name << "'." << std::endl;
	  return false;
	}
  }
  return true;
}

/// <summary>
/// Returns next element.
/// It takes advantage of knowledge that last element of a line should not be parsed by this method for used input files.
/// </summary>
/// <param name="pos">Start position of the current element</param>
/// <param name="line">Elements to be separated</param>
/// <param name="separator">Separator for splitting the line</param>
/// <param name="element">The current element (output)</param>
/// <param name="error">Error message to be printed if not enough elements</param>
/// <returns>TRUE if no error occured</returns>
bool get_element(size_t &pos, const std::string &line, const char separator, std::string &element, const std::string &error) {
  size_t next = line.find(separator, pos);
  if (next == line.npos) {
	std::cerr << "" << error << ": '" << line << "'." << std::endl;
	return false;
  }
  element = line.substr(pos, next - pos);
  pos = next + 1;
  return true;
}

/// <summary>
/// Returns next tab-separated element.
/// It takes advantage of knowledge that last element of a line should not be parsed by this method for used input files.
/// </summary>
/// <param name="pos">Start position of the current element</param>
/// <param name="line">Tab-separated elements</param>
/// <param name="element">The current element (output)</param>
/// <returns>TRUE if no error occured</returns>
bool get_element(size_t& pos, const std::string& line, std::string& element) {
  return get_element(pos, line, '\t', element, "Not enough columns in a line within annotations file");
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
	std::cout << "gc_content <genome> <annotations>\t Compute GC content for each feature type and gene id\n";
	std::cout << "                                 \t based on <genome> in FASTA format and\n";
	std::cout << "                                 \t its <annotations> in GTF file format." << std::endl;
	return (argc == 1) ? 0 : 1;
  }

  // Dictionary <chromosome name; sequence>
  std::map<std::string, std::string> sequences;
  { // Processing input FASTA file
	// Name of the current chromosome
	std::string header;
	// Sequence of the current chromosome
	std::stringstream buffer;
	// Input FASTA file
	std::ifstream sequences_input(argv[1]);
	for (std::string line; std::getline(sequences_input, line); ) {
	  if (line.empty()) {
		std::cerr << "Unexpected empty line within sequences file '" << argv[1] << "'." << std::endl;
		return 2;
	  }
	  if (line[0] == '>') { // Header line
		if (!add_sequence(sequences, header, buffer)) return 3;
		// Trim '>' from header line and take first space-separated part (Ensembl format)
		size_t sep = line.find(' ', 1);
		if (sep == line.npos) {
		  header = line.substr(1);
		} else {
		  header = line.substr(1, sep - 1);
		}
	  } else {
		buffer << line;
	  }
	}
	// Finalization
	if (!add_sequence(sequences, header, buffer)) return 3;
	sequences_input.close();
  }
  
  // UTR5, CDS, etc.
  std::set<std::string> features;
  // Chromosome => gene => feature type => base => count
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<char, size_t>>>> stats;
  { // Processing input GTF file
	// Input GTF file
	std::ifstream annotations_input(argv[2]);
	for (std::string line; std::getline(annotations_input, line); ) {
	  if (line.empty()) {
		std::cerr << "Unexpected empty line within annotations file '" << argv[2] << "'." << std::endl;
		return 4;
	  }
	  if (line[0] != '#') { // It is not a comment
		size_t position = 0;
		std::string element;
		// Parse chromosome
		if (!get_element(position, line, element)) return 5;
		std::string chromosome = element;
		// Parse source
		if (!get_element(position, line, element)) return 5;
		// Parse feature type
		if (!get_element(position, line, element)) return 5;
		std::string feature = element;
		if (feature == "gene" || feature == "transcript") continue;
		features.emplace(feature);
		// Parse start 1-based index
		if (!get_element(position, line, element)) return 5;
		size_t from = std::stoull(element);
		// Parse end 1-based index
		if (!get_element(position, line, element)) return 5;
		size_t to = std::stoull(element);
		// Parse score
		if (!get_element(position, line, element)) return 5;
		// Parse strand
		if (!get_element(position, line, element)) return 5;
		if (element.size() != 1 || (element != "+" && element != "-")) {
		  std::cerr << "Unexpected strand format in a line within annotations file: '" << line << "'." << std::endl;
		  return 34;
		}
		bool strand = element[0] == '+';
		// Parse phase
		if (!get_element(position, line, element)) return 5;

		// Parse gene_id
		position = line.find("gene_id \"", position);
		if (position == line.npos) {
		  std::cerr << "Missing 'gene_id' field in a line within annotations file: '" << line << "'." << std::endl;
		  return 8;
		}
		position += 9;
		if (!get_element(position, line, '"', element, "Unenclosed 'gene_id' field in a line within annotations file")) return 13;
		std::string gene = element;

		// Stats for the current gene and feature
		std::map<char, size_t>& stat = stats[chromosome][gene][feature];
		// Sequence of the chromosome where the current gene is
		std::string& sequence = sequences[chromosome];
		for (size_t i = from - 1; i < to; ++i) {
		  char base = sequence[i];
		  if (!strand) { // Complementary base if the gene is in a reverse strand
			switch (sequence[i]) {
			  case 'A':
				base = 'T';
				break;
			  case 'C':
				base = 'G';
				break;
			  case 'G':
				base = 'C';
				break;
			  case 'T':
			  case 'U':
				base = 'A';
				break;
			  case 'N':
				break;
			  default:
				std::cerr << "Unsuported base code: '" << sequence[i] << "'." << std::endl;
				return 30;
				break;
			}
		  }
		  ++stat[base];
		}
	  }
	}
  }
  
  // Print header
  std::cout << "gene_id";
  for (auto headers_it = features.begin(); headers_it != features.end(); ++headers_it) {
	std::cout << '\t' << *headers_it;
  }
  // Print stats
  for (auto stats_it = stats.begin(); stats_it != stats.end(); ++stats_it) { // For each chromosome
	for (auto stats_it_it = stats_it->second.begin(); stats_it_it != stats_it->second.end(); ++stats_it_it) { // For each gene
	  std::cout << '\n' << stats_it_it->first;
	  for (auto headers_it = features.begin(); headers_it != features.end(); ++headers_it) { // For each feature type
		auto substats = stats_it_it->second.find(*headers_it);
		if (substats == stats_it_it->second.end()) { // The current gene has no region of the current feature type
		  std::cout << "\tNA";
		} else {
		  size_t gc = substats->second['C'] + substats->second['G'];       // Returns default (0) if not exists, so check of existency is not necessary.
		  size_t all = gc + substats->second['A'] + substats->second['T']; // Ns are ignored for the stats.
		  if (substats->second.find('U') != substats->second.end()) {      // However default have some overheads, so check is recommended if probability of existency is low.
			all += substats->second['U'];
		  }
		  std::cout << '\t' << 1.0 * gc / all;
		}
	  }
	}
  }
  std::cout << std::endl;
  std::cout.flush();
  
  return 0;
}

