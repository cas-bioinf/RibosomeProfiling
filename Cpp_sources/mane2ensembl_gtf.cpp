// Created by Jan Jelínek (jan.jelinek@biomed.cas.cz)
// Last update: 2021-08-26
// Released under Apache License 2.0

#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>

/// <summary>
/// Determine, whether it is 5'UTR, or 3'UTR
/// </summary>
/// <param name="five_lower">UTR end for forward strand; start codon end for reverse strand</param>
/// <param name="five_upper">Start codon begin for forward strand; UTR begin for reverse strand</param>
/// <param name="three_lower">Stop codon begin for forward strand; UTR end for reverse strand</param>
/// <param name="three_upper">UTR begin for forward strand; stop codon end for reverse strand</param>
/// <param name="parts">Splitted line</param>
/// <param name="transcript">The current transcript</param>
void classify_utr(const size_t five_lower, const size_t five_upper, const size_t three_lower, const size_t three_upper, std::vector<std::string>& parts, const std::string& transcript) {
  if (five_lower < five_upper) {
    parts[2] = "five_prime_utr";
  } else if (three_lower <= three_upper) {
    parts[2] = "three_prime_utr";
  } else {
    std::cerr << "Unexpected file format - UTR region occures between start and stop codons for " << transcript << std::endl;
  }
}

/// <summary>
/// Remove stop codon from 3'UTR region
/// </summary>
/// <param name="parts">The current line</param>
/// <param name="trimmed">What length of 3'UTR was stripped of stop codon</param>
/// <param name="range">Boundaries of the line</param>
/// <param name="index">'0' for forward strand, '1' for reverse strand</param>
/// <param name="lower">Boundaries of the current line for forward strand; stop codon boundaries for reverse strand</param>
/// <param name="upper">Stop codon boundaries for forward strand; boundaries of the current line for reverse strand</param>
/// <param name="new_boundary">New boundary of the current line if stop codon is removed</param>
/// <param name="transcript">The current transcript</param>
/// <returns>TRUE if line should be further processed; FALSE if it is completely covered by a stop codon</returns>
bool adjust_boundaries(std::vector<std::string>& parts, size_t& trimmed, size_t range[2], 
                       const size_t index, const size_t lower[2], const size_t upper[2], const size_t new_boundary, const std::string& transcript) {
  if (lower[0] <= upper[1]) { // The current exon is not all after stop codon
    if (trimmed >= 3) { // Inconsistency in exons 
      std::cerr << transcript << " contains stop_codon longer than 3 bases";
    }
    trimmed += range[1] - range[0] + 1;
    if (lower[1-index] <= upper[1-index]) { // The line is just a duplication of a stop_codon, can be omitted
      return false;
    }
    // Update boundaries of the 3'UTR
    range[index] = new_boundary;
    parts[index+3] = std::to_string(range[index]);
  }
  return true;
}

/// <summary>
/// Make final check for a previous transcript and initialize values for a new gene/ transcript
/// </summary>
/// <param name="transcript">The previous transcript</param>
/// <param name="value">What value should be used for initialization of lengths and trimmed</param>
/// <param name="start_codon">Start codon boundaries</param>
/// <param name="start_codon_length">Start codon length</param>
/// <param name="stop_codon">Stop codon boundaries</param>
/// <param name="stop_codon_length">Stop codon length</param>
/// <param name="trimmed">Length of 3'UTR regions that were stripped of stop codon</param>
void init(const std::string& transcript, const size_t value, size_t start_codon[2], size_t& start_codon_length, size_t stop_codon[2], size_t& stop_codon_length, size_t& trimmed) {
  if (start_codon_length < 3) {
    std::cerr << transcript << " does not contains a complete start_codon." << std::endl;
  }
  if (stop_codon_length < 3) {
    std::cerr << transcript << " does not contains a complete stop_codon." << std::endl;
  }
  if (trimmed < 3) {
    std::cerr << transcript << " does not have whole stop_codon duplicated in UTR region." << std::endl;
  }
  start_codon[0] = -1;
  start_codon[1] = -1;
  start_codon_length = value;
  stop_codon[0] = -1;
  stop_codon[1] = -1;
  stop_codon_length = value;
  trimmed = value;
}

/// <summary>
/// Process start/stop codon and update corresponding variables
/// </summary>
/// <param name="boundaries">Boundaries of the codon</param>
/// <param name="length">Length of the codon</param>
/// <param name="line">The current line</param>
/// <param name="parts">Splitted current line</param>
/// <returns>TRUE if line should be further processed; FALSE if it should be skipped</returns>
bool process_codon(size_t boundaries[2], size_t& length, const std::string& line, const std::vector<std::string>& parts, const std::string &codon) {
  if (length == 3) { // It is already complete
    std::cerr << "Unexpected file format - multiple " << codon << " codons: " << line << std::endl;
    return false;
  }
  // Update boundaries
  size_t range[] = { std::stoull(parts[3]), std::stoull(parts[4]) };
  if (length == 0) {
    boundaries[0] = range[0];
    boundaries[1] = range[1];
  } else { // To avoid distinguishing forward and reverse strand
    boundaries[0] = std::min(boundaries[0], range[0]);
    boundaries[1] = std::max(boundaries[1], range[1]);
  }
  // Update length
  length += range[1] - range[0] + 1;
  if (length > 3) {
    std::cerr << "Unexpected file format - strange " << codon << "_codon length: " << line << std::endl;
  }
  return true;
}

int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "mane2ensembl_gtf <input> <output>\tTakes MANE's annotations in GTF format for Ensembl identifiers from file <input>,\n";
    std::cout << "                                 \ttransform them to be consistent with annotations in GTF format provided by Ensembl\n";
    std::cout << "                                 \tand store them in file <output>.\n\n";
    std::cout << "Transformations are:\n";
    std::cout << "1. 'chr' is removed from beginning of seqname;\n";
    std::cout << "2. 'UTR' feature is classified as 'five_prime_utr' or 'three_prime_utr';\n";
    std::cout << "3. stop_codon is not considered to be a part of 3'UTR;\n";
    std::cout << "4. gene_id attribute is splitted into gene_id and gene_version, the same for transcript_id etc.;\n";
    std::cout << "5. 'gene_type' tag is replaced by 'gene_biotype', the same for 'transcript_type'." << std::endl;
    return (argc == 1) ? 0 : 1;
  }

  // Input GTF file
  std::ifstream input(/*/"C:\\Working\\Valasek\\_references\\genome\\ensembl\\MANE.GRCh38.v0.95.select_ensembl_genomic.gtf" /*/argv[1]/**/);
  // Output GTF file
  std::ofstream output(/*/"C:\\Working\\Valasek\\_references\\genome\\ensembl\\MANE.GRCh38.v0.95.select_ensembl_genomic-104.gtf" /*/argv[2]/**/);
  // Transcript_id of the currently processed transcript
  std::string current_transcript = "";
  // Boundaries of start codon
  size_t start_codon[2] = { (size_t)-1, (size_t)-1 };
  // Length of start codon
  size_t start_codon_length = 3;
  // Boundaries of stop codon
  size_t stop_codon[2] = { (size_t)-1, (size_t)-1 };
  // Length of stop codon
  size_t stop_codon_length = 3;
  // Length of 3'UTR regions that were stripped of stop codon
  size_t trimmed = 3;
  for (std::string line; std::getline(input, line); ) {
    if (line.empty()) { // Not expected
      std::cerr << "Unexpected empty line." << std::endl;
    } else if (line[0] == '#') { // Comments
      output << line << '\n';
    } else { // The interesting part
      // Parsing columns
      size_t from = 0;
      
      // Parsing first 8 features fields
      std::vector<std::string> parts;
      parts.reserve(8);
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
      // The ninth attributes field
      std::string attributes = line.substr(from);

      // Solution of problem #1 - trim 'chr' from beginning of seqid
      if (parts[0].rfind("chr", 0) != 0) {
        std::cerr << "Unexpected line format - seqname does not start with 'chr': " << line << std::endl;
      }
      parts[0] = parts[0].substr(3);

      if (parts[2] == "gene") {
        // Initialize and set values that do not endanger transcript checks
        init(current_transcript, 3, start_codon, start_codon_length, stop_codon, stop_codon_length, trimmed);
        current_transcript = "";
      } else if (parts[2] == "transcript") {
        // Initialize
        init(current_transcript, 0, start_codon, start_codon_length, stop_codon, stop_codon_length, trimmed);
        size_t from = attributes.find("transcript_id \"");
        if (from == attributes.npos) {
          std::cerr << "Unexpected line format - transcript row does not contains attribute transcript_id: " << line << std::endl;
        }
        from += 15;
        size_t to = attributes.find('"', from);
        if (to == attributes.npos) {
          std::cerr << "Unexpected line format - transcript row does not contains unfinished transcript_id: " << line << std::endl;
        }
        current_transcript = attributes.substr(from, to-from);
      } else if (parts[2] == "start_codon") {
        if (!process_codon(start_codon, start_codon_length, line, parts, "start")) continue;
      } else if (parts[2] == "stop_codon") {
        if (!process_codon(stop_codon, stop_codon_length, line, parts, "stop")) continue;
      } else if (parts[2] == "UTR") {
        if (start_codon[0] == -1 || start_codon[1] == -1) {
          std::cerr << "Unexpected file format - start_codon line is missing or is not prior an UTR line: " << line << std::endl;
          continue;
        }
        if (stop_codon[0] == -1 || stop_codon[1] == -1) {
          std::cerr << "Unexpected file format - stop_codon line is missing or is not prior an UTR line: " << line << std::endl;
          continue;
        }
        size_t range[] = { std::stoull(parts[3]), std::stoull(parts[4]) };

        // Solution of problem #2 - classify UTR as 5'UTR, or 3'UTR
        if (parts[6] == "+") {
          classify_utr(range[1], start_codon[0], stop_codon[0], range[0], parts, current_transcript);
        } else if (parts[6] == "-") {
          classify_utr(start_codon[1], range[0], range[1], stop_codon[1], parts, current_transcript);
        } else {
          std::cerr << "Unexpected line format - unsupported strain identifier: " << line << std::endl;
        }

        // Solution of problem #3 - exclude stop codon from 3'UTR
        if (parts[2] == "three_prime_utr") {
          if (parts[6] == "+") {
            if (!adjust_boundaries(parts, trimmed, range, 0, range, stop_codon, stop_codon[1]+1, current_transcript)) continue;
          } else if (parts[6] == "-") {
            if (!adjust_boundaries(parts, trimmed, range, 1, stop_codon, range, stop_codon[0]-1, current_transcript)) continue;
          } else {
            std::cerr << "Unexpected line format - unsupported strain identifier: " << line << std::endl;
          }
        }
      }

      // Solution of problem #4 - split '*_id "id.version"' into '*_id "id"; *_version "version"'
      for (size_t i = attributes.find("_id \""); i != attributes.npos; i = attributes.find("_id \"", i+5)) {
        size_t begin = attributes.rfind(' ', i);
        if (begin == attributes.npos) {
          begin = 0;
        }
        size_t dot = attributes.find('.', i + 5);
        size_t quote = attributes.find('"', i + 5);
        if (quote == attributes.npos) {
          std::cerr << "Uncompleted value of '" << attributes.substr(begin, i + 3 - begin) << "': " << line << std::endl;
        }
        if (dot == attributes.npos || quote < dot) {
          std::cerr << "Unexpected format of '" << attributes.substr(begin, i + 3 - begin) << "': " << line << std::endl;
        }
        attributes.replace(dot, 1, "\"; " + attributes.substr(begin, i - begin) + "_version \"");
      }

      // Solution of problem #5 - replace '*_type' with '*_biotype'
      for (size_t i = attributes.find("_type \""); i != attributes.npos; i = attributes.find("_type \"", i + 10)) {
        attributes.insert(i + 1, "bio");
      }

      // Print repaired line
      for (size_t i = 0; i < parts.size(); i++) {
        output << parts[i] << '\t';
      }
      output << attributes << '\n';
    }
  }
  output.flush();
  output.close();
  input.close();

  return 0;
}
