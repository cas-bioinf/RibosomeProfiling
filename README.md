# RibosomeProfiling
Methods and scripts that we use for analyzing RiboSeq data

- `. 0-variables.sh` sets important paths into variables (programs, databases, inputs, outputs and logs).
- `1-references.sh "$references"` downloads almost all necessary databases and prints instructions on how to download the one that cannot be downloaded automatically.
- `2-install.sh "$programs"` installs all used programs. Some commands need sudo so check the script whether you are OK with it. It is recommended to restart the terminal after this step so that `PATH` is correctly set.
- `3-initialize.sh "$programs" "$references"` initializes the downloaded databases for usage by various programs. There is a bug in the current version of sort in Ubuntu [`(uutils coreutils) 0.8.0`] causing an error `join: ...: is not sorted:`. In such case, use `gnusort` instead.
- `4-basic.sh "$programs" "$references" "$input" "$input_big" "$output" "$output_big" "$logs"` performs a basic processing of input files - it aligns reads on a genome and a transcriptome and counts reads per genes.
- `5-differential_expression.sh "$programs" "$output"` performs a differential expression analysis. Due to the current changes in Ensembl, it is recommended to run it with an argument `-s https://sep2025.archive.ensembl.org`.
- `Cpp_sources` contains various C++ programs for computationally demanding tasks.
- `R_scripts` contains various R scripts mostly for visualization and statistical analysis.

Citations:
- Anna Herrmannová, Jan Jelínek, Klára Pospíšilová, Farkas Kerényi, Tomáš Vomastek, Kathleen Watt, Jan Brábek, Mahabub Pasha Mohammad, Susan Wagner, Ivan Topisirovic, Leoš Shivaya Valášek (2024) **Perturbations in eIF3 subunit stoichiometry alter expression of ribosomal proteins and key components of the MAPK signaling pathways** *eLife* **13**:RP95846 (https://doi.org/10.7554/eLife.95846.3)
