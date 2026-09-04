#!/bin/bash

########################################################################################################################
# A script to install programs used during analysis and their prerequisities.                                          #
#                                                                                                                      #
# Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-09-04; license: Apache License 2.0             #
########################################################################################################################

release_fastqc=0.12.1
if [ $# -eq 2 ]; then
  release_fastqc="$2"
elif [ $# -ne 1 ]; then
  echo "A script to install programs used during analysis and their prerequisities.";
  echo;
  echo "2-install.sh <location> [<fastqc_version>]     Install used programs and their prerequisities; and save them to <location> if necessary. In the case of FastQC, version no. <fastqc_version> is used (${release_fastqc} by default).";
  echo;
  echo "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-09-04; license: Apache License 2.0";
  exit;
fi

# Consistency with the variable name from '0-variables.sh' so that commands work even if using copy-paste
programs="$( realpath "$1" )/"

# prerequisities
sudo apt-get install autoconf default-jre g++ gawk git libbz2-dev libcurl4-openssl-dev liblzma-dev libncurses5-dev libxml2-dev make pipx python3-dev r-base xxd
# parallel python3-pip libssl-dev

mkdir -p "$programs"

# FastQC
echo "#### FastQC ####"
# sudo apt-get install fastqc # manual installation is better as it provides a newer version
file=fastqc_v${release_fastqc}.zip
wget -P "${programs}" "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/$file"
unzip "${programs}$file" -d "${programs}"
rm "${programs}$file"
chmod 755 "${programs}FastQC/fastqc"

# MultiQC
echo "#### MultiQC ####"
pipx install multiqc

# Cutadapt
echo "#### Cutadapt ####"
# sudo apt-get install cutadapt # pip installation is better as it provides a significantly newer version
pipx install cutadapt

# Bowtie2
echo "#### Bowtie ####"
# sudo apt-get install bowtie # manual installation is better as it provide a newer version
git clone https://github.com/BenLangmead/bowtie2.git "${programs}bowtie2/"
pushd "${programs}bowtie2/"
make
sudo make install
popd

# STAR
echo "#### STAR ####"
git clone https://github.com/alexdobin/STAR.git "${programs}STAR/"
pushd "${programs}STAR/source/"
make STAR -j1 # Bypassing a reported bug https://github.com/alexdobin/STAR/issues/2672 till the corrective pull request will be accepted
popd

# SAMtools
echo "#### SAMtools ####"
# sudo apt-get install samtools # manual installation is better as it provide a newer version
git clone https://github.com/samtools/samtools.git "${programs}samtools/"
git clone https://github.com/samtools/htslib.git "${programs}htslib/"
pushd "${programs}htslib/"
git submodule update --init --recursive
cd "${programs}samtools/"
autoheader
autoconf -Wno-syntax
./configure
make
sudo make install
popd

# HTSeq
echo "#### HTSeq ####"
pipx install numpy HTSeq

# Custom programs
git clone https://github.com/cas-bioinf/RibosomeProfiling.git "${programs}RibosomeProfiling/"
chmod u+x ${programs}RibosomeProfiling/[1-5]-*.sh ${programs}RibosomeProfiling/R_scripts/Differential_analysis.R
g++ -O3 -o ${programs}mane2ensembl_gtf       "${programs}RibosomeProfiling/Cpp_sources/mane2ensembl_gtf.cpp"
g++ -O3 -o ${programs}filter_reverse_reads   "${programs}RibosomeProfiling/Cpp_sources/filter_reverse_reads.cpp"
g++ -O3 -o ${programs}filter_ambiguous_genes "${programs}RibosomeProfiling/Cpp_sources/filter_ambiguous_genes.cpp"
g++ -O3 -o ${programs}select_transcripts     "${programs}RibosomeProfiling/Cpp_sources/select_transcripts.cpp"
g++ -O3 -o ${programs}select_features        "${programs}RibosomeProfiling/Cpp_sources/select_features.cpp"
g++ -O3 -o ${programs}gc_content             "${programs}RibosomeProfiling/Cpp_sources/gc_content.cpp"

# R libraries
sudo R -q -e 'install.packages(c("BiocManager","dbplyr","dplyr","ggplot2","ggpointdensity","ggrepel","RSQLite","gplots","reshape2","xml2")); BiocManager::install(c("biomaRt","DESeq2","SummarizedExperiment"))'

pipx ensurepath

# Print versions of all installed programs
echo "Installation completed. Installed versions are:";
"${programs}FastQC/fastqc" --version
multiqc --version
echo "cutadapt $( cutadapt --version )"
bowtie2 --version | head -n1
echo "STAR $( "${programs}STAR/source/STAR" --version )"
samtools --version | head -n2 | paste -sd' '
echo "htseq-count $( htseq-count --version )"
echo "mane2ensembl_gtf$(       ${programs}mane2ensembl_gtf       | tail -n1 | cut -d';' -f2 )"
echo "filter_reverse_reads$(   ${programs}filter_reverse_reads   | tail -n1 | cut -d';' -f2 )"
echo "filter_ambiguous_genes$( ${programs}filter_ambiguous_genes | tail -n1 | cut -d';' -f2 )"
echo "select_transcripts$(     ${programs}select_transcripts     | tail -n1 | cut -d';' -f2 )"
echo "select_features$(        ${programs}select_features        | tail -n1 | cut -d';' -f2 )"
echo "gc_content$(             ${programs}gc_content             | tail -n1 | cut -d';' -f2 )"
