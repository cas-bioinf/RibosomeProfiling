#!/bin/bash

########################################################################################################################
# A script to download reference data from source databases.                                                           #
#                                                                                                                      #
# Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-07-29; license: Apache License 2.0             #
########################################################################################################################

release_ensembl=115
release_mane=1.4
if [ $# -eq 3 ]; then
  release_ensembl="$2"
  release_mane="$3"
elif [ $# -ne 1 ]; then
  echo "A script to download reference data from source databases.";
  echo;
  echo "1-reference.sh <location> [<ensembl_version> <MANE_version>]	 Downloads reference data and save them to <location>. Ensembl is downloaded in version no. <ensembl_version> (${release_ensembl} by default), and MANE is downloaded in version no. <MANE_version> (${release_mane} by default).";
  echo;
  echo "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-07-29; license: Apache License 2.0";
  exit;
fi

# Consistency with the variable name from '0-variables.sh' so that commands work even when using copy-paste
references="$( realpath "$1" )/"

#### Genome
mkdir -p ${references}genome/ensembl
# Ensembl
wget -P ${references}genome/ensembl/ "ftp://ftp.ensembl.org/pub/release-${release_ensembl}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
wget -P ${references}genome/ensembl/ "ftp://ftp.ensembl.org/pub/release-${release_ensembl}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release_ensembl}.gtf.gz"
gunzip ${references}genome/ensembl/{Homo_sapiens.GRCh38.dna.primary_assembly.fa,Homo_sapiens.GRCh38.${release_ensembl}.gtf}.gz &
# MANE
wget -P ${references}genome/ensembl/ "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${release_mane}/MANE.GRCh38.v${release_mane}.ensembl_genomic.gtf.gz"
wget -P ${references}genome/RefSeq/ "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_${release_mane}/MANE.GRCh38.v${release_mane}.summary.txt.gz"
gunzip "${references}genome/ensembl/MANE.GRCh38.v${release_mane}.ensembl_genomic.gtf.gz" &
gunzip "${references}genome/RefSeq/MANE.GRCh38.v${release_mane}.summary.txt.gz" &
# uORFdb
wget -P ${references}genome/uORFdb/ https://www.bioinformatics.uni-muenster.de/tools/uorfdb/download/uORF_dump_uORFdb.tsv

#### rRNA
mkdir -p ${references}rRNA
# NCBI
wget -O- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_023363.1&rettype=fasta' > ${references}rRNA/ncbi.fa
wget -O- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003285.3&rettype=fasta' >>${references}rRNA/ncbi.fa
wget -O- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003286.4&rettype=fasta' >>${references}rRNA/ncbi.fa
wget -O- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NR_003287.4&rettype=fasta' >>${references}rRNA/ncbi.fa
# Ensembl
wget -O- 'https://rest.ensembl.org/sequence/id/ENSG00000211459?content-type=text/x-fasta' > ${references}rRNA/ensembl.fa
wget -O- 'https://rest.ensembl.org/sequence/id/ENSG00000210082?content-type=text/x-fasta' >>${references}rRNA/ensembl.fa

#### tRNA
mkdir -p ${references}tRNA
# GtRNAdb; the old web is currently still valid and provide the same content so I'm keeping it here for backup.
wget -P ${references}tRNA/ https://gtrnadb.org/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa
# wget -P ${references}tRNA/ https://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa
# tRNAdb
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
echo "Reference tRNA from tRNAdb must be downloaded manually."
echo "Go to https://tdb.bioinf.uni-leipzig.de/search"
echo "Select/ keep selected 'tRNA Genes' database."
echo "Select 'Homo sapiens' organism."
echo "Click on 'Search'."
echo "Click on '+ Add All'."
echo "Click on 'Add Selection to Downloads' and confirm."
echo "Click on the download icon (to the right of View Downloads)."
echo "Select 'FASTA' file type."
echo "Select/ keep selected 'Sequence' type."
echo "Select/ keep selected 'ID' or 'tRNA_ID' on the first position of Identifier Format."
echo "Click on 'Confirm Download' and save the resulting file in ${references}."
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
echo;
