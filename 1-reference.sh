#!/bin/bash

if [ $# -eq 1 ]; then
  release_ensembl=104
elif [ $# -ne 2 ]; then
  echo "A script to download reference data from source databases.";
  echo;
  echo "reference-local <location> [<ensembl_version>]	 Downloads reference data and save them to <location>. In the case of ensembl, version no. <ensembl_version> is used (104 by default).";
  exit;
else
  release_ensembl="$2"
fi

references="$( realpath "$1" )/"

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
# tRNAdb
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
echo "Reference tRNA from tRNAdb must be downloaded manually."
echo "Go to http://trna.bioinf.uni-leipzig.de/DataOutput/Search?vOrg=Homo+sapiens&vTax=9606"
echo "keep the form prefilled and search in tRNA genes database."
echo "Check 'select all sequences of search', select 'Download FASTA converted to DNA' and submit."
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
echo;
# GtRNAdb
wget -P ${references}tRNA/ http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.fa

#### Genome
mkdir -p ${references}genome/ensembl
# Ensembl
wget -P ${references}genome/ensembl/ "ftp://ftp.ensembl.org/pub/release-${release_ensembl}/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
wget -P ${references}genome/ensembl/ "ftp://ftp.ensembl.org/pub/release-${release_ensembl}/gtf/homo_sapiens/Homo_sapiens.GRCh38.${release_ensembl}.gtf.gz"
gunzip ${references}genome/ensembl/{Homo_sapiens.GRCh38.dna.primary_assembly.fa,Homo_sapiens.GRCh38.${release_ensembl}.gtf}.gz
