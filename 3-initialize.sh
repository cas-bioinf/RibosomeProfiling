#!/bin/bash

########################################################################################################################
# A script to initialize databases based on previously downloaded references.                                          #
#                                                                                                                      #
# Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-14; license: Apache License 2.0             #
########################################################################################################################

# Default values (used for the current experiment)
release_ensembl=115
release_mane=1.4
overlap=149
threads=1

help() {
  echo "A script to initialize databases based on previously downloaded references.";
  echo;
  echo "3-initialize.sh -h                                	   Prints this help.";
  echo "3-initialize.sh [OPTIONS...] <programs> <references>	 Using programs installed globally (bowtie2-build and coreutils) and in <programs> (STAR and custom programs), initialize databases Bowtie2 rRNA and DNA, STAR and uORFdb databases in the <references> directory.";
  echo;
  echo "OPTIONS";
  echo    "-e NUM\t Version of Ensembl to be used (${release_ensembl} by default).";
  echo    "-m STR\t Version of MANE to be used (${release_mane} by default).";
  echo -e "-o NUM\t The value for STAR's '--sjdbOverhang' argument. It should be -e \e[1mmax(ReadLength)-1\e[0m.";
  echo    "-t NUM\t The number of threads to be used to generate the databases.";
  echo;
  echo "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-14; license: Apache License 2.0";
}

while [ $# -gt 2 ]; do
  case $1 in
    -h* | --h* ) >&2 echo "Unexpected number of arguments: '$1'";
                 help;
                 exit 1;
                 ;;
    -e* | --e* ) shift;
                 release_ensembl="$1";
                 ;;
    -m* | --m* ) shift;
                 release_mane="$1";
                 ;;
    -o* | --o* ) shift;
                 overlap="$1";
                 ;;
    -t* | --t* ) shift;
                 threads="$1";
                 ;;
    * )          >&2 echo "Unknown argument: '$1'";
                 help;
                 exit 1;
  esac
  shift
done
if [ $# -eq 2 ]; then
  # Consistency with the variable names from '0-variables.sh' so that commands work even if using copy-paste
  programs="$( realpath "$1" )/"
  references="$( realpath "$2" )/"
else
  if [ $# -eq 1 ] && ! [[ "$1" =~ '^-?-h' ]]; then
    >&2 echo "Unexpected argument: '$1'";
    e=1
  else
    e=0
  fi
  help;
  exit "$e";
fi

# Bowtie2
for type in rRNA tRNA; do
  mkdir -p ${references}${type}/bowtie2/
  bowtie2-build --threads "$threads" "$( ls ${references}${type}/*.f* | paste -sd',' )" "${references}${type}/bowtie2/${type}";
done

# STAR
mkdir -p ${references}genome/ensembl/STAR;
"${programs}STAR/source/STAR" --runThreadN "$threads" --runMode genomeGenerate --genomeDir "${references}genome/ensembl/STAR/" --outFileNamePrefix "${references}genome/ensembl/STAR/" --genomeFastaFiles "${references}genome/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa" --sjdbGTFfile "${references}genome/ensembl/Homo_sapiens.GRCh38.${release_ensembl}.gtf" --sjdbOverhang "$overlap"

# MANE (just transformation to compatible format with Ensembl)
# gene range is not repaired after Plus_Clinical removal, however it is not currently used
${programs}mane2ensembl_gtf <( grep -v 'tag "MANE_Plus_Clinical";' "${references}genome/ensembl/MANE.GRCh38.v${release_mane}.ensembl_genomic.gtf" ) "${references}genome/ensembl/MANE.GRCh38.v${release_mane}.ensembl_genomic-select-ensembl.gtf"

# uORFdb (just extract Homo sapiens rows, and add ensembl identifiers)
# Records assigned to the same ensembl ID probably could be paired
sed -n '1p;/^Homo sapiens\thg38\t/p' ${references}genome/uORFdb/uORF_dump_uORFdb.tsv >${references}genome/uORFdb/uORF_dump_uORFdb-Homo_sapiens-hg38.tsv
join -t$'\t' -12 -21 --header <( cut -f2,6,8 ${references}/genome/RefSeq/MANE.GRCh38.v${release_mane}.summary.txt | sed -e 's/\.[0-9]*//g' | { sed -u 1q; sort -k2,2; } ) <( awk -F'\t' 'BEGIN{OFS=FS}{id=$9;sub("\\..*","",id);print id,$0}' ${references}genome/uORFdb/uORF_dump_uORFdb-Homo_sapiens-hg38.tsv | { sed -u 1q; sort -t$'\t' -k1,1 -k10,10; } ) >${references}genome/uORFdb/uORF_dump_uORFdb-Homo_sapiens-hg38-MANE.tsv
awk -F$'\t' 'BEGIN{OFS=FS;SUBSEP=FS} NR==1{header=$1 OFS $2 OFS $3;codons["complete"]=0;codons["ATG_TRR"]=0} NR>1{if($24==""){stats[$1,$2,$3][""]=0;next} else if($27=="ATG" && ($28=="TAA" || $28=="TAG" || $28=="TGA")) {stats[$1,$2,$3]["ATG_TRR",$34]+=1};types[$34]=0;stats[$1,$2,$3]["complete",$34]+=1}END{printf header; for (codon in codons){for (type in types){printf OFS codon "-" type}; printf OFS codon "-all"}; print""; for (key in stats) {printf key;for (codon in codons){sum=0;for (type in types){sum+=stats[key][codon,type];printf OFS (stats[key][codon,type]==""?0:stats[key][codon,type])}; printf OFS sum}; print ""}}' ${references}genome/uORFdb/uORF_dump_uORFdb-Homo_sapiens-hg38-MANE.tsv >${references}genome/uORFdb/uORF_dump_uORFdb-Homo_sapiens-hg38-MANE-stats.tsv
