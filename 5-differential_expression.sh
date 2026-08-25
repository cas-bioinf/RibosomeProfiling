#!/bin/bash

########################################################################################################################
# A script to perform differential expression analysis.                                                                #
#                                                                                                                      #
# Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-25; license: Apache License 2.0             #
########################################################################################################################

help() {
  echo "A script to perform differential expression analysis.";
  echo;
  echo "5-differential_expression.sh -h                 	 Prints this help.";
  echo "5-differential_expression.sh <programs> <output>	 Using custom scripts, performs differential expression analysis.";
  echo;
  echo "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-25; license: Apache License 2.0";
}

if [ $# -eq 2 ]; then
  # Consistency with the variable names from '0-variables.sh' so that commands work even if using copy-paste
  programs="$( realpath "$1" )/"
  output="$( realpath "$2" )/"
else
  if [ $# -eq 1 ] && ! [[ "$1" =~ '^-?-h' ]]; then
    >&2 echo "Unexpected argument: '$1'";
    e=1
  elif [ $# -gt 1 ]; then
    >&2 echo "Missing $( echo "4-$#" | bc) argument$( [ $# -lt 3 ] && echo s )";
    e=1
  else
    e=0
  fi
  help;
  exit "$e";
fi

shopt -s extglob

# Create table with gene counts
mkdir -p ${output}_deseq/
table=${output}_deseq/r01M1_sE-Ensembl-CDS.tsv;
i=0;
>${table}.$i;
for file in ${output}*_R1*/!(*[^-_]){FP,mRNA}[-_]*_CDS.tsv; do
  j=$(($i+1));
  join -t$'\t' -a1 -a2 -e0 -o auto ${table}.$i <( echo -ne 'gene\t'; basename ${file%_R1[_-]*} | awk -F'[-_]' 'BEGIN{OFS="_"} $0~/^S25-075/{ print $3, "HEK", substr($4,1,2)=="WT" ? $3== "FP" ? "WT" : "wt" : "Triple", substr($4,length($4),1); next} $0~/^(c|nt)_[1-4]_mRNA$/{ print "mRNA", $1 == "nt" ? $1 : "eIF3cKD", $2; next} $0~/^mRNA_nt_1$/{ print "mRNA_nt_5"; next} $0~/^KP3_mRNA_nt_2$/{ print "mRNA_nt_6"; next} $0~/^KP3_FP_c_4$/{ print "FP_eIF3cKD_1"; next} $0~/^FP_c[5-7]$/{ print "FP_eIF3cKD", substr($2,2)-3; next} $0~/^FP_nt_1$/{ print "FP_nt_1"; next} $0~/^KP[23]_FP_nt_[24]$/{ print "FP_nt", $4/2+1; next} $0~/^FP_nt_5[5-7]$/{ print "FP_nt", $3-51; next} {print "Unexpected line: " $0 > "/dev/stderr"}'; head -n -5 $file ) >${table}.$j;
  rm ${table}.$i;
  i=$j;
done;
mv ${table}.$i ${table};


# Create design tables
head -n1 "$table" | cut -f2- | tr '\t' '\n' | sort | awk -F'_' 'BEGIN{OFS="\t"; print "","assay","siRNA","id"} $2=="HEK"{print $0,$1,$3=="wt"?toupper($3):$3, $4; next} {print $0,$1,$2,$3}' >${output}_deseq/design-all.tbl
grep -v '_HEK_' ${output}_deseq/design-all.tbl >${output}_deseq/design-eIF3cKD.tbl
sed '1p;/_HEK_/!d' ${output}_deseq/design-all.tbl >${output}_deseq/design-eIF3cKI.tbl
rm ${output}_deseq/design-all.tbl

# Differential expression etc.
for design in ${output}_deseq/design-*.tbl; do
  dir=${table%.tsv}-$( basename -s.tbl "$design" | cut -d- -f2- )/;
  echo $table $dir;
  mkdir -p "$dir";
  ${programs}RibosomeProfiling/R_scripts/Differential_analysis.R $design $table --output "$dir" --reference2 $( basename -s.tbl ${design#*-} | sed -e 's/^eIF3cKD$/nt/; s/^eIF3cKI$/WT/' ) --pca_ids T --use_bm_cache F &>"${dir}log.txt";
done;

# Add name to the first column (it is easier than repair it within R)
sed -i -e '1s/^""/"gene_id"/' ${output}_deseq/*/*.tbl

# Aggregation of files
for dir in ${output}_deseq/*/; do
  i=1;
  cut -f1,2 $( ls ${dir}!(sig*).tbl | head -n1 ) > ${dir}all.tbl.$i;
  for file in ${dir}!(sig*).tbl; do
    j=$(($i+1));
    join -t $'\t' --header ${dir}all.tbl.$i <( cut -f1,3-7 $file | sed "1 s/\t\"/\t\"$( basename -s .tbl $file )./g" ) >${dir}all.tbl.$j;
    rm ${dir}all.tbl.$i;
    i=$j;
  done;
  join -t $'\t' --header ${dir}all.tbl.$i <( cut -f1,8- $file ) > ${dir}all.tbl;
  rm ${dir}all.tbl.$i;
done;

# translation only
for file in ${output}_deseq/*/sig0.05-TE-!(*[-_]*)_vs_!(*[-_]*).tbl; do
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]++;next} FNR==1{print;next} !a[$1]' ${file/_TE_/_mRNA_} $file >${file%.tbl}-translation_only.tbl;
done;
