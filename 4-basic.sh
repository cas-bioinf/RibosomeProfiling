#!/bin/bash

########################################################################################################################
# A script to make a basic processing of input files.                                                                  #
#                                                                                                                      #
# Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-18; license: Apache License 2.0             #
########################################################################################################################

# Default values (used for the current experiment)
default_references='${references}'
default_ensembl="${default_references}genome/ensembl/Homo_sapiens.GRCh38.115.gtf"
default_mane="${default_references}genome/ensembl/MANE.GRCh38.v1.4.ensembl_genomic-select-ensembl.gtf"
threads=1

help() {
  mnt=$( realpath . | cut -d'/' -f3 );
  if [[ "$mnt" == "" ]]; then
    mnt=c;
  fi
  echo "A script to make a basic processing of input files.";
  echo;
  echo "4-basic.sh -h                                                                                     Prints this help.";
  echo "4-basic.sh [OPTIONS...] <programs> <references> <input> <input_big> <output> <output_big> <logs>  Using programs installed globally (cutadapt, bowtie2, samtools, HTSeq-count and coreutils) and in <programs> (FastQC, STAR and custom programs), and using databases installed in <references> process all '*_R1*.*.gz' files from <input_big>. The files should be in the fastq format. The files will be temporally copied into <input> directory. Output files will be in <output> and intermediate files will be preserved in <output_big>, each file will have its own subdirectory named after its name without extension. Logs will be stored in <logs> directory. This structure was chosen because it allows you to have large files on a slow but large HDD and process them on a small but fast SSD";
  echo;
  echo "OPTIONS";
  echo -e "-e NUM\t Path to the Ensembl annotations in GTF format to be used ('${default_ensembl}' by default).";
  echo -e "-m STR\t Path to the MANE annotations in GTF format to be used ('${default_mane}' by default).";
  echo -e "-t NUM\t The number of threads to be used to generate the databases.";
  echo;
  echo "REMARKS";
  echo "If you got 'Exiting because of FATAL ERROR: could not create FIFO file...' error, you run the script on a disc without named pipes allowed (it is usual in the Windows WSL). In such case, you must enable them first using the following commands (move out of the current disc; unmount the disc; remount it with named pipes allowed; return back):";
  # It is formatted this way to allow copy-paste from this file without additional edit
  echo '--------------------------------------------
pushd /;
sudo umount /mnt/'"$mnt"';
sudo mount -t drvfs '"${mnt^^}"':\\ /mnt/'"$mnt"' -o metadata;
popd;
--------------------------------------------';
  echo "If you wish to set it permanently, you can add the following lines in a file '/etc/wsl.conf' instead:";
  echo '--------------------
[automount]
options = "metadata"
--------------------';
  echo;
  echo "Created by Jan Jelínek (jan.jelinek@biomed.cas.cz); last update: 2026-08-18; license: Apache License 2.0";
}

while [ $# -gt 7 ]; do
  case $1 in
    -h* | --h* ) >&2 echo "Unexpected number of arguments: '$1'";
                 help;
                 exit 1;
                 ;;
    -e* | --e* ) shift;
                 ensembl="$1";
                 ;;
    -m* | --m* ) shift;
                 mane="$1";
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

if [ $# -eq 7 ]; then
  # Consistency with the variable names from '0-variables.sh' so that commands work even if using copy-paste
  programs="$( realpath "$1" )/"
  references="$( realpath "$2" )/"
  input="$( realpath "$3" )/"
  input_big="$( realpath "$4" )/"
  output="$( realpath "$5" )/"
  output_big="$( realpath "$6" )/"
  logs="$( realpath "$7" )/"
else
  if [ $# -eq 1 ] && ! [[ "$1" =~ '^-?-h' ]]; then
    >&2 echo "Unexpected argument: '$1'";
    e=1
  elif [ $# -gt 1 ]; then
    >&2 echo "Missing $( echo "7-$#" | bc) argument$( [ $# -lt 6 ] && echo s )";
    e=1
  else
    e=0
  fi
  help;
  exit "$e";
fi

# Default values (used for the current experiment)
: ${ensembl:=${default_ensembl/$default_references/$references}}
: ${mane:=${default_mane/$default_references/$references}}

mkdir -p "$input" "$logs"

# STAR database can be preloaded and used multiple times to save some time
"${programs}STAR/source/STAR" --genomeDir "${references}genome/ensembl/STAR" --genomeLoad LoadAndExit;
for file in ${input_big}*_R1*.*.gz; do
  # Init
  id=$( basename "${file%%.*}" );
  echo "$id        $file";
  mkdir -p "${output}$id/" "${output_big}$id/";
  input_file="${input}$( basename "$file" )";
  cp "$file" "$input_file";
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "$input_file" >"${logs}FastQC.${id}.out" 2>"${logs}FastQC.${id}.err";

  # Cutadapt
  id_trimmed="${id}-u3a10m15t";
  cutadapt -u 3 -a AAAAAAAAAA -m 15 $( [[ "${id^^}" == *FP* ]] && [[ "$id" != S25-075-* ]] && echo "-M 45" ) -j $threads --trim-n --trimmed-only -o "${output}$id/${id_trimmed}.fastq.gz" "$input_file" >"${logs}cutadapt.${id_trimmed}.out" 2>"${logs}cutadapt.${id_trimmed}.err";
  rm "$input_file" &
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "${output}$id/${id_trimmed}.fastq.gz" >"${logs}FastQC.${id_trimmed}.out" 2>"${logs}FastQC.${id_trimmed}.err";

  # Bowtie2 - rRNA
  id_rRNA="${id_trimmed}-rRNA";
  "bowtie2" -p $threads --un-gz "${output}$id/${id_rRNA}.fastq.gz" -x "${references}rRNA/bowtie2/rRNA" -U "${output}$id/${id_trimmed}.fastq.gz" -S /dev/null >"${logs}bowtie2.${id_rRNA}.out" 2>"${logs}bowtie2.${id_rRNA}.err";
  mv "${output}$id/${id_trimmed}.fastq.gz" "${output_big}$id/" &
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "${output}$id/${id_rRNA}.fastq.gz" >"${logs}FastQC.${id_rRNA}.out" 2>"${logs}FastQC.${id_rRNA}.err";

  # Bowtie2 - tRNA
  id_tRNA="${id_rRNA}-tRNA";
  "bowtie2" -p $threads --un-gz "${output}$id/${id_tRNA}.fastq.gz" -x "${references}tRNA/bowtie2/tRNA" -U "${output}$id/${id_rRNA}.fastq.gz"    -S /dev/null >"${logs}bowtie2.${id_tRNA}.out" 2>"${logs}bowtie2.${id_tRNA}.err";
  mv "${output}$id/${id_rRNA}.fastq.gz" "${output_big}$id/" &
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "${output}$id/${id_tRNA}.fastq.gz" >"${logs}FastQC.${id_tRNA}.out" 2>"${logs}FastQC.${id_tRNA}.err";

  # STAR
  id_aligned_partial="${id_tRNA}-r01M1_sE_";
  id_aligned_sorted="${id_aligned_partial}Aligned.sortedByCoord.out";
  id_aligned_transcriptome="${id_aligned_partial}Aligned.toTranscriptome.out";
  "${programs}STAR/source/STAR" --runThreadN $threads --genomeDir "${references}genome/ensembl/STAR" --genomeLoad LoadAndKeep --readFilesIn "${output}$id/${id_tRNA}.fastq.gz" --readFilesCommand zcat --limitBAMsortRAM  20000000000 --outFileNamePrefix "${output}$id/${id_aligned_partial}" --outSAMtype BAM SortedByCoordinate --outFilterMismatchNoverLmax 0.1 --outFilterMultimapNmax 1 --chimSegmentMin 12 --quantTranscriptomeSAMoutput BanSingleEnd --quantMode TranscriptomeSAM GeneCounts --outSAMattributes All >"${logs}STAR.${id_aligned_partial}.out" 2>"${logs}STAR.${id_aligned_partial}.err";
  mv "${output}$id/${id_tRNA}.fastq.gz" "${output_big}$id/" &
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "${output}$id/${id_aligned_sorted}.bam"        >"${logs}FastQC.${id_aligned_sorted}.out"        2>"${logs}FastQC.${id_aligned_sorted}.err";
  samtools view -F 256 -h "${output}$id/${id_aligned_transcriptome}.bam" | samtools fastq | "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" stdin:"${id_aligned_transcriptome}-1" >"${logs}FastQC.${id_aligned_transcriptome}-1.out" 2>"${logs}FastQC.${id_aligned_transcriptome}-1.err";

  # Discards alignments aligned to a reverse transcript in the alignment to transcriptome;
  # The strange thing with /dev/fd/3 is because in the case >( samtools... ) the bash does not wait for its finish
  id_forward=${id_aligned_transcriptome}-forward;
  { "${programs}filter_reverse_reads" <( samtools view -h "${output}$id/${id_aligned_transcriptome}.bam" ) /dev/fd/3 >"${logs}filter_reverse_reads.${id_forward}.out" 2>"${logs}filter_reverse_reads.${id_forward}.err"; } 3>&1 | samtools view -h -o "${output}$id/${id_forward}.bam" -;
  mv "${output}$id/${id_aligned_transcriptome}.bam" "${output_big}$id/" &
  samtools view -F 256 -h "${output}$id/${id_forward}.bam" | samtools fastq | "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" stdin:"${id_forward}-1" >"${logs}FastQC.${id_forward}-1.out" 2>"${logs}FastQC.${id_forward}-1.err";

  # Discards reads aligned to multiple transcripts
  id_unique=${id_forward}-unique;
  { "${programs}filter_ambiguous_genes" "$ensembl" <( samtools view -h "${output}$id/${id_forward}.bam" ) /dev/fd/3; } 3>&1 >"${logs}filter_ambiguous_genes.${id_unique}.out" 2>"${logs}filter_ambiguous_genes.${id_unique}.err" | samtools view -h -o "${output}$id/${id_unique}.bam" -;
  mv "${output}$id/${id_forward}.bam" "${output_big}$id/"&
  samtools view -F 256 -h "${output}$id/${id_unique}.bam" | samtools fastq | "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" stdin:"${id_unique}-1" >"${logs}FastQC.${id_unique}-1.out" 2>"${logs}FastQC.${id_unique}-1.err";

  # Preserve only alignments to MANE transcripts
  id_mane=${id_unique}-mane;
  { "${programs}select_transcripts" <( grep -o 'transcript_id "[^"]*";' "$mane" | cut -d'"' -f2 | sort -u ) <( samtools view -h "${output_big}$id/${id_unique}.bam" ) /dev/fd/3; } 3>&1 >"${logs}select_transcripts.${id_mane}.out" 2>"${logs}select_transcripts.${id_mane}.err" | samtools view -h -o "${output}$id/${id_mane}.bam" -;
  mv "${output}$id/${id_unique}.bam" "${output_big}$id/"&
  "${programs}FastQC/fastqc" -t $threads -o "${output}$id/" "${output}$id/${id_mane}.bam" >"${logs}FastQC.${id_mane}.out" 2>"${logs}FastQC.${id_mane}.err";

  # HTSeq-count
  samtools index -@ "$threads" "${output}$id/${id_aligned_sorted}.bam";
  id_genes_cds=${id_aligned_sorted}-Ensembl_CDS;
  id_transcripts_cds=${id_mane}-cds;
  htseq-count -f bam -m intersection-strict -t CDS -c "${output}$id/${id_genes_cds}.tsv" "${output}$id/${id_aligned_sorted}.bam" "$ensembl" >"${logs}htseq_count.${id_genes_cds}.out" 2>"${logs}htseq_count.${id_genes_cds}.err" &
  "${programs}select_features" <( samtools view -h "${output}$id/${id_mane}.bam" ) "${mane}" five_prime_utr 1000 mi 5 mi s >( samtools view -h -b -o "${output}$id/${id_mane}-5utr.bam" - ) start_codon 3 le 3 le s >( samtools view -h -b -o "${output}$id/${id_mane}-start.bam" - ) CDS 2 mi 2 mi s >( tee >( samtools view -h -b -o "${output}$id/${id_mane}-cds.bam" - ) | samtools view - | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2,$1}' >"${output}$id/${id_transcripts_cds}.tsv" 2>"${logs}transcripts_counts.${id_transcripts_cds}.err" ) stop_codon 6 le 0 le s >( samtools view -h -b -o "${output}$id/${id_mane}-stop.bam" - ) three_prime_utr 8 mi 1000 mi s >( samtools view -h -b -o "${output}$id/${id_mane}-3utr.bam" - ) >"${logs}/select_features.${id_mane}.out" 2>"${logs}/select_features.${id_mane}.err";
  if [[ ${id,,} == *"mrna"* ]]; then
    id_genes_exon=${id_aligned_sorted}-Ensembl_exon;
    id_transcripts_exon=${id_mane};
    htseq-count -f bam -m intersection-strict -t exon -c "${output}$id/${id_genes_exon}.tsv" "${output}$id/${id_aligned_sorted}.bam" "$ensembl" >"${logs}htseq_count.${id_genes_exon}.out" 2>"${logs}htseq_count.${id_genes_exon}.err" &
    samtools view "${output}$id/${id_transcripts_exon}.bam" | cut -f3 | sort | uniq -c | awk 'BEGIN{OFS="\t"} {print $2,$1}' >"${output}$id/${id_transcripts_exon}.tsv" 2>"${logs}transcripts_counts.${id_transcripts_exon}.err" &
  fi;
done;
"${programs}STAR/source/STAR" --genomeDir "${references}genome/ensembl/STAR" --genomeLoad Remove;
