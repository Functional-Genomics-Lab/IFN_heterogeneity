#! /bin/bash
# load STAR
module load star/2.7.3a

# Read file names
while IFS=',' read -ra array; do
  #ar1+=("${array[0]}")
  ar2+=("${array[1]}")
  #ar3+=("${array[2]}")
  #ar4+=("${array[3]}")
done < ~/project/Gozde_data/RNA-seq/raw_data/Fastq\ Report_ds.ff48cf3bb65844319cff97f8b85b628a/Reports/Quality_Metrics.csv                                    
IFS=" " read -r -a ar2 <<< "$(tr ' ' '\n' <<< "${ar2[@]:1:48}" | sort -u | tr '\n' ' ')"
#printf '%s\n' "${ar2[@]}"

# read the file locations
STAR_DIR="~/project/Gozde_data/humanhg38plusSeV/star"
FASTQ_DIR="~/project/Gozde_data/RNA-seq/raw_data"

# Run it for all the samples of interest:
for SAMPLE in "${ar2[@]}"
do
# Define the list of fastq files per sample
FILES=`'ls' ${FASTQ_DIR}/${SAMPLE}*/*.fastq.gz | paste -s -d ' ' -`

#printf '%s\n' $SAMPLE
#printf '%s\n' $FILES
# Run STAR
STAR --runThreadN 20 \
--genomeDir ${STAR_DIR} --readFilesIn $FILES \
--readFilesCommand gunzip -c --outFileNamePrefix ~/project/Gozde_data/RNA-seq/aligned_reads/${SAMPLE}_ \
--outFilterMultimapNmax 1 \
--outSAMtype BAM SortedByCoordinate \
--alignIntronMin 1 --alignIntronMax 3000
done
