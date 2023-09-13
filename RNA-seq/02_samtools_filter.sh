#! /bin/bash
# load
module load samtools/gcc/1.10


# Read file names
while IFS=',' read -ra array; do
  #ar1+=("${array[0]}")
  ar2+=("${array[1]}")
  #ar3+=("${array[2]}")
  #ar4+=("${array[3]}")
done < ~/project/Gozde_data/RNA-seq/raw_data/Fastq\ Report_ds.ff48cf3bb65844319cff97f8b85b628a/Reports/Quality_Metrics.csv                                    
IFS=" " read -r -a ar2 <<< "$(tr ' ' '\n' <<< "${ar2[@]:1:48}" | sort -u | tr '\n' ' ')"
#printf '%s\n' "${ar2[@]}"

# read the file location(s)
BAM_DIR=~/project/Gozde_data/RNA-seq/aligned_reads

# Run it for all the samples of interest:
for SAMPLE in "${ar2[@]}"
do
# Define the list of fastq files per sample
FILE=`'ls' ${BAM_DIR}/${SAMPLE}*_marked_duplicates.bam`
echo ${FILE}
samtools view -h -b -F 4 -f 0x2 -q 30 ${FILE} > ${BAM_DIR}/${SAMPLE}_aligned_filtered.bam


done
