# Install the human and Sendai virus (SeV) reference genome and merge 
- Install human genome (GRCh38):
http://ftp.ensembl.org/pub/release-105/gtf/homo_sapiens/


```bash
-Download the reference genome and human gene list (.gtf) file:
# the latest genome fasta file
$wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/latest/hg38.fa.gz

# gene annotation, "gtf", file:
$wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz
$gunzip hg38.ncbiRefSeq.gtf.gz

# check the line number in the gtf
wc -l ~/Gozde/Interferon_project/reference_genomes/human/hg38/hg38.ncbiRefSeq.gtf
#4344986

# Create the merged genome for human and SeV
cat ~/Gozde/Interferon_project/reference_genomes/human/hg38/hg38.fa ~/Gozde/Interferon_project/reference_genomes/SeV_Cantell/SeV_Cantell_strain_genome.fasta > ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genome.fasta

# Create the merged gtf file for gene info for human and SeV
cat ~/Gozde/Interferon_project/reference_genomes/human/hg38/hg38.ncbiRefSeq.gtf ~/Gozde/Interferon_project/reference_genomes/SeV_Cantell/SeV_Cantell_strain_genes.gtf > ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes.gtf

# check the line number in the new gtf
wc -l ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes.gtf
#4344998
```

# cellranger mkref
-Resource: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/count

- Need to filter out the non-polyA regions in the .gtf file, to prevent multi mapping of the protein coding reads.
Multi-mapped reads are not counted.
https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf

```bash
~/Downloads/softwares/cellranger-6.0.1/cellranger mkgtf ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes.gtf ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes_filtered.gtf --attribute=key:allowable_value
```

```bash
# Create the correct format of human+SeV genome using cellranger mkref
cd ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV
~/Downloads/softwares/cellranger-6.0.1/cellranger mkref --genome=humanhg38plusSeV --fasta=~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genome.fasta  --genes=~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes_filtered.gtf
```

#There was an error due to chromosome names not matching in fasta vs gtf.
#We need to add 'chr' to the beginning of each row containing gene info in gtf.
#```bash
#awk 'OFS="\t" {if (NR > 5) $1="chr"$1; print}' ~/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes.gtf > /thk-#linux/home/Gozde/Interferon_project/reference_genomes/merged/Hg38+SeV/human_hg38_SeV_genes_chr.gtf
#```

# make directory

```bash
# make new directory for mapping & counting results:
mkdir ~/interferon_project/APH_scRNA/Analysis_Results
cd ~/interferon_project/APH_scRNA/Analysis_Results
```

# map each sample of both runs to the custom humanplusSeV genome using cellranger in STARship:
```bash
# DMSO sample
~/interferon_project/softwares/cellranger-6.0.2/cellranger count --include-introns --id GM_DMSO_SeV_counts --localcores=12 --localmem=100  --fastqs ~/interferon_project/APH_scRNA/fastq/first_run/,~/interferon_project/APH_scRNA/fastq/second_run/ --sample GM-DMSO,2DMSO --transcriptome ~/interferon_project/APH_scRNA/humanhg38plusSeV/
```

```bash
# APH sample
~/interferon_project/softwares/cellranger-6.0.2/cellranger count --include-introns --id GM_APH_SeV_counts --localcores=12 --localmem=100  --fastqs ~/interferon_project/APH_scRNA/fastq/first_run/,~/interferon_project/APH_scRNA/fastq/second_run/ --sample GM-APH,1APH --transcriptome ~/interferon_project/APH_scRNA/humanhg38plusSeV/
```
