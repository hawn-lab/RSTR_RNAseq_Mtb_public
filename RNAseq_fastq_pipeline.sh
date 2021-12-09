#!/bin/bash
# Example fastq processing

#### Set parameters ####
sample="sample1"
read1="rstr-raw-fastq/sample1.read1.fastq.gz"
read2="rstr-raw-fastq/sample1.read2.fastq.gz"
threads=90
ref_release="102"
adapter_length="7"

##### Quality assessment 1 ##### 

fastqc $read1 -o fastqc/ -t $threads
fastqc $read2 -o fastqc/ -t $threads

##### Adapter removal ##### 
## Remove adapters 
## Remove reads with > 1 ambiguous base
## Trim ends until reach base with quality 30+
## Remove reads < 15 bp

AdapterRemoval --file1 $read1 --file2 $read2 \
    --basename fastq_trim/sample1 --gzip \
    --trim5p $adapter_length --maxns 1 --minlength 15 \
    --trimqualities --minquality 30 \
    --threads $threads

##### Quality assessment 2 ##### 

fastqc fastq_trim/"$sample".pair1.truncated.gz -o fastqc/ -t $threads
fastqc fastq_trim/"$sample".pair2.truncated.gz -o fastqc/ -t $threads

##### Alignment ##### 

mkdir STARref/
cd STARref/
sudo curl -O ftp://ftp.ensembl.org/pub/release-$ref_release/gtf/homo_sapiens/Homo_sapiens.GRCh38.$ref_release.gtf.gz
sudo curl -O ftp://ftp.ensembl.org/pub/release-$ref_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Make genome index
mkdir STARindex
STAR --runMode genomeGenerate \
     --genomeDir STARindex \
     --genomeFastaFiles STARref/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile STARref/Homo_sapiens.GRCh38.$ref_release.gtf \
     --sjdbOverhang 99 \
     --runThreadN $threads

cp ./Log.out STARindex

## Align with STAR
STAR --genomeDir STARindex \
         --readFilesIn fastq_trim/"$sample".pair1.truncated.gz fastq_trim/"$sample".pair2.truncated.gz \
         --readFilesCommand zcat \
         --outFileNamePrefix bam/$sample \
         --outSAMtype BAM SortedByCoordinate \
         --runThreadN $threads \
         --runRNGseed 8756

##### Assess alignments ##### 
#Get Picard ref 
mkdir PICARDref
cd PICARDref
sudo curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz
gunzip refFlat.txt.gz
## Remove chr in chromosome name to match ensembl alignment
sed 's/chr//' refFlat.txt > refFlat.ensembl.txt
cd ../..

# median CV of gene model coverage
java -XX:ParallelGCThreads=$threads \
        -jar ~/project/apps/anaconda3/share/picard-2.25.2-0/picard.jar \
        CollectRnaSeqMetrics \
        REF_FLAT=PICARDref/refFlat.ensembl.txt \
        I=bam/"$sample".bam  \
        O=bam_metrics/"$sample".picard.tsv \
        ASSUME_SORTED=true \
        STRAND_SPECIFICITY=NONE \
        MINIMUM_LENGTH=500

samtools flagstat -@ $threads bam/"$sample".bam > bam_metrics/"$sample".stat.tsv

##### Quality filter BAM ##### 
## -h: keep header
## -f 3: keeps paired reads where both mapped
## -F 1284: removes unmapped reads, non-primary alignments, and PCR duplicates
## -q 30: min MAPQ of 30

## Paired reads
samtools view bam/"$sample".bam \
      -h -f 3 -F 1284 -q 30 \
      -@ $threads \
      > bam/"$sample".filter.bam

##### Assess filtered alignments #####

samtools flagstat -@ $threads bam/"$sample".filter.bam \
    >> bam_metrics/"$sample".filter.stat.tsv

##### Count reads in genes ##### 

featureCounts -T 64 -g gene_id -t exon -p \
  -a STARref/Homo_sapiens.GRCh38.$ref_release.gtf \
  -o counts/"$sample".featurecounts.paired.tsv \
  bam/"$sample".filter.bam

##### END ##### 