#!/bin/bash
input="./samples.txt"
while IFS= read -r line
do
	java -jar ~/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 ${line}_R1_001.fastq.gz ${line}_R2_001.fastq.gz ${line}_R1_paired.fastq ${line}_R1_unpaired.fastq ${line}_R2_paired.fastq ${line}_R2_unpaired.fastq ILLUMINACLIP:~/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
	bowtie2 -p 32 -q -x $1 -1 ${line}_R1_paired.fastq -2 ${line}_R2_paired.fastq -U ${line}_R1_unpaired.fastq,${line}_R2_unpaired.fastq -S ${line}_bowtie2.sam
	samtools view -b -@ 32 ${line}_bowtie2.sam > ${line}_bowtie2.bam
	samtools sort -@ 32 -o ${line}_bowtie2.sort.bam ${line}_bowtie2.bam
	samtools index -@ 32 -b ${line}_bowtie2.sort.bam ${line}_bowtie2.sort.bam.bai
	~/bbmap/pileup.sh in=${line}_bowtie2.sam basecov=${line}_bowtie2_coverage.txt delcoverage=f 32bit=t -Xmx64g
	lofreq call-parallel --pp-threads 32 -f $2 -d 100000 -o ${line}.vcf ${line}_bowtie2.sort.bam
done < "$input"