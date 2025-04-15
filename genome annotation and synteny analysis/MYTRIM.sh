#!/bin/bash

# 29.01.2022

# trim the raw RNA seq reads, remove adapters and trim the left 10 biased bps

# first make a name list of all the fastq files
# ls NextSeq500_20180307_NS163_o4112_DataDelivery/2018* | cut -f2 -d'/' | cut -f2 -d'.' | cut -f2 -d'_' > tissues.txt
# ls NextSeq500_20180307_NS163_o4112_DataDelivery/2018* | cut -f2 -d'/' | cut -f2 -d'.' | cut -f3 -d'_' > replicate.txt
# paste tissues.txt replicate.txt -d'_' | uniq > fastq.names
# rm tissues.txt replicate.txt

# make a new folder to store trimmed reads
# mkdir yutang_trimmomatic

# start trimmiing with trimmomatic
# cat fastq.names | parallel -j 12 -k "trimmomatic PE -validatePairs -threads 5 -phred33 NextSeq500_20180307_NS163_o4112_DataDelivery/20180307.A-New_{}_R1.fastq.gz NextSeq500_20180307_NS163_o4112_DataDelivery/20180307.A-New_{}_R2.fastq.gz -baseout yutang_trimmomatic/{}_trimmed.fq.gz ILLUMINACLIP:True_seq_adapter.fa:2:30:10 HEADCROP:10 MINLEN:60"

# do fastqc with trimmed reads
# mkdir yutang_trimmomatic_fastqc
# cat fastq.names | parallel -j 12 -k "fastqc -t 5 -o yutang_trimmomatic_fastqc yutang_trimmomatic/{}_trimmed_1P.fq.gz yutang_trimmomatic/{}_trimmed_2P.fq.gz" 

# trimmomatic can remove adaptors but there are still some poly G or poly X sequencing remaining in the reads
# try fastp to trim the raw reads

# mkdir yutang_fastp
cat fastq.names | parallel -j 12 -k "fastp -i NextSeq500_20180307_NS163_o4112_DataDelivery/20180307.A-New_{}_R1.fastq.gz -I NextSeq500_20180307_NS163_o4112_DataDelivery/20180307.A-New_{}_R2.fastq.gz -o yutang_fastp/{}_trimmed_R1.fq.gz -O yutang_fastp/{}_trimmed_R2.fq.gz -g -y -x -f 10 -F 10 -w 5 -Q -l 40 -5 -3"

# do fastqc with trimmed reads
# mkdir yutang_fastp_fastqc
cat fastq.names | parallel -j 12 -k "fastqc -t 5 -o yutang_fastp_fastqc yutang_fastp/{}_trimmed_R1.fq.gz yutang_fastp/{}_trimmed_R2.fq.gz"

