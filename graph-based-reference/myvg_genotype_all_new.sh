#!/bin/bash

# 17.02.2023

# genotype SNPs using vg 

myfa="rabiosa_unphased_chr_1_7.fa"
myvcf="pangenie_snp_rename.vcf.gz"
myprefix="Raball"

# first autoindex the graph for alignment
# vg autoindex --workflow giraffe --prefix ${myprefix} --ref-fasta ../${myfa} --vcf ../${myvcf} --threads 60 --tmp-dir ../tmp 

# giraffe mapping GBS reads to graph
# mkdir mygam
# cat ../sample.lst | parallel -j 5 -k "vg giraffe -Z ${myprefix}.giraffe.gbz -m ${myprefix}.min -d ${myprefix}.dist -f ../GBS/{}.fastq.trim.fq -t 12 -N {} > mygam/{}.gam"

# count read support
# mkdir mypack
cat ../sample.lst | parallel -j 5 -k "vg pack -x ${myprefix}.giraffe.gbz -t 12 -g mygam/{}.gam -o mypack/{}.pack -Q 5 -s 5"

# calculate snarl
vg snarls ${myprefix}.giraffe.gbz > ${myprefix}.snarls

# genotype the vcf
mkdir myvcf
cat ../sample.lst | parallel -j 10 -k "vg call ${myprefix}.giraffe.gbz -z -a -k mypack/{}.pack -s {} -t 6 -r ${myprefix}.snarls > myvcf/{}.vcf"

# compress vcf for merging
mkdir compressvcf
cat ../sample.lst | parallel -j 60 -k "bgzip -c myvcf/{}.vcf > compressvcf/{}.vcf.gz"

# index vcf
cat ../sample.lst | parallel -j 60 -k "tabix -p vcf compressvcf/{}.vcf.gz"

# concatenate vcf
bcftools merge compressvcf/rabiosa.vcf.gz compressvcf/lmGbsJP*.vcf.gz > ${myprefix}_merge.vcf

# filter vcf file
bcftools view -g het -v snps -m2 -M2 -q 0.03 -i 'F_MISSING<0.5' ${myprefix}_merge.vcf > ${myprefix}_merge_filtered.vcf

