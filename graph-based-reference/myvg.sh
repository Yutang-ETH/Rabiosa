#!/bin/bash

# 05.10.2022

# map short reads to pangenome graph using giraffe

mygfa="/scratch/yutang/giraffe/pangenome"
myvcfbub="/home/yutachen/public/Yutangchen/vcfbub/target/release"


# first autoindex the graph for alignment
# vg autoindex --workflow giraffe --prefix Rabiosa --gfa ${mygfa} --threads 48 --tmp-dir ./tmp -V 1 

# convert giraffe.gbz format to xg format as suggested by a thread on vg github issue page #3694
# export TMPDIR=/scratch/yutang/giraffe/tmp 
# vg convert -x -t 48 Rabiosa.giraffe.gbz > Rabiosa.xg

# now map short reads to graph using giraffe
# mkdir mygam
# cat sample.lst | parallel -j 4 -k "vg giraffe -Z Rabiosa.giraffe.gbz -m Rabiosa.min -d Rabiosa.dist -f GBS/{}.fastq.trim.fq -t 12 > mygam/{}.gam"

# mkdir mybam
# cat sample.lst | parallel -j 3 -k "vg giraffe -Z Rabiosa.giraffe.gbz -m Rabiosa.min -d Rabiosa.dist -f GBS/{}.fastq.trim.fq -t 10 -o BAM --sample {} rabiosa_unphased_chr.fa > mybam/{}.bam"

# mkdir mybam 
# cat sample.lst | parallel -j 3 -k "vg giraffe -Z Rabiosa.giraffe.gbz -m Rabiosa.min -d Rabiosa.dist -f GBS/{}.fastq.trim.fq -t 12 -o BAM --ref-paths ref_paths.txt > mybam/{}.bam"

# now call variants for every gam file
# vg pack converts alignment to coverage
# mkdir mypack
# cat sample.lst | parallel -j 3 -k "vg pack -x Rabiosa.xg -g mygam/{}.gam -Q 5 -s 5 -o mypack/{}.pack -t 5"

# calculate snarls
# vg snarls -t 48 Rabiosa.xg > mysnarls.pb

# now vg call to output variants to vcf
# mkdir myvcf
# cat sample.lst | parallel -j 3 -k "vg call Rabiosa.xg -k mypack/{}.pack -a -t 5 -r mysnarls.pb -g Rabiosa.giraffe.gbz -s {} > myvcf/{}.vcf"

# unfortunately, above pipeline doesn't work efficiently, try to deconstruct the graph first to make vcf and then genotype snps in vcf with the samples from GBS data
# vg deconstruct -p rabiosa_0_chr7 -a -e -t 48 -v ${mygfa}/rabiosa_chr7.gfa | bgzip -c -@ 16 > rabiosa_deconstruct_chr7.vcf.gz

# vg deconstruct -p rabiosa_0_chr1 -p rabiosa_0_chr2 -p rabiosa_0_chr3 -p rabiosa_0_chr4 -p rabiosa_0_chr5 -p rabiosa_0_chr6 -p rabiosa_0_chr7 -a -e -t 48 --ploidy 1 -v ${mygfa}/Rabiosa.gfa > rabiosa_deconstruct_unphased.vcf
# vg deconstruct -p rabiosa_0_chr1 -p rabiosa_0_chr2 -p rabiosa_0_chr3 -p rabiosa_0_chr4 -p rabiosa_0_chr5 -p rabiosa_0_chr6 -p rabiosa_0_chr7 -H '?' -a -e -t 48 -v ${mygfa}/Rabiosa.gfa > rabiosa_deconstruct_all_vg140.vcf 
# vg deconstruct -p rabiosa_0_chr1 -p rabiosa_0_chr2 -p rabiosa_0_chr3 -p rabiosa_0_chr4 -p rabiosa_0_chr5 -p rabiosa_0_chr6 -p rabiosa_0_chr7 -H '?' -a -e -t 48 -v ${mygfa}/Rabiosa.gfa > rabiosa_deconstruct_all_vg143.vcf

# bgzip -c -@ 48 rabiosa_deconstruct_all_vg143.vcf > rabiosa_deconstruct_all_vg143.vcf.gz
${myvcfbub}/vcfbub -l 0 -r 10000 -i rabiosa_deconstruct_all_vg143.vcf.gz > rabiosa_deconstruct_all_vg143_vcfbub_filtered.vcf
