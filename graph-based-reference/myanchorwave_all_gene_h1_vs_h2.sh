#!/bin/bash

# 04.09.2022

# whole genome alignment is the key to compare the two haplotypes and it is crutial for conscruting pangenome
# considering the high divergince between haplotypes, using anchorwave to do whole genome alignment as it relies on gene collinearity

# prepare data
# get chr from gff3
# copy unphased gene.gff3 here in this folder
# cp /home/yutachen/public/Yutangchen/Rabiosa_annotation/EVM/filter/EVM.gene.gff3 ./
# grep '^chr' EVM.gene.gff3 > unphased.chr.gene.gff3

# prepare unphased.chr.fa and hap1.chr.fa
# mkdir assembly
# cp /home/yutachen/public/Yutangchen/Rabiosa_anchorwave/assembly/rabiosa_hap2_chr.fa ./assembly/
# cp /home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap1/assembly/rabiosa_unphased_chr.fa ./assembly/

# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate seqtk

# now lift over gff

source ~/anaconda3/etc/profile.d/conda.sh
conda activate anchorwave

# anchorwave gff2seq -i unphased.chr.gene.gff3 -r assembly/rabiosa_unphased_chr.fa -m 20 -o cds.fa
# gmap_build --dir=./hap2 --genomedb=hap2 assembly/rabiosa_hap2_chr.fa
# gmap -t 50 -A -f samse -d hap2 -D hap2/ cds.fa > Rabiosa_hap2.sam

# cp /home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap1/Rabiosa_unphased.sam ./Rabiosa_unphased.sam

# whole genome alignment
# anchorwave genoAli -i unphased.chr.gene.gff3 -as cds.fa -r assembly/rabiosa_unphased_chr.fa -a Rabiosa_hap2.sam -ar Rabiosa_unphased.sam -s assembly/rabiosa_hap2_chr.fa -t 30 -IV -n Rabiosa_hap.anchors -o Rabiosa_hap.maf -f Rabiosa_hap.f.maf -m 20 > Rabiosa_hap.log

# now convert maf to sam 
python2 maf-convert sam Rabiosa_hap.maf | sed 's/[0-9]\+H//g' > Rabiosa_hap_anchorwave.sam
cat Rabiosa_hap_anchorwave.sam | samtools view -O BAM --reference assembly/Rabiosa_unphased_chr.fa - | samtools sort -@ 60 - > Rabiosa_hap_anchorwave.bam
samtools index Rabiosa_hap_anchorwave.bam

# convert sam to paf for constructing the pangenome 
mkdir chr_sam chr_paf
for i in {1..7}
do 
    samtools view -h --reference assembly/Rabiosa_unphased_chr.fa Rabiosa_hap_anchorwave.bam chr${i} > chr_sam/chr${i}.sam
    paftools.js sam2paf chr_sam/chr${i}.sam > chr_paf/chr${i}.paf
done 

cat Rabiosa_hap.maf | python2 anchorwave-maf-swap.py > anchorwave_swap.maf
maf-convert sam anchorwave_swap.maf | sed 's/[0-9]\+H//g' > anchorwave_swap.sam
samtools view -O BAM --reference assembly/Rabiosa_hap2_chr.fa anchorwave_swap.sam | samtools sort -@ 60 - > anchorwave_swap.bam
samtools index -@ 60 anchorwave_swap.bam

for i in {1..7}
do
    samtools view -h --reference assembly/Rabiosa_hap2_chr.fa anchorwave_swap.bam chr${i} > chr_sam/chr${i}_swap.sam
    paftools.js sam2paf chr_sam/chr${i}_swap.sam > chr_paf/chr${i}_swap.paf
done

source ~/anaconda3/etc/profile.d/conda.sh
conda activate tritex

# add query coordinate to paf
mkdir chr_paf_new
for i in {1..7}
do
    Rscript Integrate_paf.R chr_paf/chr${i}.paf chr_paf/chr${i}_swap.paf "rabiosa_2_" "rabiosa_0_" chr_paf_new/chr${i}.paf 
done

# combine every chr paf to one paf
cat chr_paf_new/*.paf > Rabiosa_unphased_hap2_chr.paf
