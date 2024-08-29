#!/bin/bash

# 29.06.2022

# using syri to call varaints based on the whole genome alignment produced by anchorwave

# source ~/anaconda3/etc/profile.d/conda.sh
# conda activate syri_env

# path to assemblies
hap1=Rabiosa_hap1_chr.fa
hap2=Rabiosa_hap2_chr.fa

# path to paf
mypaf=Rabiosa_all_chr.paf

# run syri
# syri -c $mypaf -r $hap1 -q $hap2 -F P -f --prefix Rabiosa_ --cigar --nc 7

# extract different types of variants from syri output
for x in $(cut -f11 Rabiosa_syri.out | sort | uniq)
do
    echo ${x}
    grep -w "${x}" Rabiosa_syri.out | cut -f1,2,3,6,7,8 > Rabiosa_variant_${x}.txt
done

