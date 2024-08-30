#!/bin/bash

# calculate distance with given order

mypath=/home/yutachen/public/Yutangchen/lepmap3/bin
vcf=rabiosa_30_final.vcf

# import data
java -cp $mypath ParentCall2 data=pedigree.txt vcfFile=${vcf} > p.call

# filter data
java -cp $mypath Filtering2 data=p.call dataTolerance=0.000001 removeNonInformative=1 > p_f.call

# make physical order
grep 'chr' p_f.call | cut -f1 > chr.txt
seq 1 $(wc -l chr.txt | cut -f1 -d " ") > physical_order.txt
paste physical_order.txt chr.txt > po.txt
rm physical_order.txt chr.txt

# calculate the distance
for i in $(seq 1 7)
do
    grep "chr$i" po.txt | cut -f1 > mypo_$i.txt 
    java -cp $mypath OrderMarkers2 evaluateOrder=mypo_$i.txt data=p_f.call improveOrder=0 numThreads=40 chromosome=$i informativeMask=2 outputPhasedData=1 >order_$i.txt 2>order_$i.err 
    echo "chr$i finished"
done

# cat chroms
cat order_*.txt > order_x.txt

# convert map to R/qtl format
awk -vfullData=1 -f map2genotypes.awk order_x.txt > genotypes_x.txt

# extract the row number of each SNP in the p.call file
cut -f 1,2 p_f.call | awk '(NR>=7)' > snps.txt

# match marker number to position
awk -vFS="\t" -vOFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' snps.txt order_x.txt > order_x.mapped

# get maker position and the map and prepare the input for allmap
grep -v "#" order_x.mapped | cut -f 1-4 > scaffold_position_map.txt

# I copied these commands from WHATIDID.txt
# subset vcf file with informative SNPs in the map
cut -f1,2 scaffold_position_map.txt > snps_in_map.txt
bgzip -c rabiosa_30_final.vcf > rabiosa_30_final.vcf.gz
tabix -p vcf rabiosa_30_final.vcf.gz
bcftools view -R snps_in_map.txt rabiosa_30_final.vcf.gz > snps_in_map.vcf

bcftools view -H snps_in_map.vcf | cut -f1,2,4,5 > phase_in_assembly.txt

