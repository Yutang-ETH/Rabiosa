#!/bin/bash

myvcf=$1

for i in {1..7}
do
    grep -v '#' ${myvcf} | grep "chr${i}" | cut -f $(expr 9 + ${i}) > chr${i}_1.geno
    grep -v '#' ${myvcf} | grep "chr${i}" | cut -f $(expr 16 + ${i}) > chr${i}_2.geno
    paste -d "|" chr${i}_1.geno chr${i}_2.geno > chr${i}.geno
    rm chr${i}_1.geno chr${i}_2.geno 
done

cat chr*.geno > mygeno
rm chr*.geno

grep '##' ${myvcf} > meta_header
grep '#' ${myvcf} | tail -n 1 | cut -f1-9 > field_header 
echo "rabiosa" > mysample
paste field_header mysample > new_field_header

grep -v "#" ${myvcf} | cut -f1-9 > variant
paste variant mygeno > newvariant

cat meta_header new_field_header newvariant > pangenie.vcf

rm meta_header field_header mysample new_field_header variant newvariant mygeno

# extract only SNPs
bcftools view -v snps -m2 -M2 pangenie.vcf > pangenie_snp.vcf
