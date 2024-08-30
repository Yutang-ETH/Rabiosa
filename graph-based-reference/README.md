Code used for building graph-based reference

So, I run anchorwave to align Rabiosa unphased, h1 and h2.

Code for h1 vs h2 alignment is in haplotype_variation. So is the R code to concatenate paf file of each chromosome.

With the pair-wise alignment, build the pangenome graph using PGGB, see mypggb.sh

mypggb.sh generates a gfa file, and this gfa file was deconstructed to result in a vcf file

```
# deconstruct the gfa
vg deconstruct -p rabiosa_0_chr1 -p rabiosa_0_chr2 -p rabiosa_0_chr3 -p rabiosa_0_chr4 -p rabiosa_0_chr5 -p rabiosa_0_chr6 -p rabiosa_0_chr7 -H '?' -a -e -t 48 -v ${mygfa}/Rabiosa.gfa > rabiosa_deconstruct_all_vg143.vcf

# compress the vcf
bgzip -c -@ 48 rabiosa_deconstruct_all_vg143.vcf > rabiosa_deconstruct_all_vg143.vcf.gz

# filter vcf
${myvcfbub}/vcfbub -l 0 -r 10000 -i rabiosa_deconstruct_all_vg143.vcf.gz > rabiosa_deconstruct_all_vg143_vcfbub_filtered.vcf

# reformat the vcf
./reformatvcf.sh rabiosa_deconstruct_all_vg143_vcfbub_filtered.vcf

# rename the vcf
cut -f1 rabiosa_unphased_chr.fa.fai > name_1
cut -f1 rabiosa_unphased_chr.fa.fai | sed 's/rabiosa_0_//' > name_2
paste name_1 name_2 > change_chr_name.txt
rm name_1 name_2
bcftools annotate --rename-chrs change_chr_name.txt pangenie_snp.vcf.gz | bgzip > pangenie_snp_rename.vcf.gz

# with the renamed vcf, do genotyping
myvg_genotype_all_new.sh
```
