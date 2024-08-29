This folder contains the code used to analyse the sequence variation between haplotypes

```
# for call variants from whole-genome alignment, see myanchorwave_all_gene.sh and mysyri.sh
# for map short or long reads to the assembly and variant calling from reads, see map_short_snakemake and map_long_snakemake

# for later parsing syri output and finding common variants, see below

# custom python script to parsy syri variant calling results
./parse_syri_out.py Rabiosa_syri.out Rabiosa_syri.snp Rabiosa_syri.indel Rabiosa_syri.pav

# filter SVs with Ns
grep -v 'NN' Rabiosa_syri.pav | cut -f1,2,3,6,7 > Rabiosa_syri.pav.filtered

# extract different types of variants from syri output
for x in $(cut -f11 Rabiosa_syri.out | sort | uniq)
do
    echo ${x}
    grep -w "${x}" Rabiosa_syri.out | cut -f1,2,3,6,7,8 > Rabiosa_variant_${x}.txt
done

# find common SV
Rscript intersect_SV_modified.R Rabiosa_syri.pav.filtered sniffles_called_sv.txt

# find common SNP and INDEL
Rscript Find_common_small_variants.R SNP_parsed.txt Rabiosa_syri.snp INDEL_parsed.txt Rabiosa_syri.indel

# find the position of gaps
./Identify_gap.py rabiosa_hap1_chr.fa hap1_gap.txt
./Identify_gap.py rabiosa_hap2_chr.fa hap2_gap.txt
```
