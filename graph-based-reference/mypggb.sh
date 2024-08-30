#!/bin/bash

mypggb=/home/yutachen/pggb-0.3.1
mydip="rabiosa_dip_chr.fa"
mypaf1="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave/chr_paf_pan"
mypaf2="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap1/chr_paf_new"
mypaf3="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap2/chr_paf_new"

# build the pangenome graph for hap1, hap2 and unphased, I think this will be the final pangenome for Rabiosa
# for i in $(seq 1 7)
# do
#     grep "rabiosa_[012]_chr${i}" ${mydip} | sed 's/>//' > name.lst
#     seqtk subseq ${mydip} name.lst > rabiosa_hap_chr${i}.fa
#     samtools faidx rabiosa_hap_chr${i}.fa
#     cat ${mypaf1}/chr${i}.paf ${mypaf2}/chr${i}.paf ${mypaf3}/chr${i}.paf > rabiosa_dip_${i}.paf
#     $mypggb/pggb -i rabiosa_hap_chr${i}.fa -n 3 -o chr_test_${i} -t 48 -T 40 -P asm20 -G 23117,23219 -a rabiosa_dip_${i}.paf 
# done

# combine pangenome graphs
# first make a list of graph files, one graph per line

# rm pangenome.lst
# touch pangenome.lst
# chmod 777 pangenome.lst

# for i in $(seq 1 7)
# do
#     ls chr_test_${i}/*.smooth.final.og >> pangenome.lst
# done    

# now combine using odgi squeeze
# mkdir chr_all
# odgi squeeze -f pangenome.lst -o chr_all/Rabiosa.og -t 48 -P 

# now convert og to gfa format 
# odgi view -i chr_all/Rabiosa.og -g -a -t 48 -P > chr_all/Rabiosa.gfa

# now 1D visualize the combined og 
# odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc.png -t 48
# odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y500_a10.png -t 48 -x 3000 -y 500 -a 10
# odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y1000_a10.png -t 48 -x 3000 -y 1000 -a 10
# odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y1500_a20.png -t 48 -x 3000 -y 1500 -a 20

# now 2D visualize the combinded og
# first generate the layout file
# odgi layout -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.og.lay -t 48 -T chr_all/Rabiosa.all.og.lay.tsv
# odgi draw -i chr_all/Rabiosa.og -c chr_all/Rabiosa.all.og.lay -p chr_all/Rabiosa.all.2D.png -C -w 20 -H 1000


# detect gene PAV with the bed derived from gff and also translate position of a path to position on other paths
# mkdir chr_pav
# for i in $(seq 1 7)
# do
#     cut -f 1,2,3 ~/public/Yutangchen/Rabiosa_compare/All/all/Rabiosa_0_chr.bed | grep "chr${i}" | sed 's/Lm0/rabiosa_0/' > chr_pav/chr${i}.multiple.references.gene.bed
#     cut -f 1,2,3 ~/public/Yutangchen/Rabiosa_compare/All/all/Rabiosa_1_chr.bed | grep "chr${i}" | sed 's/Lm1/rabiosa_1/' >> chr_pav/chr${i}.multiple.references.gene.bed
#     cut -f 1,2,3 ~/public/Yutangchen/Rabiosa_compare/All/all/Rabiosa_2_chr.bed | grep "chr${i}" | sed 's/Lm2/rabiosa_2/' >> chr_pav/chr${i}.multiple.references.gene.bed
#     odgi pav -i chr_test_${i}/*.smooth.final.og -b chr_pav/chr${i}.multiple.references.gene.bed -t 48 > chr_pav/chr${i}.gene.pavs.txt
    # odgi position -i chr_test_${i}/*.smooth.final.og -b chr_pav/chr${i}.multiple.references.gene.bed -r rabiosa_0_chr${i} -t 48 > chr_pav/rabiosa_0_chr${i}.position.bed       
    # odgi position -i chr_test_${i}/*.smooth.final.og -b chr_pav/chr${i}.multiple.references.gene.bed -r rabiosa_1_chr${i} -t 48 > chr_pav/rabiosa_1_chr${i}.position.bed
    # odgi position -i chr_test_${i}/*.smooth.final.og -b chr_pav/chr${i}.multiple.references.gene.bed -r rabiosa_2_chr${i} -t 48 > chr_pav/rabiosa_2_chr${i}.position.bed
# done

# touch chr_pav/all.gene.pavs.txt
# chmod 777 chr_pav/all.gene.pavs.txt
# for i in {1..7}
# do
#     grep '^rabiosa_0_' chr_pav/chr${i}.gene.pavs.txt >> chr_pav/all.gene.pavs.txt 
# done

for i in {1..7}
do
    bedtools makewindows -g <(grep "rabiosa_0_chr${i}" rabiosa_dip_chr.fa.fai | cut -f 1,2) -w 100000 > chr_pav/rabiosa.chr${i}.w100kbp.bed
    odgi pav -i chr_test_${i}/*.smooth.final.og -b chr_pav/rabiosa.chr${i}.w100kbp.bed -t 48 > chr_pav/chr${i}.window.pavs.txt
done

touch chr_pav/all.window.pavs.txt
chmod 777 chr_pav/all.window.pavs.txt
for i in {1..7}
do
    grep '^rabiosa_0_' chr_pav/chr${i}.window.pavs.txt >> chr_pav/all.window.pavs.txt
done

cp chr_pav/all.window.pavs.txt ~/public/Yutangchen/pandip/chr_pav


# extract pav ratio based on rabiosa_0_chr*
# for i in $(seq 1 7)
# do
#     echo chr${i}
#     grep "^rabiosa_0_chr${i}" chr_pav/chr${i}.gene.pavs.matrix.txt >> chr_pav/chrx.hap0.gene.pavs.matrix.txt
#     grep "^rabiosa_0_chr${i}" chr_pav/rabiosa_1_chr${i}.position.bed >> chr_pav/rabiosa_1_chrx.position.bed
#     grep "^rabiosa_0_chr${i}" chr_pav/rabiosa_2_chr${i}.position.bed >> chr_pav/rabiosa_2_chrx.position.bed
# done

# extract only hap1 and hap2 paths and save in a new graph
# rm hap1_hap2_graph.lst
# touch hap1_hap2_graph.lst
# chmod 777 hap1_hap2_graph.lst

# for i in $(seq 1 7)
# do
#     odgi paths -i chr${i}/*.smooth.final.og -L | tail -n 2 > chr${i}/mypath.txt
#     odgi extract -i chr${i}/*.smooth.final.og -p chr${i}/mypath.txt -o chr${i}/chr${i}_hap1_hap2_graph.og -t 48 -c 0 -E
#     odgi viz -i chr${i}/chr${i}_hap1_hap2_graph.og -o chr${i}/chr${i}_hap1_hap2_graph_1D.png
#     ls chr${i}/chr${i}_hap1_hap2_graph.og >> hap1_hap2_graph.lst
# done

# odgi squeeze -f hap1_hap2_graph.lst -o chr_all/Rabiosa_hap1_hap2.og -t 48 -P
# odgi view -i chr_all/Rabiosa_hap1_hap2.og -g -a -t 48 -P > chr_all/Rabiosa_hap1_hap2.gfa



