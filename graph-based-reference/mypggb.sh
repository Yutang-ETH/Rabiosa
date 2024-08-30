#!/bin/bash

mypggb=/home/yutachen/pggb-0.3.1
mydip="rabiosa_dip_chr.fa"
mypaf1="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave/chr_paf_pan"
mypaf2="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap1/chr_paf_new"
mypaf3="/home/yutachen/public/Yutangchen/Rabiosa_anchorwave_unphased_hap2/chr_paf_new"

# build the pangenome graph for hap1, hap2 and unphased, I think this will be the final pangenome for Rabiosa
for i in $(seq 1 7)
do
    grep "rabiosa_[012]_chr${i}" ${mydip} | sed 's/>//' > name.lst
    seqtk subseq ${mydip} name.lst > rabiosa_hap_chr${i}.fa
    samtools faidx rabiosa_hap_chr${i}.fa
    cat ${mypaf1}/chr${i}.paf ${mypaf2}/chr${i}.paf ${mypaf3}/chr${i}.paf > rabiosa_dip_${i}.paf
    $mypggb/pggb -i rabiosa_hap_chr${i}.fa -n 3 -o chr_test_${i} -t 48 -T 40 -P asm20 -G 23117,23219 -a rabiosa_dip_${i}.paf 
done

# combine pangenome graphs
# first make a list of graph files, one graph per line

rm pangenome.lst
touch pangenome.lst
chmod 777 pangenome.lst

for i in $(seq 1 7)
do
    ls chr_test_${i}/*.smooth.final.og >> pangenome.lst
done    

# now combine using odgi squeeze
mkdir chr_all
odgi squeeze -f pangenome.lst -o chr_all/Rabiosa.og -t 48 -P 

# now convert og to gfa format 
odgi view -i chr_all/Rabiosa.og -g -a -t 48 -P > chr_all/Rabiosa.gfa

# now 1D visualize the combined og 
odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc.png -t 48
odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y500_a10.png -t 48 -x 3000 -y 500 -a 10
odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y1000_a10.png -t 48 -x 3000 -y 1000 -a 10
odgi viz -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.final.og.viz_multiqc_x3000_Y1500_a20.png -t 48 -x 3000 -y 1500 -a 20

# now 2D visualize the combinded og
# first generate the layout file
odgi layout -i chr_all/Rabiosa.og -o chr_all/Rabiosa.all.og.lay -t 48 -T chr_all/Rabiosa.all.og.lay.tsv
odgi draw -i chr_all/Rabiosa.og -c chr_all/Rabiosa.all.og.lay -p chr_all/Rabiosa.all.2D.png -C -w 20 -H 1000
