#!/bin/bash

SNP=$1
missing_value=$2

grep "^Missing" $SNP | cut -f2 -d ":" > missing

grep "^Sample" $SNP | cut -f2 -d ":" | tail -n +2 > name

grep "^Sample" $SNP | cut -f2 -d ":" | head -n 1 > parent

paste name missing > missing_rate.txt

rm name missing

awk -v var="$missing_value" '{ if($2 <= var) { print }}' missing_rate.txt | cut -f1 > my_final_sample

cat parent my_final_sample > my_final_sample.txt

sed -i 's/ //g' my_final_sample.txt

rm parent my_final_sample 

wc -l my_final_sample.txt
