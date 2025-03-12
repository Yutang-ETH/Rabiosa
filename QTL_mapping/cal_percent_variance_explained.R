# calculate %var of the QTL 

# set working directory
setwd("P:/Yutangchen/GBS_snakemake/stemrustqtl")

# import genotype and phenotype table
data_unphased <- read.cross("csvr",".","rabiosa_unphased_pheno_sr_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
data_vg <- read.cross("csvr",".","rabiosa_vg_pheno_sr_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))

data_jitter <- jittermap(data_unphased)
data_jitter <- jittermap(data_vg)
crossInt1 <- calc.genoprob(data_jitter, step = 1, error.prob = 0.001)

Cim<-cim(crossInt1, n.marcovar = 5, window=10)

# the SNPs / map position with the highest LOD, chr2_231403834 / 102.05846, chr7_249653003 / 150.7636
# the SNPs / map position with the highest LOD, chr2_231403793 / 77.86411, chr7_242856925 / 145.63713

sr <- sim.geno(data_jitter)
qtl <- makeqtl(sr, chr = c(2, 7), pos = c(102.05846, 150.7636))
plot(qtl)
out.fq <- fitqtl(sr, qtl = qtl, formula = y~Q1+Q2)
summary(out.fq)

# results based on SNPs from single reference
'
Drop one QTL at a time ANOVA table: 
----------------------------------  
        df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
2@102.1  1       34.38  3.907  5.931   18.92            0  2.92e-05 ***
7@150.8  1      273.47 21.552 47.179  150.54            0   < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'

# results based on SNPs from graph-based reference
'
Drop one QTL at a time ANOVA table: 
  ----------------------------------  
  df Type III SS    LOD   %var F value Pvalue(Chi2) Pvalue(F)    
2@77.9   1       61.22  4.829  9.127   23.67            0  3.27e-06 ***
  7@145.6  1      236.79 15.373 35.302   91.55            0   < 2e-16 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
'