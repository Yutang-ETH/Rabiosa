# add phenotype to genetic linkage map

# setwd
setwd("P:/Yutangchen/GBS_snakemake")

# import stemrust genotype from Jenny to get stem rust phenotypes
phenotype_jenny <- read.csv("Jenny/Mergeallphenoemmeans.csv", header = T, stringsAsFactors = F)
phenotype_jenny <- phenotype_jenny[, -1]


# function to add phenotypes to the genotype dataset
add_phenotype <- function(genotype = genotype, 
                          sample_name = sample_name,
                          pheno = mypheno,
                          output = output){
  
  # remove parents columns, col 4 and 5
  genotype <- genotype[, -c(4, 5)]
  
  # remove rabiosa from sample list
  sample_name <- as.data.frame(sample_name[-1, ])
  colnames(sample_name) <- "sample"
  
  # add ID to sample_pan
  sample_name$ID <- gsub("_", "", sample_name$sample)
  sample_name$ID <- as.numeric(sub("lmGbsJP", "", sample_name$ID))
  
  # get corresponding phenotypes based on the ID
  sample_pheno <- na.omit(pheno[match(sample_name$ID, pheno$PlantID), ])
  
  # add col names to genotype table
  colnames(genotype) <- c("pheno", "chr", "map", sample_name$ID)
  
  # subset genotype table based on the samples with phenotypes
  genotype_sub <- genotype[colnames(genotype) %in% c("pheno", "chr", "map", sample_pheno$PlantID)]
  
  # add phenotypes to genotype_unphased_sub
  pheno_all <- as.data.frame(t(rbind.data.frame(colnames(sample_pheno)[2:ncol(sample_pheno)],
                                                rep("", ncol(sample_pheno)-1), 
                                                rep("", ncol(sample_pheno)-1),
                                                sample_pheno[, 2:ncol(sample_pheno)])))
  colnames(pheno_all) <- colnames(genotype_sub)
  genotype_qtl <- rbind.data.frame(pheno_all, genotype_sub[-1, ])  
  
  # save genotype_unphased_qtl as csv file for Rqtl
  write.table(genotype_qtl, output, sep = ",", quote = F, row.names = F, col.names = F)
  
}

#-----------------------------------------------------------------------------------------------------#

# output the unphased map for qtl mapping, this is rabiosa only map

# import genotype 
genotype_unphased <- read.csv("rabiosa_f/lepmap_phase/mycsvr_rabiosa.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_unphased <- read.table("rabiosa_f/lepmap/my_final_sample.txt", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_unphased, 
              sample_name = sample_unphased, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_unphased_rabiosa_only_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#

# output the unphased map for qtl mapping, the consensus map containing informative markers of both parents (sikem and rabiosa)

# import genotype 
genotype_unphased <- read.csv("rabiosa_f/lepmap_consensus_map/mycsvr_rabiosa_consensus.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_unphased <- read.table("rabiosa_f/lepmap_consensus_map/my_final_sample.txt", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_unphased, 
              sample_name = sample_unphased, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_unphased_consensus_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#

# output the unphased map for qtl mapping, the sikem only map

# import genotype 
genotype_unphased <- read.csv("rabiosa_f/lepmap_sikem_map/mycsvr_sikem_only.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_unphased <- read.table("rabiosa_f/lepmap_sikem_map/my_final_sample.txt", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_unphased, 
              sample_name = sample_unphased, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_unphased_sikem_only_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#

# output the hap1 map for qtl mapping

# import genotype 
genotype_hap1 <- read.csv("hap1_f/lepmap_phase/mycsvr_rabiosa.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_hap1 <- read.table("hap1_f/lepmap/my_final_sample.txt", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_hap1, 
              sample_name = sample_hap1, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_hap1_pheno_all_qtl.csv")



#-----------------------------------------------------------------------------------------------------#

# output the hap2 map for qtl mapping

# import genotype 
genotype_hap2 <- read.csv("hap2_f/lepmap_phase/mycsvr_rabiosa.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_hap2 <- read.table("hap2_f/lepmap/my_final_sample.txt", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_hap2, 
              sample_name = sample_hap2, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_hap2_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#

# output the pangenie map for qtl mapping

# import genotype 
genotype_pangenie <- read.csv("pan/pangenie/mycsvr_rabiosa.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_pangenie <- read.table("pan/pangenie/mysample.lst", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_pangenie, 
              sample_name = sample_pangenie, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_pangenie_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#

# output the vg map for qtl mapping

# import genotype 
genotype_vg <- read.csv("pan/vg_giraffe_genotype/vg_all/mycsvr_rabiosa.csv", header = F, stringsAsFactors = F)

# get sample name in the unphased rabiosa genotype table
sample_vg <- read.table("pan/vg_giraffe_genotype/mysample.lst", header = F, stringsAsFactors = F)

# output the csvr tabel
add_phenotype(genotype = genotype_vg, 
              sample_name = sample_vg, 
              pheno = phenotype_jenny,
              output = "stemrustqtl/rabiosa_vg_pheno_all_qtl.csv")

#-----------------------------------------------------------------------------------------------------#