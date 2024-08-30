# find the peak and define the candidate region

# set working directory
setwd("P:/Yutangchen/GBS_snakemake/stemrustqtl")

# import lod table
lod_unphased <- read.table("unphased/Rabiosa_unphased_qtl_lod_all.txt", header = T, stringsAsFactors = F)
lod_hap1 <- read.table("hap1/Rabiosa_hap1_qtl_lod_all.txt", header = T, stringsAsFactors = F)
lod_hap2 <- read.table("hap2/Rabiosa_hap2_qtl_lod_all.txt", header = T, stringsAsFactors = F)
lod_vg <- read.table("vg/Rabiosa_vg_qtl_lod_all.txt", header = T, stringsAsFactors = F)

define_qtl_loci <- function(lod_data = lod_unphased, 
                            qtltable = "unphased/Rabiosa_unphased_qtl_loci.txt",
                            lod = "cim_lod"){
  
  mytrait <- unique(lod_data$trait)
  mychr <- unique(lod_data$chr)
  
  qtltrait <- NULL
  
  for(j in 1:length(mytrait)){
    
    # choose a trait
    mytest <- lod_data[lod_data$trait == mytrait[j], ]
    
    qtlchr <- NULL
    
    for(i in 1:length(mychr)){
      
      # subset significant markers for one chr
      mysig <- mytest[mytest$chr == mychr[i] & mytest[[lod]] >= unique(mytest$threshold), ]
      
      # find the gap between loci
      loci_border <- mysig[which(diff(mysig[mysig$chr == i, 2]) >= 10) + 1, ]
      
      # define loci
      loci_start <- c()
      loci_end <- c()
      
      for(x in 1:nrow(loci_border)){
        
        loci_start <- c(loci_start, which(mysig$loc == loci_border$loc[x]))
        loci_end <- c(loci_end, which(mysig$loc == loci_border$loc[x])-1)
        
      }
      
      loci_start <- c(1, loci_start)
      loci_end <- c(loci_end, nrow(mysig))
      
      qtl <- c()
      for(x in 1:length(loci_start)){
        
        qtl <- c(qtl, rep(x, loci_end[x]-loci_start[x]+1))
        
      }
      
      mysig$qtl <- qtl
      qtlchr <- rbind.data.frame(qtlchr, mysig)
      
    }
    
    qtltrait <- rbind.data.frame(qtltrait, qtlchr)
    
  }
  
  qtltrait$loc <- sub("chr.*_", "", qtltrait$loc)
  write.table(qtltrait, qtltable, quote = F, sep = "\t", row.names = F, col.names = T)
  
}

define_qtl_loci(lod_data = lod_unphased, qtltable = "unphased/Rabiosa_unphased_qtl_loci.txt")
define_qtl_loci(lod_data = lod_hap1, qtltable = "hap1/Rabiosa_hap1_qtl_loci.txt")
define_qtl_loci(lod_data = lod_hap2, qtltable = "hap2/Rabiosa_hap2_qtl_loci.txt")
define_qtl_loci(lod_data = lod_vg, qtltable = "vg/Rabiosa_vg_qtl_loci.txt")

