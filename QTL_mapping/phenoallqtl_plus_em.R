# doing qtl mapping for all phenotypes Jenny shared with me

# set working directory
# setwd("P:/Yutangchen/GBS_snakemake/stemrustqtl")

# load library
library(qtl)

# import genotype and phenotype table
# data_unphased <- read.cross("csvr",".","rabiosa_unphased_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
# data_hap1 <- read.cross("csvr",".","rabiosa_hap1_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
# data_hap2 <- read.cross("csvr",".","rabiosa_hap2_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
# data_pangenie <- read.cross("csvr",".","rabiosa_pangenie_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
# data_vg <- read.cross("csvr",".","rabiosa_vg_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))

# data_unphased_consensus <- read.cross("csvr",".","rabiosa_unphased_consensus_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))
data_unphased_sikem <- read.cross("csvr",".","rabiosa_unphased_sikem_only_pheno_all_qtl.csv",na.strings=NA,genotypes=c("A","H","B"))

myqtl <- function(output_datasummary = output_datasummary,
                  myqtldata = myqtldata,
                  output_qtlplot = output_qtlplot,
                  output_qtllod = output_qtllod,
                  dataname = "unphased"
){
  
  pdf(output_datasummary, 15, 10)
  plot(myqtldata)
  dev.off()
  
  data_jitter <- jittermap(myqtldata)
  
  crossInt1 <- calc.genoprob(data_jitter, step = 1, error.prob = 0.001)
  
  myqtlout <- NULL
  
  pdf(output_qtlplot, 10, 5)
  for(i in 1:length(phenames(data_jitter))){
    
    print(paste("phenotype", i, phenames(data_jitter)[i], sep = " "))
    
    Em<- scanone(crossInt1, method="em", pheno.col = i)
    
    Cim<-cim(crossInt1, n.marcovar = 5, window=10, pheno.col = i)
    
    operm <- scanone(crossInt1, n.perm=200, verbose=FALSE, pheno.col = i)
    
    plot(Em, 
         ylab="LOD score", 
         col = "red", 
         ylim=c(0, max(Em$lod, Cim$lod)+5),
         main = paste(dataname, phenames(data_jitter)[i], sep = "_"))
    
    # plot the LOD score
    plot(Cim, 
         ylab="LOD score", 
         col="blue", 
         ylim=c(0,max(Em$lod, Cim$lod)+5),
         add = TRUE)
    
    legend("topright", c("cim", "em"), col = c("blue", "red"), lty = 1, cex = 0.5, title = "QTL mapping methods")
    
    # use the 95% percentile of the permutation test as the suggested threshold
    abline(h= round(as.numeric(sort(operm[1:200]))[190], 2) ,col="black", lty=3, lwd = 2)
    
    print(round(as.numeric(sort(operm[1:200]))[190], 2))
    
    Cim$em <- Em$lod
    Cim$loc <- rownames(Cim)
    Cim$trait <- rep(phenames(data_jitter)[i], nrow(Cim))
    Cim$threshold <- rep(round(as.numeric(sort(operm[1:200]))[190], 2), nrow(Cim))
    
    colnames(Cim) <- c("chr", "pos", "cim_lod", "em_lod", "loc", "trait", "threshold")
    
    myqtlout <- rbind.data.frame(myqtlout, Cim)
    
  }
  dev.off()
  
  write.table(myqtlout, output_qtllod, quote = F, row.names = F, col.names = T, sep = "\t")
  
}



# myqtl(output_datasummary = "unphased/Rabiosa_unphased_pheno_all.pdf",
#       myqtldata = data_unphased,
#       output_qtlplot = "unphased/Rabiosa_unphased_qtl_plot_all.pdf",
#       output_qtllod = "unphased/Rabiosa_unphased_qtl_lod_all.txt", 
#       dataname = "Unphased")
# 
# myqtl(output_datasummary = "hap1/Rabiosa_hap1_pheno_all.pdf",
#       myqtldata = data_hap1,
#       output_qtlplot = "hap1/Rabiosa_hap1_qtl_plot_all.pdf",
#       output_qtllod = "hap1/Rabiosa_hap1_qtl_lod_all.txt",
#       dataname = "Hap1")
# 
# myqtl(output_datasummary = "hap2/Rabiosa_hap2_pheno_all.pdf",
#       myqtldata = data_hap2,
#       output_qtlplot = "hap2/Rabiosa_hap2_qtl_plot_all.pdf",
#       output_qtllod = "hap2/Rabiosa_hap2_qtl_lod_all.txt",
#       dataname = "Hap2")
# 
# myqtl(output_datasummary = "vg/Rabiosa_vg_pheno_all.pdf",
#       myqtldata = data_vg,
#       output_qtlplot = "vg/Rabiosa_vg_qtl_plot_all.pdf",
#       output_qtllod = "vg/Rabiosa_vg_qtl_lod_all.txt",
#       dataname = "Graph_ref")


# myqtl(output_datasummary = "unphased_consensus/Rabiosa_unphased_consensus_pheno_all.pdf",
#       myqtldata = data_unphased_consensus,
#       output_qtlplot = "unphased_consensus/Rabiosa_unphased_consensus_qtl_plot_all.pdf",
#       output_qtllod = "unphased_consensus/Rabiosa_unphased_consensus_qtl_lod_all.txt", 
#       dataname = "Unphased_consensus")

myqtl(output_datasummary = "unphased_sikem_only/Rabiosa_unphased_sikem_only_pheno_all.pdf",
      myqtldata = data_unphased_sikem,
      output_qtlplot = "unphased_sikem_only/Rabiosa_unphased_sikem_only_qtl_plot_all.pdf",
      output_qtllod = "unphased_sikem_only/Rabiosa_unphased_sikem_only_qtl_lod_all.txt", 
      dataname = "Unphased_sikem_only")
