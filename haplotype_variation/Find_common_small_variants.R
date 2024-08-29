# find the common SNPs, INDELs between BCFtools and Syri
myargs <- commandArgs(trailingOnly=TRUE)

# myargs[1] = BCFSNP
# myargs[2] = SYRISNP
# myargs[3] = BCFINDEL
# myargs[4] = SYRIINDEL

BCFSNP <- read.table(myargs[1], header = F, stringsAsFactors = F)
SYRISNP <- read.table(myargs[2], header = F, stringsAsFactors = F)
BCFINDEL <- read.table(myargs[3], header = F, stringsAsFactors = F)
SYRIINDEL <- read.table(myargs[4], header = F, stringsAsFactors = F)

mychr <- sort(unique(BCFSNP$V1))

SNP_common <- NULL
INDEL_common <- NULL

for(x in mychr){
  
  BCFSNP_chr <- BCFSNP[BCFSNP$V1 == x, ]
  SYRISNP_chr <- SYRISNP[SYRISNP$V1 == x, ]
  
  BCFINDEL_chr <- BCFINDEL[BCFINDEL$V1 == x, ]
  SYRIINDEL_chr <- SYRIINDEL[SYRIINDEL$V1 == x, ]
  
  SNP_common <- rbind.data.frame(SNP_common, SYRISNP_chr[SYRISNP_chr$V2 %in% BCFSNP_chr$V2, ])
  INDEL_common <- rbind.data.frame(INDEL_common, SYRIINDEL_chr[SYRIINDEL_chr$V2 %in% BCFINDEL_chr$V2, ])
  
}

write.table(SNP_common, "SNP_common.txt", quote = F, sep = "\t", row.names = F, col.names = F)
write.table(INDEL_common, "INDEL_common.txt", quote = F, sep = "\t", row.names = F, col.names = F)


