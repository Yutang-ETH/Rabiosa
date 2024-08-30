# make rf plot
# load Rqtl
library(qtl)

# import the genotype data from lepmap3
mytest <- read.table("genotypes_x.txt", header = F, sep = "\t", stringsAsFactors = F)

# make the input file for allmaps
spm <- read.table("scaffold_position_map.txt", header = F, stringsAsFactors = F)

# remove prefix in chr names
spm$V1 <- sub("rabiosa_0_", "", spm$V1)

# convert the format to csvr format for R/qtl
dim(mytest) # 100240 rows and 310 columns

# exclude the first 7 rows and first 4 columns
mygeno <- mytest[-c(1:6), -c(1:4)]

# in my genome, convert "1 2" and "2 1" to H, "1 1" to A and "2 2" to B

for (i in 1:dim(mygeno)[2]){
  mygeno[, i][mygeno[, i] == "1 2"] <- "H"
  mygeno[, i][mygeno[, i] == "2 1"] <- "H"
  mygeno[, i][mygeno[, i] == "1 1"] <- "A"
  mygeno[, i][mygeno[, i] == "2 2"] <- "B"
}

# clolumn combine marker name, linkage group, position and genotype
# mycsvr <- cbind.data.frame(paste(spm$V1, spm$V2, sep = "_"),
#                            sub("chr", "", spm$V1),
#                            spm$V3,
#                            mygeno)
# mycsvr <- type.convert(mycsvr, as.is = T)

# clolumn combine marker name, linkage group, position and genotype for Rabiosa, use spm$V4
mycsvr <- cbind.data.frame(paste(spm$V1, spm$V2, sep = "_"),
                           sub("chr", "", spm$V1),
                           spm$V4,
                           mygeno)
mycsvr <- type.convert(mycsvr, as.is = T)

# add one row of phenotype
mycsvr <- rbind.data.frame(c("pheno", "", "", rep(0, dim(mygeno)[2])), mycsvr)


# write it out as a comma delimited file
# write.table(mycsvr, "mycsvr.csv", sep = ",", row.names = F, col.names = F, quote = F)
write.table(mycsvr, "mycsvr_rabiosa.csv", sep = ",", row.names = F, col.names = F, quote = F)

# now import mycsvr with read.cross funtion in R/qtl
mydata <- read.cross(format = "csvr", file = "mycsvr_rabiosa.csv")

# png("rf_plot_selected.png", 1000, 1000)
# layout(matrix(1:9, 3, 3, byrow = T))
# for (i in 1:7){
#   plotRF(data_jitter, i, what = "rf")
# }
# dev.off()


png("rf_plot_vg_all.png")
plotRF(mydata, what = "rf")
dev.off()
