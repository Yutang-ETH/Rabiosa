# compare genetic linkage maps between unphased and vg graph reference

# 21.02.2023

# need fasta.fai, QTL mapping results (LOD score)

# set working directory
setwd("P:/Yutangchen/GBS_snakemake")

# load fasta.fai
myfai <- read.table("rabiosa_f/genome/assembly.fa.fai", header = F, stringsAsFactors = F)
myfai <- myfai[1:7, 1:2]

# load mapping results
qtl_unphased <- read.table("stemrustqtl/unphased_rabiosa_only/Rabiosa_unphased_qtl_lod_all.txt", header = T, stringsAsFactors = F)
qtl_vg <- read.table("stemrustqtl/vg/Rabiosa_vg_qtl_lod_all.txt", header = T, stringsAsFactors = F)

# select rows with "_", "chr1_1111"
qtl_unphased <- qtl_unphased[grep("_", qtl_unphased$loc), ]
qtl_vg <- qtl_vg[grep("_", qtl_vg$loc), ]

# only focus on stemrust data
qtl_unphased_sr <- qtl_unphased[qtl_unphased$trait == "emmeanSR", ]
qtl_vg_sr <- qtl_vg[qtl_vg$trait == "emmeanSR", ]

# get the base pair position
qtl_unphased_sr$bp <- as.numeric(sub('.*_', '', qtl_unphased_sr$loc))
qtl_vg_sr$bp <- as.numeric(sub('.*_', '', qtl_vg_sr$loc))

# start plotting
mygap <- 30*10^6
myx <- cumsum(as.numeric(myfai$V2) + mygap)

mybegin <- c(0, myx[1:6]) + mygap
myend <- myx - mygap + mygap

qtl_unphased_sr$newbp <- qtl_unphased_sr$bp + rep(mybegin, summary(as.factor(qtl_unphased_sr$chr))) 
qtl_vg_sr$newbp <- qtl_vg_sr$bp + rep(mybegin, summary(as.factor(qtl_vg_sr$chr)))

par(mar = c(5, 6, 2, 2))

plot(x = seq(0, max(myx), 10^6),
     y = seq(0, max(myx), 10^6), 
     ylim = c(0, 60), 
     type = "n", axes = F, ann = F)

# draw the map
segments(x0 = mybegin, x1 = myend, y0 = 10, y1 = 10, col = "gray40", lwd = 15)
segments(x0 = mybegin, x1 = myend, y0 = 20, y1 = 20, col = "gray40", lwd = 15)


# add marker postion to the map
segments(x0 = qtl_unphased_sr$newbp, x1 = qtl_unphased_sr$newbp, y0 = 9, y1 = 11, lwd = 0.2, col = "darkorange")
segments(x0 = qtl_vg_sr$newbp, x1 = qtl_vg_sr$newbp, y0 = 19, y1 = 21, lwd = 0.2, col = "darkorange")

# find common SNPs between two map
common_marker <- intersect(qtl_unphased_sr$newbp, qtl_vg_sr$newbp)
segments(x0 = common_marker, x1 = common_marker, y0 = 12, y1 = 18, col = "black")

# plot QTL mappinp results
for(i in 1:7){
  
  points(qtl_unphased_sr[qtl_unphased_sr$chr == i, 9], qtl_unphased_sr[qtl_unphased_sr$chr == i, 3] + 25, type = "l", col = "red", lwd = 2)
  points(qtl_vg_sr[qtl_vg_sr$chr == i, 9], qtl_vg_sr[qtl_vg_sr$chr == i, 3] + 25, type = "l", col = "blue", lwd = 2)
  
}

# add significant threshold
segments(x0 = mygap, x1 = max(myx) - mygap + mygap, y0 = 3.08 + 25, y1 = 3.08 + 25, col = "black", lty = 2, lwd = 1)

# add x axis coordinate
segments(x0 = mybegin, x1 = myend, y0 = 5, y1 = 5)

for(i in 1:7){
  
  segments(x0 = seq(mybegin[i], myend[i], 50*10^6), x1 = seq(mybegin[i], myend[i], 50*10^6), y0 = 4.5, y1 = 5)
  text(x = seq(mybegin[i], myend[i], 50*10^6), y = 4.5, 
       labels = (seq(mybegin[i], myend[i], 50*10^6) - mybegin[i])/10^6,
       cex = 1, pos = 1)
  text(x = (mybegin[i] + myend[i])/2, y = 1,
       labels = i, 
       cex = 1.5, pos = 1, xpd = NA)
  
}

# add yaxis labels
text(x = mygap, y = 10, labels = "Single\nreference", pos = 2, cex = 1.2, xpd = NA)
text(x = mygap, y = 20, labels = "Graph-based\nreference", pos = 2, cex = 1.2, xpd = NA)

# add xaxis labels
text(x = (mygap + max(myend))/2, y = 0, labels = "Chromosome (Mb)", pos = 1, cex = 1.5, xpd = NA, offset = 1.5)

# add yaxis lod score
segments(x0 = 20*10^6, x1 = 20*10^6, y0 = 25, y1 = 25 + 30, col = "black", xpd = NA)
segments(x0 = 15*10^6, x1 = 20*10^6, y0 = seq(25, 25 + 30, 5), y1 = seq(25, 25 + 30, 5), xpd = NA)
text(x = 15*10^6, y = seq(25, 25 + 30, 5), labels = seq(25, 25 + 30, 5) - 25, cex = 1.2, pos = 2, xpd = T)
text(x = 0, y = (25 + 25 + 30)/2, labels = "LOD", cex = 1.5, srt = 90, pos = 2, offset = 2.5, xpd = NA)

# add legened to the QTL mapping plot
text(x = mygap, y = 25 + 25, labels = "Graph-based reference", col = "blue", cex = 1.2, pos = 4)
text(x = mygap, y = 25 + 22, labels = "Single reference", col = "red", cex = 1.2, pos = 4)

# annotate the peak
# get the coordinate of the highest dot in chr7
qtl_unphased_chr7 <- qtl_unphased_sr[qtl_unphased_sr$chr == 7, ]
unphased_highest_x <- qtl_unphased_chr7[which.max(qtl_unphased_chr7$cim_lod), 9]
unphased_highest_y <- max(qtl_unphased_chr7$cim_lod) + 25
qtl_vg_chr7 <- qtl_vg_sr[qtl_vg_sr$chr == 7, ]
vg_highest_x <- qtl_vg_chr7[which.max(qtl_vg_chr7$cim_lod), 9]
vg_highest_y <- max(qtl_vg_chr7$cim_lod) + 25

segments(x0 = unphased_highest_x - 5*mygap, x1 = unphased_highest_x - 0.5*mygap, y0 = unphased_highest_y + 2, y1 = unphased_highest_y)
text(x = unphased_highest_x - 5*mygap, y = unphased_highest_y + 2, labels = qtl_unphased_chr7[which.max(qtl_unphased_chr7$cim_lod), 8], cex = 1.2, pos = 2, offset = 0.2)

segments(x0 = vg_highest_x - 5*mygap, x1 = vg_highest_x - 0.2*mygap, y0 = vg_highest_y + 2, y1 = vg_highest_y)
text(x = vg_highest_x - 5*mygap, y = vg_highest_y + 2, labels = qtl_vg_chr7[which.max(qtl_vg_chr7$cim_lod), 8], cex = 1.2, pos = 2, offset = 0.2)

# get the coordinate of the highest dot in chr2 
qtl_unphased_chr2 <- qtl_unphased_sr[qtl_unphased_sr$chr == 2, ]
unphased_highest_x <- qtl_unphased_chr2[which.max(qtl_unphased_chr2$cim_lod), 9]
unphased_highest_y <- max(qtl_unphased_chr2$cim_lod) + 25
qtl_vg_chr2 <- qtl_vg_sr[qtl_vg_sr$chr == 2, ]
vg_highest_x <- qtl_vg_chr2[which.max(qtl_vg_chr2$cim_lod), 9]
vg_highest_y <- max(qtl_vg_chr2$cim_lod) + 25

segments(x0 = unphased_highest_x + 3*mygap, x1 = unphased_highest_x, y0 = unphased_highest_y + 10, y1 = unphased_highest_y)
text(x = unphased_highest_x + 3*mygap, y = unphased_highest_y + 10, labels = qtl_unphased_chr2[which.max(qtl_unphased_chr2$cim_lod), 8], cex = 1.2, pos = 4, offset = 0.2)

segments(x0 = vg_highest_x + 3*mygap, x1 = vg_highest_x, y0 = vg_highest_y + 5, y1 = vg_highest_y)
text(x = vg_highest_x + 3*mygap, y = vg_highest_y + 5, labels = qtl_vg_chr2[which.max(qtl_vg_chr2$cim_lod), 8], cex = 1.2, pos = 4, offset = 0.2)
