# set working directory
setwd("P:/Yutangchen/Rabiosa_canu/winnowmap")

myblock <- read.table("phased.list", header = F, stringsAsFactors = F)
chrlength <- read.table("../reference/211209_Rabiosa_new_pseudomolecules_v1.fasta.fai", header = F, stringsAsFactors = F)
chrlength <- chrlength[, 2]

png('phase_block.png', width = 8, height = 5, units = 'in', res = 480)
par(mar=c(5, 5, 5, 4))
plot(seq(0, (max(chrlength) + 50*10^6), 10^6),
     seq(0, (max(chrlength) + 50*10^6), 10^6),
     ylim = c(0, 30), 
     ann = F, axes = F, type = "n")

for(i in 1:7){
  
  block <- myblock[myblock$V2 == paste("chr", i, sep = ""), ]
  block_largest <- block[which.max(block$V6), ]
  block_small <- block[-which.max(block$V6), ]
  
  rect(xleft = 1, ybottom = i*4 - 0.5, xright = block_largest[1, 5], ytop = i*4, col = "red2", border = NA)
  
  for(j in 1:dim(block_small)[1]){
    rect(xleft = block_small[j, 4], 
         ybottom = i*4 + 0.5, 
         xright = block_small[j, 5], 
         ytop = i*4 + 0.8, 
         col = "black")
  }
  
  segments(x0 = chrlength[i]+5, 
           y0 = i*4 - 0.8, 
           x1 = chrlength[i]+5, 
           y1 = i*4 + 0.3, 
           col = "blue", 
           lwd = 2)
  
  mtext(text = paste("Chr", i, sep = ""),
        side = 2,
        at = i*4, 
        las = 2,
        cex = 0.8)
}


axis(side = 1, 
     at = seq(0, max(chrlength) + 50*10^6, 50*10^6), 
     labels = c(0, paste(seq(50*10^6, max(chrlength) + 50*10^6, 50*10^6)/10^6, "", sep = "")), 
     cex.axis = 0.8, 
     line = -1)

legend("topright", 
       legend = c("Largest phase block", "Small phase block"),
       fill = c("red", "black"),
       cex = 0.6, bty = "n",  
       x.intersp = 0.5, y.intersp = 1)

title(xlab = 'Base pair positioin (Mb)', line = 1)
dev.off()

##########################################################################################################
# draw how much data were tagged
tag_stats <- read.table("tag_reads.stats", header = T, stringsAsFactors = F)

tag_length <- ceiling(as.numeric(gsub(",", "",tag_stats$sum_len))/10^9)

png('total_length_reads_assigned.png', width = 5, height = 8, units = 'in', res = 480)
par(mar=c(6, 6, 4, 6))
plot(1:sum(tag_length), 1:sum(tag_length), axes = F, ann = F, type = "n", xlim = c(0, 20))

rect(xleft = 1, 
     ybottom = 0,
     xright = 1 + 3,
     ytop = tag_length[1],
     col = "indianred2")

rect(xleft = 1, 
     ybottom = tag_length[1],
     xright = 1 + 3,
     ytop = sum(tag_length[1:2]),
     col = "lightblue2")

rect(xleft = 7, 
     ybottom = 0,
     xright = 7 + 3,
     ytop = tag_length[3],
     col = "gray80")

rect(xleft = 13, 
     ybottom = 0,
     xright = 13 + 3,
     ytop = tag_length[4],
     col = "gray10")

text(x = c(2, 8, 14), y = c(0, 0, 0), srt = 70, adj = c(1, 1), xpd = T, offset = 1,
     labels = c("Tagged", "Untagged", "Unaligned"), cex = 0.8)

axis(side = 2, 
     at = seq(0, sum(tag_length), 20), 
     labels = c(0, paste(seq(20, sum(tag_length), 20), "", sep = "")), 
     las = 2, cex.axis = 0.8)

title(ylab = "Total length of reads assigned (Gb)", cex = 0.8)

legend("topright", legend = c("Haplotype1", "Haplotype2", "Untagged", "Unaligned"),
       fill = c("indianred2", "lightblue3", "gray80", "gray10"),
       x.intersp = 0.5, y.intersp = 0.8, bty = "n",
       yjust = 0.5, xpd=T, inset=c(-0.2,0.4))
dev.off()