# install package and load the library
install.packages("circlize")
library(circlize)

# set the working directory
setwd("P:/Yutangchen/Rabiosa_assembly/canu/circos_plot")

# import chromosome information and make the required input file for circlize
myfai <- read.table("../rabiosa_assembly_final/rabiosa_dip_chr.fa.fai", header = F, stringsAsFactors = F)

mychrinfo <- cbind.data.frame(myfai$V1, rep(0, nrow(myfai)), myfai$V2)
colnames(mychrinfo) <- c("chr", "start", "end")
mychrinfo$chr <- factor(mychrinfo$chr, levels = mychrinfo$chr)
head(mychrinfo)

# import gene and repeat 
HC_hap_0 <- read.table("hap_0/EVM.hc.bed", header = F, stringsAsFactors = F)
LC_hap_0 <- read.table("hap_0/EVM.lc.bed", header = F, stringsAsFactors = F)
REPEAT_hap_0 <- read.table("hap_0/repeat.bed", header = F, stringsAsFactors = F)

HC_hap_1 <- read.table("hap_1/EVM.hc.bed", header = F, stringsAsFactors = F)
LC_hap_1 <- read.table("hap_1/EVM.lc.bed", header = F, stringsAsFactors = F)
REPEAT_hap_1 <- read.table("hap_1/repeat.bed", header = F, stringsAsFactors = F)

HC_hap_2 <- read.table("hap_2/EVM.hc.bed", header = F, stringsAsFactors = F)
LC_hap_2 <- read.table("hap_2/EVM.lc.bed", header = F, stringsAsFactors = F)
REPEAT_hap_2 <- read.table("hap_2/repeat.bed", header = F, stringsAsFactors = F)

HC <- rbind.data.frame(HC_hap_0[order(HC_hap_0$V1), ], HC_hap_1[order(HC_hap_1$V1), ], HC_hap_2[order(HC_hap_2$V1), ])
HC$V4 <- floor((HC$V3 + HC$V2)/2)
HC$V1 <- paste(c(rep("rabiosa_0", nrow(HC_hap_0)), 
                 rep("rabiosa_1", nrow(HC_hap_1)),
                 rep("rabiosa_2", nrow(HC_hap_2))), HC$V1, sep = "_")

LC <- rbind.data.frame(LC_hap_0[order(LC_hap_0$V1), ], LC_hap_1[order(LC_hap_1$V1), ], LC_hap_2[order(LC_hap_2$V1), ])
LC$V4 <- floor((LC$V3 + LC$V2)/2)
LC$V1 <- paste(c(rep("rabiosa_0", nrow(LC_hap_0)), 
                 rep("rabiosa_1", nrow(LC_hap_1)),
                 rep("rabiosa_2", nrow(LC_hap_2))), LC$V1, sep = "_")

REPEAT <- rbind.data.frame(REPEAT_hap_0[order(REPEAT_hap_0$V1), ], REPEAT_hap_1[order(REPEAT_hap_1$V1), ], REPEAT_hap_2[order(REPEAT_hap_2$V1), ])
REPEAT$V4 <- floor((REPEAT$V3 + REPEAT$V2)/2)
REPEAT$V5 <- floor(REPEAT$V3 - REPEAT$V2)
REPEAT$V1 <- paste(c(rep("rabiosa_0", nrow(REPEAT_hap_0)), 
                 rep("rabiosa_1", nrow(REPEAT_hap_1)),
                 rep("rabiosa_2", nrow(REPEAT_hap_2))), REPEAT$V1, sep = "_")

# import genetic map
spm_unphased <- read.table("P:/Yutangchen/GBS_snakemake/rabiosa_f/lepmap/scaffold_position_map.txt", header = F, stringsAsFactors = F)
spm_hap1 <- read.table("P:/Yutangchen/GBS_snakemake/hap1_f/lepmap/scaffold_position_map.txt", header = F, stringsAsFactors = F)
spm_hap2 <- read.table("P:/Yutangchen/GBS_snakemake/hap2_f/lepmap/scaffold_position_map.txt", header = F, stringsAsFactors = F)

spm_unphased$V1 <- paste("rabiosa_0", spm_unphased$V1, sep = "_")
spm_hap1$V1 <- paste("rabiosa_1", spm_hap1$V1, sep = "_")
spm_hap2$V1 <- paste("rabiosa_2", spm_hap2$V1, sep = "_")

spm <- rbind.data.frame(spm_unphased, spm_hap1, spm_hap2)

# GC content
GC_unphased <- read.table("hap_0/rabiosa_unphased_chr.gc", header = F, stringsAsFactors = F)
GC_hap1 <- read.table("hap_1/rabiosa_hap1_chr.gc", header = F, stringsAsFactors = F)
GC_hap2 <- read.table("hap_2/rabiosa_hap2_chr.gc", header = F, stringsAsFactors = F)
GC_unphased$V1 <- sub("_sliding:.*", "", GC_unphased$V1)
GC_hap1$V1 <- sub("_sliding:.*", "", GC_hap1$V1)
GC_hap2$V1 <- sub("_sliding:.*", "", GC_hap2$V1)
GC <- rbind.data.frame(GC_unphased, GC_hap1, GC_hap2)

#-------------------------------------------------------------------------------------------------------#

# make the circos plot

png('rabiosa_circos_new.png', width = 10, height = 10, units = 'in', res = 480)

# add gap between genomes
circos.par(gap.after = c(rep(1, 6), 1, rep(1, 6), 1, rep(1, 6), 15), 
           track.margin = c(0.01, 0.01),
           start.degree = 90)

# initialize the circos plot
circos.initialize(mychrinfo$chr, xlim = mychrinfo[, 2:3])

# first track, highlight each assembly
circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(0.3),
             bg.col = c(rep("gray80", 7), rep(adjustcolor("red", alpha.f = 0.5), 7), rep(adjustcolor("blue", alpha.f = 0.5), 7)),
             bg.border = "black" ) # c(rep("gray60", 7), rep("tomato", 7), rep("lightblue", 7)))

for(i in 1:nrow(mychrinfo)){
  
  circos.text((mychrinfo[i, 3] - mychrinfo[i, 2])/2, 1 + mm_y(4), 
              labels = gsub(".*chr", "", levels(mychrinfo$chr)[i]), 
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              facing = "bending.inside",
              cex = 1)
  
  circos.axis(h = "top", 
              major.at = seq(0, mychrinfo[i, 3] + 50*10^6, by = 50*10^6),
              labels = seq(0, floor((mychrinfo[i, 3] + 50*10^6)/10^6), by = 50),
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              labels.cex = 0.5,
              lwd = 1,
              minor.ticks = 0,
              major.tick.length = mm_y(0.2),
              direction = "outside",
              labels.facing = "outside",
              labels.pos.adjust = F)
              # col = c(rep("gray80", 7), rep("tomato", 7), rep("lightblue", 7))[i])

}

# second track, genetic linkage map vs physical map
circos.track(sector = mychrinfo$chr, ylim = c(0, floor(max(spm$V4))), cell.padding = c(0, 0, 0, 0), track.height = mm_h(5))

# make the plot
for(i in 1:nrow(mychrinfo)){
  
  circos.points(x = spm[spm$V1 == mychrinfo[i, 1], 2], 
                y = spm[spm$V1 == mychrinfo[i, 1], 4],
                sector.index = levels(mychrinfo$chr)[i],
                track.index = 2,
                pch = 1,
                cex = 0.1,)
}

circos.yaxis(
  side = "left",
  at = seq(0, 150, 30),
  labels = seq(0, 150, 30),
  labels.cex = 0.5,
  track.index = 2,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 2),
  sector.index = levels(mychrinfo$chr)[1])

# third track, distribution of HC gene, heat map or hist gram

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), cell.padding = c(0, 0, 0, 0), track.height = mm_h(4), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count how many genes in every window
  mynumber <- c()
  for(j in 1:(length(mywindow)-1)){
    mynumber <- c(mynumber, sum(HC[HC$V1 == mychrinfo[i, 1], 4] >= mywindow[j] & HC[HC$V1 == mychrinfo[i, 1], 4] < mywindow[j+1]))
  }
  
  # plot polygon
  circos.polygon(x = c(0, mywindow), 
                 y = c(0, mynumber, 0),
                 border = "red", # "gray30",
                 col = adjustcolor("red", alpha.f = 0.4), # "seagreen4",
                 track.index = 3,
                 sector.index = levels(mychrinfo$chr)[i])

}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 3),
  sector.index = levels(mychrinfo$chr)[1])

# fourth track, distribution of LC gene

circos.track(sector = mychrinfo$chr, ylim = c(0, 100), cell.padding = c(0, 0, 0, 0), track.height = mm_h(4), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count how many genes in every window
  mynumber <- c()
  for(j in 1:(length(mywindow)-1)){
    mynumber <- c(mynumber, sum(LC[LC$V1 == mychrinfo[i, 1], 4] >= mywindow[j] & LC[LC$V1 == mychrinfo[i, 1], 4] < mywindow[j+1]))
  }
  
  # plot polygon
  circos.polygon(x = c(0, mywindow), 
                 y = c(0, mynumber, 0),
                 border = "blue", # "gray30",
                 col = adjustcolor("blue", alpha.f = 0.4), # "greenyellow",
                 track.index = 4,
                 sector.index = levels(mychrinfo$chr)[i])
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 100, 25),
  labels = seq(0, 100, 25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 4),
  sector.index = levels(mychrinfo$chr)[1])

# fifth track, distribution of repeat

circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(3), bg.border = NA)

# use a line to show the abundance of repeats

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  mydf <- cbind.data.frame(mywindow[1:(length(mywindow)-1)], mywindow[2:length(mywindow)])
  names(mydf) <- c("win1", "win2")
  
  # count how many repeats in every window
  myfun <- function(myrepeat = myrepeat, x){
    return(round(sum(myrepeat[myrepeat[ ,4] >= x[1] & myrepeat[ ,4] < x[2], 5])/10^6, 2))
  }
  
  mypercent <- unlist(apply(mydf, 1, myfun, myrepeat = REPEAT[REPEAT$V1 == mychrinfo[i, 1], ]))
  
  # plot polygon
  circos.polygon(x = c(0, mywindow), 
                 y = c(0, mypercent, 0),
                 border = "green4", # "gray30",
                 col = adjustcolor("green4", alpha.f = 0.4), # "steelblue4",
                 track.index = 5,
                 sector.index = levels(mychrinfo$chr)[i])
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 1, 0.25),
  labels = seq(0, 1, 0.25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 5),
  sector.index = levels(mychrinfo$chr)[1])

# six track, distribution of GC contene

circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(3), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # 10^6, window size 1 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^6)
  
  # count how many GC in every window
  mypercent <- round(GC[GC$V1 == mychrinfo[i, 1], 2]/100, 4)
  
  # plot polygon
  # circos.polygon(x = c(0, mywindow, max(mywindow)), 
  #                y = c(0, mypercent, 0),
  #                border = "orange4", # "gray30",
  #                col = adjustcolor("orange4", alpha.f = 0.4), # "steelblue4",
  #                track.index = 6,
  #                sector.index = levels(mychrinfo$chr)[i])
  
  # plot line 
  circos.lines(x = mywindow-5.0e+05,
               y = mypercent,
               sector.index = levels(mychrinfo$chr)[i],
               col = "orange4",
               type = "s")
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0, 1, 0.25),
  labels = seq(0, 1, 0.25),
  labels.cex = 0.5,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 6),
  sector.index = levels(mychrinfo$chr)[1])

circos.clear()

# add labels to each track
text(0, 1, labels = "A", pos = 2, offset = 1.5)
text(0, 0.9, labels = "B", pos = 2, offset = 1.5)
text(0, 0.8, labels = "C", pos = 2, offset = 1.5)
text(0, 0.7, labels = "D", pos = 2, offset = 1.5)
text(0, 0.6, labels = "E", pos = 2, offset = 1.5)
text(0, 0.5, labels = "F", pos = 2, offset = 1.5)

text(-0.2, 0.2, labels = "A: assembly and chromosome", cex = 0.6, font = 2, pos = 4)
text(-0.2, 0.1, labels = "B: genetic map position vs chromosome base position", cex = 0.6, font = 2, pos = 4)
text(-0.2, 0, labels = "C: high confidence gene density", cex = 0.6, font = 2, pos = 4)
text(-0.2, -0.1, labels = "D: low confidence gene density", cex = 0.6, font = 2, pos = 4)
text(-0.2, -0.2, labels = "E: TE density", cex = 0.6, font = 2, pos = 4)
text(-0.2, -0.3, labels = "F: GC content", cex = 0.6, font = 2, pos = 4)
dev.off()

