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
spm_unphased <- read.table("P:/Yutangchen/GBS_snakemake/rabiosa_f/lepmap_consensus_map/scaffold_position_map.txt", header = F, stringsAsFactors = F)
spm_hap1 <- read.table("P:/Yutangchen/GBS_snakemake/hap1_f/lepmap/scaffold_position_map.txt", header = F, stringsAsFactors = F)
spm_hap2 <- read.table("P:/Yutangchen/GBS_snakemake/hap2_f/lepmap/scaffold_position_map.txt", header = F, stringsAsFactors = F)

spm_unphased$V1 <- paste("rabiosa_0", spm_unphased$V1, sep = "_")
spm_unphased <- spm_unphased[!duplicated(spm_unphased$V4), ]
spm_hap1$V1 <- paste("rabiosa_1", spm_hap1$V1, sep = "_")
spm_hap1 <- spm_hap1[!duplicated(spm_hap1$V4), ]
spm_hap2$V1 <- paste("rabiosa_2", spm_hap2$V1, sep = "_")
spm_hap2 <- spm_hap2[!duplicated(spm_hap2$V4), ]

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

# add gap between genomes
circos.par(gap.after = c(rep(1, 6), 1, rep(1, 6), 1, rep(1, 6), 15), 
           track.margin = c(0.01, 0.01),
           start.degree = 90)

# initialize the circos plot
circos.initialize(mychrinfo$chr, xlim = mychrinfo[, 2:3])

# first track, highlight each assembly
circos.track(sector = mychrinfo$chr, ylim = c(0, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(1),
             bg.col = c(rep("gray60", 7), rep(adjustcolor("tomato", alpha.f = 0.5), 7), rep(adjustcolor("lightblue", alpha.f = 0.5), 7)),
             bg.border = c(rep("gray60", 7), rep("tomato", 7), rep("lightblue", 7))) # c(rep("gray60", 7), rep("tomato", 7), rep("lightblue", 7)))

for(i in 1:nrow(mychrinfo)){
  
  circos.text((mychrinfo[i, 3] - mychrinfo[i, 2])/2, 3 + mm_y(4), 
              labels = gsub(".*chr", "", levels(mychrinfo$chr)[i]), 
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              facing = "bending.inside",
              cex = 2)
  
  circos.axis(h = "top", 
              major.at = seq(0, mychrinfo[i, 3] + 50*10^6, by = 50*10^6),
              labels = seq(0, floor((mychrinfo[i, 3] + 50*10^6)/10^6), by = 50),
              sector.index = levels(mychrinfo$chr)[i],
              track.index = 1,
              labels.cex = 0.7,
              lwd = 1,
              minor.ticks = 0,
              major.tick.length = mm_y(0.2),
              direction = "outside",
              labels.facing = "outside",
              labels.pos.adjust = F)
  # col = c(rep("gray80", 7), rep("tomato", 7), rep("lightblue", 7))[i])
  
}

# second track, genetic linkage map vs physical map
circos.track(sector = mychrinfo$chr, ylim = c(0, floor(max(spm$V4))), cell.padding = c(0, 0, 0, 0), track.height = mm_h(8))

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
  at = seq(0, 150, 50),
  labels = seq(0, 150, 50),
  labels.cex = 0.8,
  track.index = 2,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 2),
  sector.index = levels(mychrinfo$chr)[1])

# set gap size between tracks
set_track_gap(mm_h(1))

# third track, distribution of HC gene, heat map or hist gram
hchist <- NULL
for(i in 1:nrow(mychrinfo)){
  
  # 10^7, window size 10 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^7)
  
  # count how many genes in every window
  mynumber <- c()
  for(j in 1:(length(mywindow)-1)){
    mynumber <- c(mynumber, sum(HC[HC$V1 == mychrinfo[i, 1], 4] >= mywindow[j] & HC[HC$V1 == mychrinfo[i, 1], 4] < mywindow[j+1]))
  }
  
  hcchr <- rep(mychrinfo$chr[i], length(mynumber))
  hc_chr_window_count <- cbind.data.frame(hcchr, mywindow[-1], mynumber)
  colnames(hc_chr_window_count) <- c('chr', 'window', 'number')
  hchist <- rbind.data.frame(hchist, hc_chr_window_count)
  
}

circos.track(sector = mychrinfo$chr, ylim = c(min(hchist$number), max(hchist$number)), cell.padding = c(0, 0, 0, 0), track.height = mm_h(8), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # plot polygon
  circos.polygon(x = c(head(hchist[hchist$chr == mychrinfo$chr[i], 2], 1), 
                       hchist[hchist$chr == mychrinfo$chr[i], 2], 
                       tail(hchist[hchist$chr == mychrinfo$chr[i], 2], 1)), 
                 y = c(min(hchist[hchist$chr == mychrinfo$chr[i], 3]), 
                       hchist[hchist$chr == mychrinfo$chr[i], 3], 
                       min(hchist[hchist$chr == mychrinfo$chr[i], 3])),
                 border = "white", # "gray30",
                 col = adjustcolor("seagreen4", alpha.f = 0.4), # "seagreen4",
                 track.index = 3,
                 sector.index = levels(mychrinfo$chr)[i])
  
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(min(hchist$number), max(hchist$number), 100),
  labels = seq(min(hchist$number), max(hchist$number), 100),
  labels.cex = 0.8,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 3),
  sector.index = levels(mychrinfo$chr)[1])

# set gap size between tracks
set_track_gap(mm_h(3))

# fourth track, distribution of LC gene
lchist <- NULL
for(i in 1:nrow(mychrinfo)){
  
  # 10^7, window size 10 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^7)
  
  # count how many genes in every window
  mynumber <- c()
  for(j in 1:(length(mywindow)-1)){
    mynumber <- c(mynumber, sum(LC[LC$V1 == mychrinfo[i, 1], 4] >= mywindow[j] & LC[LC$V1 == mychrinfo[i, 1], 4] < mywindow[j+1]))
  }
  
  lcchr <- rep(mychrinfo$chr[i], length(mynumber))
  lc_chr_window_count <- cbind.data.frame(lcchr, mywindow[-1], mynumber)
  colnames(lc_chr_window_count) <- c('chr', 'window', 'number')
  lchist <- rbind.data.frame(lchist, lc_chr_window_count)
  
}

circos.track(sector = mychrinfo$chr, ylim = c(min(lchist$number), max(lchist$number)), cell.padding = c(0, 0, 0, 0), track.height = mm_h(8), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  

  # plot polygon
  circos.polygon(x = c(head(lchist[lchist$chr == mychrinfo$chr[i], 2], 1), 
                       lchist[lchist$chr == mychrinfo$chr[i], 2], 
                       tail(lchist[lchist$chr == mychrinfo$chr[i], 2], 1)), 
                 y = c(min(lchist[lchist$chr == mychrinfo$chr[i], 3]), 
                       lchist[lchist$chr == mychrinfo$chr[i], 3], 
                       min(lchist[lchist$chr == mychrinfo$chr[i], 3])),
                 border = "white", # "gray30",
                 col = adjustcolor("greenyellow", alpha.f = 0.4), # "greenyellow",
                 track.index = 4,
                 sector.index = levels(mychrinfo$chr)[i])
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(min(lchist$number), max(lchist$number), 100),
  labels = seq(min(lchist$number), max(lchist$number), 100),
  labels.cex = 0.8,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 4),
  sector.index = levels(mychrinfo$chr)[1])

# set gap size between tracks
set_track_gap(mm_h(3))

# fifth track, distribution of repeat
rephist <- NULL
for(i in 1:nrow(mychrinfo)){
  
  # 10^7, window size 10 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^7)
  mydf <- cbind.data.frame(mywindow[1:(length(mywindow)-1)], mywindow[2:length(mywindow)])
  names(mydf) <- c("win1", "win2")
  
  # count how many repeats in every window
  myfun <- function(myrepeat = myrepeat, x){
    return(round(sum(myrepeat[myrepeat[ ,4] >= x[1] & myrepeat[ ,4] < x[2], 5])/10^7, 2))
  }
  
  mypercent <- unlist(apply(mydf, 1, myfun, myrepeat = REPEAT[REPEAT$V1 == mychrinfo[i, 1], ]))
  
  repchr <- rep(mychrinfo$chr[i], length(mypercent))
  rep_chr_window_count <- cbind.data.frame(repchr, mywindow[-1], mypercent)
  colnames(rep_chr_window_count) <- c('chr', 'window', 'percent')
  rephist <- rbind.data.frame(rephist, rep_chr_window_count)
  
}


circos.track(sector = mychrinfo$chr, ylim = c(0.5, 1), cell.padding = c(0, 0, 0, 0), track.height = mm_h(6), bg.border = NA)

# use a line to show the abundance of repeats

for(i in 1:nrow(mychrinfo)){
  
  # plot line 
  circos.lines(x = rephist[rephist$chr == mychrinfo$chr[i], 2],
               y = rephist[rephist$chr == mychrinfo$chr[i], 3],
               sector.index = levels(mychrinfo$chr)[i],
               col = "blue4",
               type = "l")
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(0.5, 1, 0.25),
  labels = seq(0.5, 1, 0.25),
  labels.cex = 0.8,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 5),
  sector.index = levels(mychrinfo$chr)[1])

# set gap size between tracks
set_track_gap(mm_h(3))

# six track, distribution of GC contente
gcwindow <- c()
for(i in 1:nrow(mychrinfo)){
  
  # 10^7, window size 10 Mb
  mywindow <- seq(mychrinfo$start[i], mychrinfo$end[i], by = 10^7)
  
  gcwindow <- c(gcwindow, mywindow)
  
}

GC$window <- gcwindow
colnames(GC) <- c('chr', 'GC', 'window')


circos.track(sector = mychrinfo$chr, ylim = c(40, 50), cell.padding = c(0, 0, 0, 0), track.height = mm_h(6), bg.border = NA)

for(i in 1:nrow(mychrinfo)){
  
  # plot line 
  circos.lines(x = GC[GC$chr == mychrinfo$chr[i], 3],
               y = GC[GC$chr == mychrinfo$chr[i], 2],
               sector.index = levels(mychrinfo$chr)[i],
               col = "orange4",
               type = "s")
}

# add axis
circos.yaxis(
  side = "left",
  at = seq(40, 50, 5),
  labels = seq(40, 50, 5),
  labels.cex = 0.8,
  tick.length = convert_x(0.1, "mm", levels(mychrinfo$chr)[1], 6),
  sector.index = levels(mychrinfo$chr)[1])

circos.clear()

# add labels to each track
text(0, 1, labels = "A", pos = 2, offset = 1.3, cex = 1.5)
text(0, 0.9, labels = "B", pos = 2, offset = 1.3, cex = 1.5)
text(0, 0.75, labels = "C", pos = 2, offset = 1.3, cex = 1.5)
text(0, 0.6, labels = "D", pos = 2, offset = 1.3, cex = 1.5)
text(0, 0.45, labels = "E", pos = 2, offset = 1.3, cex = 1.5)
text(0, 0.3, labels = "F", pos = 2, offset = 1.3, cex = 1.5)


