# make the synteny plot based on collinear blocks

# install.packages('yarrr')
# library(yarrr)

# set working directory
setwd("P:/Yutangchen/Rabiosa_compare/mysyri")

# import bed and anchor file
myfai <- read.table("Rabiosa_dip_chr.fa.fai", header = F, stringsAsFactors = F)
# mypaf <- read.table("Rabiosa_all_chr_nocg.paf", header = F, stringsAsFactors = F)
mysyn <- read.table("Rabiosa_variant_SYN.txt", header = F, stringsAsFactors = F)
mysyn$V1 <- gsub("chr", "rabiosa_1_chr", mysyn$V1)
mysyn$V4 <- gsub("chr", "rabiosa_2_chr", mysyn$V4)

mynotal <- read.table("Rabiosa_variant_NOTAL.txt", header = F, stringsAsFactors = F)
mynotal$V1 <- gsub("chr", "rabiosa_1_chr", mynotal$V1)
mynotal$V4 <- gsub("chr", "rabiosa_2_chr", mynotal$V4)

myinv <- read.table("Rabiosa_variant_INV.txt", header = F, stringsAsFactors = F)
myinv$V1 <- gsub("chr", "rabiosa_1_chr", myinv$V1)
myinv$V4 <- gsub("chr", "rabiosa_2_chr", myinv$V4)

mysnp <- read.table("SNP_common.txt", header = F, stringsAsFactors = F)
mysnp$V1 <- gsub("chr", "rabiosa_1_chr", mysnp$V1)

myindel <- read.table("INDEL_common.txt", header = F, stringsAsFactors = F)
myindel$V1 <- gsub("chr", "rabiosa_1_chr", myindel$V1)

mypav <- read.table("PAV_overlap.txt", header = F, stringsAsFactors = F)
mypav$V1 <- gsub("chr", "rabiosa_1_chr", mypav$V1)

mygap_hap1 <- read.table("hap1_gap.txt", header = F, stringsAsFactors = F)
mygap_hap1$V1 <- gsub("chr", "rabiosa_1_chr", mygap_hap1$V1)

mygap_hap2 <- read.table("hap2_gap.txt", header = F, stringsAsFactors = F)
mygap_hap2$V1 <- gsub("chr", "rabiosa_2_chr", mygap_hap2$V1)
#--------------------------------------------------------------------------------------------#

# change paf to 1-based coordinate and calculate blast-like alignment identity
# mypaf$V3 <- mypaf$V3 + 1
# mypaf$V8 <- mypaf$V8 + 1
# blast identity is V10 divided by V11
# mypaf$V12 <- round(mypaf$V10/mypaf$V11, 2)

# calculate indel size, the indel size should be the size of the ref genotype
myindel$V6 <- nchar(myindel$V4)

# select only common SV between syri and sniffles and calculate the length of SV for the reference
mypav <- mypav[mypav$V6 == "common", ]
mypav$V4 <- mypav$V3 - mypav$V2 + 1

###
mycount_fun <- function(x, myvar = myvar){
  return(length(myvar[myvar >= as.numeric(x[1]) & myvar < as.numeric(x[2])]))
}

Intersect_variant <- function(x, myvar = myvar){
  
  myvar$V2 <- as.numeric(myvar$V2)
  myvar$V3 <- as.numeric(myvar$V3)
  
  # the variant is contained in the window
  var1 <- myvar[myvar[, 2] >= as.numeric(x[1]) & myvar[, 2] < as.numeric(x[2]) & myvar[, 3] >= as.numeric(x[1]) & myvar[, 3] < as.numeric(x[2]), ]
  
  # left overlap window
  var2 <- myvar[myvar[, 2] < as.numeric(x[1]) & myvar[, 3] >= as.numeric(x[1]) & myvar[, 3] < as.numeric(x[2]), ]
  
  # window contained in the variant
  var3 <- myvar[myvar[, 2] < as.numeric(x[1]) & myvar[, 3] >= as.numeric(x[2]), ]
  
  # right overlap window
  var4 <- myvar[myvar[, 2] >= as.numeric(x[1]) & myvar[, 2] < as.numeric(x[2]) & myvar[, 3] >= as.numeric(x[2]), ]
  
  # 
  
  if(nrow(var1) >= 1){
    var1_overlap <- sum(var1[, 3] - var1[, 2] + 1)
  }else{
    var1_overlap <- 0
  }
  
  if(nrow(var2) >= 1){
    var2_overlap <- sum(var2[, 3] - as.numeric(x[1]) + 1)
  }else{
    var2_overlap <- 0
  }
  
  
  if(nrow(var3) >= 1){
    var3_overlap <- sum(as.numeric(x[2]) - as.numeric(x[1]))
  }else{
    var3_overlap <- 0
  }
  
  if(nrow(var4) >= 1){
    var4_overlap <- sum(as.numeric(x[2]) - var4[, 2])
  }else{
    var4_overlap <- 0
  }
  
  return(sum(var1_overlap, var2_overlap, var3_overlap, var4_overlap))
}

#--------------------------------------------------------------------------------------------#
# make a function to draw synteny map
Draw_synteny <- function(mysyn = mysyn, 
                         myinv = myinv, 
                         chr = chr,
                         chr2 = chr2,
                         mygap_hap1 = mygap_hap1,
                         mygap_hap2 = mygap_hap2,
                         hap1_length = hap1_length,
                         hap2_length = hap2_length,
                         mynotal = mynotal,
                         mysnp = mysnp,
                         myindel = myindel,
                         mypav = mypav,
                         mymar = c(6, 6, 3, 6),
                         mylable = T){
  
  # make synteny plot
  par(xpd = T, mar=mymar, cex.lab = 1.5)
  plot(x = seq(1, max(hap1_length, hap2_length), 10^6), 
       y = seq(1, max(hap1_length, hap2_length), 10^6), 
       axes = F,
       xlim = c(0, max(hap1_length, hap2_length)), 
       ylim = c(0, 55), 
       type = "n",
       xlab = paste(sub(".*chr", "Chr", chr), '(Mb)', sep = ' '), ylab = "")
  
  # syntenic region
  mysyn_chr <- mysyn[mysyn$V1 == chr, ]
  mysyn_chr$V2 <- mysyn_chr$V2 
  mysyn_chr$V3 <- mysyn_chr$V3 
  mysyn_chr$V5 <- mysyn_chr$V5 
  mysyn_chr$V6 <- mysyn_chr$V6 
  
  for(j in 1:nrow(mysyn_chr)){
    
    
    # draw segments from ref to query
    segments(x0 = mysyn_chr$V2[j], 
             y0 = 20,
             x1 = mysyn_chr$V5[j],
             y1 = 2,
             col = "lightgray")
    
    segments(x0 = mysyn_chr$V3[j], 
             y0 = 20,
             x1 = mysyn_chr$V6[j],
             y1 = 2,
             col = "lightgray")
    
    polygon(c(mysyn_chr$V2[j], mysyn_chr$V3[j], mysyn_chr$V6[j], mysyn_chr$V5[j]), 
            c(20, 20, 2, 2),
            col = "lightgray",
            border = NA)
    
  }
  
  # inverted region
  myinv_chr <- myinv[myinv$V1 == chr, ]
  myinv_chr$V2 <- myinv_chr$V2 
  myinv_chr$V3 <- myinv_chr$V3 
  myinv_chr$V5 <- myinv_chr$V5 
  myinv_chr$V6 <- myinv_chr$V6 
  
  for(j in 1:nrow(myinv_chr)){
    
    # draw segments from ref to query
    segments(x0 = myinv_chr$V2[j], 
             y0 = 20,
             x1 = myinv_chr$V6[j],
             y1 = 2,
             col = "gold")
    
    segments(x0 = myinv_chr$V3[j], 
             y0 = 20,
             x1 = myinv_chr$V5[j],
             y1 = 2,
             col = "gold")
    
    polygon(c(myinv_chr$V2[j], myinv_chr$V3[j], myinv_chr$V5[j], myinv_chr$V6[j]), 
            c(20, 20, 2, 2),
            col = "gold",
            border = NA)
    
  }
  
  
  # make the chr using a rect, ref frist then query
  segments(x0 = 0, 
           y0 = 20, 
           x1 = hap1_length, 
           y1 = 20, 
           lwd = 3,
           col = "red") # transparent("red", trans.val = 0.3))
  
  segments(x0 = 0, 
           y0 = 2, 
           x1 = hap2_length, 
           y1 = 2, 
           lwd = 3,
           col = "steelblue") #transparent("steelblue", trans.val = 0.3))
  
  axis(side = 1,
       at = seq(0, max(hap1_length, hap2_length) + 25*10^6, 25*10^6),
       labels = paste(floor(seq(0, max(hap1_length, hap2_length) + 25*10^6, 25*10^6)/10^6), "", sep = ""),
       cex.axis = 1.5) 
  
  # add gap to both haplotypes
  gap_hap1_chr <- mygap_hap1[mygap_hap1$V1 == chr, ]
  
  for(j in 1:nrow(gap_hap1_chr)){
    
    segments(x0 = (gap_hap1_chr[j, 2] + gap_hap1_chr[j, 3])/2, 
             y0 = 20.2, 
             x1 = (gap_hap1_chr[j, 2] + gap_hap1_chr[j, 3])/2, 
             y1 = 20.7, 
             lwd = 1,
             col = "black")
    
  }
  
  gap_hap2_chr <- mygap_hap2[mygap_hap2$V1 == chr2, ]
  
  for(j in 1:nrow(gap_hap2_chr)){
    
    segments(x0 = (gap_hap2_chr[j, 2] + gap_hap2_chr[j, 3])/2, 
             y0 = 1.8, 
             x1 = (gap_hap2_chr[j, 2] + gap_hap2_chr[j, 3])/2, 
             y1 = 1.3, 
             lwd = 1,
             col = "black")
    
  }
  
  if(mylable){
    mtext("Rabiosa h1", side = 2, line = 0, las = 2, at = 20, cex = 1)
    mtext("Rabiosa h2", side = 2, line = 0, las = 2, at = 2, cex = 1)
  }
  
  # make bin of 100 Kb size
  mybin <- seq(1, hap1_length, 10^5)
  mybin_M <- cbind(mybin[1:(length(mybin)-1)], mybin[2:length(mybin)])
  if(mybin_M[dim(mybin_M)[1], 2] < hap1_length){
     
    # if the last value is lower than the total length, then append a new bin to the matrix
    mybin_M <- rbind(mybin_M, c(mybin_M[dim(mybin_M)[1], 2], hap1_length))
    
  }else{
    
    # otherwise, don't do anything
    mybin_M <- mybin_M
    
  }
  
  # mycount_fun

  
  # calculate how many variants in every bin
  mysnp_chr <- mysnp[mysnp$V1 == chr, ]
  myindel_chr <- myindel[myindel$V1 == chr, ]
  mypav_chr <- mypav[mypav$V1 == chr, ]
  
  mycount_SNP <- unlist(apply(mybin_M, 1, mycount_fun, myvar = mysnp_chr$V2))
  mycount_INDEL <- unlist(apply(mybin_M, 1, mycount_fun, myvar = myindel_chr$V2))
  mycount_pav <- unlist(apply(mybin_M, 1, mycount_fun, myvar = mypav_chr$V2))
  
  # draw SNP distribution
  myfit <- round(max(mycount_SNP)/5)
  polygon(x = c(0, mybin_M[, 1], max(mybin_M[, 1])), y = c(24, 24 + mycount_SNP/myfit, 24), border = NA, col = "darkmagenta")
  
  if(mylable){
    axis(side = 2, at = seq(24, 24 + 5, 5),
         labels = NA, line = -1,
         las = 2, cex.axis = 1, tck = -0.01)
    
    for(j in 1:length(seq(24, 24 + 5, 5))){
      
      mtext(text = paste(round(myfit*(seq(0, 5, 5)[j])/1000, 1), "k", sep = ""), 
            at = seq(24, 24 + 5, 5)[j],
            side = 2, cex = 1, line = -0.7, las = 2)
    }
  }  
  
  # add a y-lable for SNP
  if(mylable){
    mtext("SNP", side = 2, line = 1, las = 2, at = (24 + 24 + 5)/2, cex = 1)
  }
  
  # add small indels
  myfit <- round(max(mycount_INDEL)/5)
  polygon(x = c(0, mybin_M[, 1], max(mybin_M[, 1])), y = c(32, 32 + mycount_INDEL/myfit, 32), border = NA, col = "darkmagenta")
  
  if(mylable){
    axis(side = 2, at = seq(32, 32 + 5, 5),
         labels = NA, line = -1,
         las = 2, cex.axis = 1, tck = -0.01)
    
    for(j in 1:length(seq(32, 32 + 5, 5))){
      
      mtext(text = paste(round(myfit*(seq(0, 5, 5)[j])/1000, 2), "k", sep = ""), 
            at = seq(32, 32 + 5, 5)[j],
            side = 2, cex = 1, line = -0.7, las = 2)
      
    }
  }
  
  # add a y-lable for indel
  if(mylable){
    mtext("INDEL", side = 2, line = 1, las = 2, at = (32 + 32 + 5)/2, cex = 1)
  }
  
  # add PAV
  myfit <- round(max(mycount_pav)/5)
  polygon(x = c(0, mybin_M[, 1], max(mybin_M[, 1])), y = c(40, 40 + mycount_pav/myfit, 40), border = NA, col = "darkmagenta")
  
  if(mylable){
    axis(side = 2, at = seq(40, 40 + 5, 5),
         labels = NA, line = -1,
         las = 2, cex.axis = 1, tck = -0.01)
    
    for(j in 1:length(seq(40, 40 + 5, 5))){
      
      mtext(text = myfit*(seq(0, 5, 5)[j]), 
            at = seq(40, 40 + 5, 5)[j],
            side = 2, cex = 1, line = -0.7, las = 2)
    }
  }
  
  # add a y-lable for PAV
  if(mylable){
    mtext("PAV", side = 2, line = 1, las = 2, at = (40 + 40 + 5)/2, cex = 1)
  }
  
  
  # add alignment identity
  # calculate sequence identity between two haplotypes in very 100 Kb region
  # identity = aligned sequence/100Kb
  # aligned sequence = 100kb - total length of SNP, INDEL, PAV, gaps and notal
  mynotal_chr <- mynotal[mynotal$V1 == chr, ]
  
  # calculated how many variants in every window
  size_SNP <- as.numeric(unlist(apply(mybin_M, 1, Intersect_variant, myvar = mysnp_chr)))
  size_INDEL <- as.numeric(unlist(apply(mybin_M, 1, Intersect_variant, myvar = myindel_chr)))
  size_PAV <- as.numeric(unlist(apply(mybin_M, 1, Intersect_variant, myvar = mypav_chr)))
  size_gap <- as.numeric(unlist(apply(mybin_M, 1, Intersect_variant, myvar = mygap_hap1)))
  size_NOTAL <- as.numeric(unlist(apply(mybin_M, 1, Intersect_variant, myvar = mynotal_chr)))
  
  myidentity <- as.numeric(100000 - c(size_SNP + size_INDEL + size_PAV + size_gap + size_NOTAL))
  myidentity[myidentity <= 0] <- 0
  myidentity <- as.numeric(myidentity/100000)
  
  points(x = c(0, sort(c(as.numeric(mybin_M[, 1]), as.numeric(mybin_M[, 2]))), max(as.numeric(mybin_M[, 2]))),
         y = c(48, rep(48 + myidentity*5, each = 2), 48), type = "s", col = "gray40")

  
  # make an axis for the identity
  if(mylable){
    axis(side = 2,
         at = seq(48, 48 + 5, 2.5),
         labels = NA,
         cex.axis = 1, line = -1,
         las = 2, tck = -0.01)
    
    for(j in 1:length(seq(48, 48 + 5, 2.5))){
      
      mtext(text = seq(0, 100, 50)[j], 
            at = seq(48, 48 + 5, 2.5)[j],
            side = 2, cex = 1, line = -0.7, las = 2)
    }
  }
  
  # add a reference line
  segments(x0 = 1, y0 = 48 + 5*0.3, x1 = hap1_length, y1 = 48 + 5*0.3, lty = 2, col = "red")
  
  # add a y-lable for alignment identity
  if(mylable){
    mtext("Sequence\nidentity", side = 2, line = 1, las = 2, at = (48 + 48 + 5)/2, cex = 1)
  }
  
  
  # add legend
  legend(x = max(hap1_length, hap2_length)*0.97, y = 15,
         bty = "n",
         legend = c("Synteny", "Inversion", "Gap", 'Not aligned'),
         pt.bg = c("lightgray", "gold", "black", 'white'), border = 'black',
         cex = 1,
         pch = 22,
         xpd = t, 
         xjust = 0,
         yjust = 1,
         x.intersp = c(0.3, 0.3, 0.3, 0.3),
         y.intersp = c(1, 1, 1, 1))
  
}

mychr7 <- Draw_synteny(mysyn,
             myinv,
             chr = paste("rabiosa_1_chr", 7, sep = ""),
             chr2 = paste("rabiosa_2_chr", 7, sep = ""),
             mygap_hap1,
             mygap_hap2,
             hap1_length = myfai[myfai$V1 == paste("rabiosa_1_chr", 7, sep = ""), 2],
             hap2_length = myfai[myfai$V1 == paste("rabiosa_2_chr", 7, sep = ""), 2],
             mynotal,
             mysnp,
             myindel,
             mypav,
             mylable = T)

# pdf("Rabiosa_hap_synteny_new.pdf", 15, 5)
# for(i in 1:7){
#   
#   Draw_synteny(mysyn,
#                myinv,
#                chr = paste("rabiosa_1_chr", i, sep = ""),
#                chr2 = paste("rabiosa_2_chr", i, sep = ""),
#                mygap_hap1,
#                mygap_hap2,
#                hap1_length = myfai[myfai$V1 == paste("rabiosa_1_chr", i, sep = ""), 2],
#                hap2_length = myfai[myfai$V1 == paste("rabiosa_2_chr", i, sep = ""), 2],
#                mynotal,
#                mysnp,
#                myindel,
#                mypav,
#                mylable = T)
#   
# }
# dev.off()


for(i in 1:7){
  # png(paste("Rabiosa_hap_synteny_chr", i, ".png", sep = ''), width = 15, height = 5, units = 'in', res = 600)
  pdf(paste("Rabiosa_hap_synteny_chr", i, ".pdf", sep = ''), width = 15, height = 5)
  Draw_synteny(mysyn,
               myinv,
               chr = paste("rabiosa_1_chr", i, sep = ""),
               chr2 = paste("rabiosa_2_chr", i, sep = ""),
               mygap_hap1,
               mygap_hap2,
               hap1_length = myfai[myfai$V1 == paste("rabiosa_1_chr", i, sep = ""), 2],
               hap2_length = myfai[myfai$V1 == paste("rabiosa_2_chr", i, sep = ""), 2],
               mynotal,
               mysnp,
               myindel,
               mypav,
               mylable = T)
  dev.off()
}



# pdf("Rabiosa_hap_synteny_new.pdf", 10, 20)
# layout(matrix(c(1:7), 7, 1, byrow = TRUE))
# for(i in 1:7){
#   
#   Draw_synteny(mysyn,
#                myinv,
#                chr = paste("rabiosa_1_chr", i, sep = ""),
#                chr2 = paste("rabiosa_2_chr", i, sep = ""),
#                mygap_hap1,
#                mygap_hap2,
#                hap1_length = myfai[myfai$V1 == paste("rabiosa_1_chr", i, sep = ""), 2],
#                hap2_length = myfai[myfai$V1 == paste("rabiosa_2_chr", i, sep = ""), 2],
#                mynotal,
#                mysnp,
#                myindel,
#                mypav,
#                mylable = T)
#   
# }
# dev.off()
