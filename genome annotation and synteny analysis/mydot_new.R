# make dot plot based on collinear blocks

library(stringr)

# set working directory
setwd("P:/Yutangchen/Rabiosa_compare/All/compare_all")

# import bed and anchor file
mybed <- read.table("assembly1.bed", header = F, stringsAsFactors = F)
myanchor <- read.table("assembly1.assembly2.anchors", header = F, stringsAsFactors = F)

# make a query and subject column
myspecies <- c("HORVU", "KYUS", "Lm0", "Lm1", "Lm2", "V3.Lp")

# make new name of species for later plot naming, order of species in the vector should be same as that in myspecies 
# myspecies_new <- c("H.vulgare Morex", 
#                    "L.perenne Kyuss",
#                    "L.multiflorum Rabiosa unphased",
#                    "L.multiflorum Rabiosa hap1",
#                    "L.multiflorum Rabiosa hap2",
#                    "L.perenne 261")

myspecies_new <- c("Morex v3", 
                   "Kyuss",
                   "Rabiosa v2",
                   "Rabiosa h1",
                   "Rabiosa h2",
                   "Lp261")

# the aim of this function is to add species column to myanchor to make it easier to subset myanchor
Add_species_myanchor <- function(myanchor = myanchor, myspecies = myspecies){
  
  xx <- myanchor$V1
  for(x in myspecies){
    
    xx <- sub(paste(x, ".*", sep = ""), x, xx)
    
  }
  
  myanchor$V3 <- xx
  
  # subject
  xx <- myanchor$V2
  for(x in myspecies){
    
    xx <- sub(paste(x, ".*", sep = ""), x, xx)
    
  }
  
  myanchor$V4 <- xx
  
  return(myanchor)
  
}

# use the mid value to represent the position of a gene in mybed
Make_bed <- function(mybed = mybed){
  
  # use the mid value to represent the position of a gene in mybed
  mybed$V5 <- floor((mybed$V3 + mybed$V2)/2)
  colnames(mybed) <- c("chr", "start", "end", "gene", "mid", "order")
  return(mybed)
  
}

# make a function to subeset myanchor by species and meanwhile can replace the species name with new names for plot naming
Subset_anchor <- function(myanchor = myanchor, 
                          myspecies = myspecies,
                          myspecies_new = myspecies_new, 
                          qsp = myspecies[3:5], 
                          ssp = myspecies[c(1, 2, 6)], 
                          mybed = mybed){
  
  myanchor_r <- myanchor[myanchor$V3 %in% qsp & myanchor$V4 %in% ssp, ] 
  myanchor_r <- myanchor_r[order(myanchor_r$V3, myanchor_r$V4), ]
  
  # based on the gene model name, find corresponding postion and chr in mybed and add that to myanchor_r
  # use merge, but first need to give column names
  
  colnames(myanchor_r) <- c("qg", "sg", "qsp", "ssp")
  
  
  myanchor_rq <- merge(myanchor_r, mybed, by.x = "qg", by.y = "gene", sort = F)
  myanchor_rq <- myanchor_rq[c(1, 2, 3, 4, 5, 8)]
  colnames(myanchor_rq) <- c("qg", "sg", "qsp", "ssp", "qchr", "qmid")
  
  myanchor_rqs <- merge(myanchor_rq, mybed, by.x = "sg", by.y = "gene", sort = F)
  myanchor_rqs <- myanchor_rqs[c(1:7, 10)]
  colnames(myanchor_rqs) <- c("sg", "qg", "qsp", "ssp", "qchr", "qmid", "schr", "smid")
  
  # change the name of species in myanchor_rqs for ploting naming
  for(i in 1:length(myspecies)){
    
    myanchor_rqs$qsp[myanchor_rqs$qsp == myspecies[i]] <- myspecies_new[i]
    myanchor_rqs$ssp[myanchor_rqs$ssp == myspecies[i]] <- myspecies_new[i]
    
  }
  
  return(myanchor_rqs)
  
}

# function to make dot plot
# xlt means x axis label and tick
make_dot <- function(df = df, 
                     qsp = "species", 
                     ssp = "species", 
                     segcol = "blue",
                     dotcol = "black",
                     xlt = T, 
                     ylt = T,
                     mymar = c(4, 5, 4, 3)){
  
  df <- df[df$qsp == qsp  & df$ssp == ssp, ]
  
  df <- df[order(df$qchr), ]
  
  qmid_max <- aggregate(qmid ~ qchr, df, max)
  smid_max <- aggregate(smid ~ schr, df, max)
  
  qmid_max$position <- cumsum(c(0, qmid_max$qmid[1:(nrow(qmid_max)-1)]))
  smid_max$position <- cumsum(c(0, smid_max$smid[1:(nrow(smid_max)-1)]))
  
  for(i in 1:nrow(qmid_max)){
    
    df$qchr[df$qchr == qmid_max[i, 1]] <- qmid_max[i, 3]
    df$schr[df$schr == smid_max[i, 1]] <- smid_max[i, 3]
    
  }
  
  df$qchr <- as.numeric(df$qchr)
  df$schr <- as.numeric(df$schr)
  
  df$qpos <- df$qchr + df$qmid
  df$spos <- df$schr + df$smid
  
  par(xpd = T, mar = mymar)
  plot(df$spos, df$qpos, 
       pch = 16, cex = 0.2,
       axes = F, ann = F,
       xlim = c(0, max(df$spos)),
       ylim = c(0, max(df$qpos)),
       col = adjustcolor(dotcol, alpha.f = 0.3))
  # type = "c")
  
  # col = mycol)
  
  for(i in 1:(nrow(qmid_max)+1)){
    
    # add vertical lines
    segments(x0 = c(0, cumsum(smid_max$smid))[i], 
             y0 = 0, 
             x1 = c(0, cumsum(smid_max$smid))[i], 
             y1 = max(cumsum(qmid_max$qmid)),
             lty = 5,
             col = segcol) # "black") # mycol
    
    #add horizontal lines
    segments(x0 = 0, 
             y0 = c(0, cumsum(qmid_max$qmid))[i], 
             x1 = max(cumsum(smid_max$smid)), 
             y1 = c(0, cumsum(qmid_max$qmid))[i],
             lty = 5,
             col = segcol) # "black") # mycol
    
  }
  
  # add text to x, y axis
  
  if(xlt){
    
    text(x = floor(cumsum(smid_max$smid) - smid_max$smid/2), y = 0, 
         labels = str_to_title(gsub("*.*_", "", smid_max$schr)),
         cex = 1, pos = 1)
    text(x = floor(max(cumsum(smid_max$smid))/2), y = 0,
         labels = ssp, cex = 1.5, pos = 1, offset = 2)
    
  }
  
  if(ylt){
    
    text(x = 0, y = floor(cumsum(qmid_max$qmid) - qmid_max$qmid/2), 
         labels = str_to_title(gsub("*.*_", "", qmid_max$qchr)),
         cex = 1, pos = 2)
    text(x = 0, y = floor(max(cumsum(qmid_max$qmid))*2/3),
         labels = qsp, cex = 1.5, pos = 2, offset = 4, srt = 90)
    
  }
  
}

#--------------------------------------------------------------------------------------------------#

myanchor <- Add_species_myanchor(myanchor = myanchor, myspecies = myspecies)
mybed <- Make_bed(mybed = mybed)

#######################################################
dotanchor <- Subset_anchor(myanchor = myanchor,
                           myspecies = myspecies,
                           myspecies_new = myspecies_new,
                           qsp = myspecies, 
                           ssp = myspecies, 
                           mybed = mybed)

pdf("Rabiosa_vs_all.pdf", width = 15, height = 12)
# png("Rabiosa_vs_all.png", width = 15, height = 12, units = 'in', res = 480)
layout(matrix(1:9, nrow = 3, ncol = 3, byrow = T))

make_dot(df = dotanchor, ssp = "Morex v3", qsp = "Rabiosa v2", segcol = "darkseagreen4", mymar = c(2, 5, 4, 2))
make_dot(df = dotanchor, ssp = "Morex v3", qsp = "Rabiosa h1", segcol = "darkseagreen4", mymar = c(2, 5, 4, 2))
make_dot(df = dotanchor, ssp = "Morex v3", qsp = "Rabiosa h2", segcol = "darkseagreen4", mymar = c(2, 5, 4, 2))

make_dot(df = dotanchor, ssp = "Lp261", qsp = "Rabiosa v2", segcol = "darkseagreen4", mymar = c(2, 5, 2, 2))
make_dot(df = dotanchor, ssp = "Lp261", qsp = "Rabiosa h1", segcol = "darkseagreen4", mymar = c(2, 5, 2, 2))
make_dot(df = dotanchor, ssp = "Lp261", qsp = "Rabiosa h2", segcol = "darkseagreen4", mymar = c(2, 5, 2, 2))


make_dot(df = dotanchor, ssp = "Kyuss", qsp = "Rabiosa v2", segcol = "darkseagreen4", mymar = c(4, 5, 2, 2))
make_dot(df = dotanchor, ssp = "Kyuss", qsp = "Rabiosa h1", segcol = "darkseagreen4", mymar = c(4, 5, 2, 2))
make_dot(df = dotanchor, ssp = "Kyuss", qsp = "Rabiosa h2", segcol = "darkseagreen4", mymar = c(4, 5, 2, 2))


dev.off()

# synteny between Rabiosa v2 and Rabiosa h1 and h2
png("Rabiosa_vs_h1h2.png", width = 10, height = 5, units = 'in', res = 480)
layout(matrix(1:2, nrow = 1, ncol = 2, byrow = T))
make_dot(df = dotanchor, ssp = "Rabiosa h1", qsp = "Rabiosa v2", segcol = "darkred", mymar = c(4, 5, 2, 3))
make_dot(df = dotanchor, ssp = "Rabiosa h2", qsp = "Rabiosa v2", segcol = "darkred", mymar = c(4, 5, 2, 3))
dev.off()

# synteny between h1 and h2
png("h1_vs_h2.png", width = 9, height = 9, units = 'in', res = 600)
make_dot(df = dotanchor, ssp = "Rabiosa h1", qsp = "Rabiosa h2", segcol = "darkseagreen4", mymar = c(4, 5, 2, 3))
dev.off()

