myargs <- commandArgs(trailingOnly=TRUE)

# myargs[1] = Rabiosa_syri.pav.filtered
# myargs[2] = sniffles_called_sv.txt

SYRIPAV <- read.table(myargs[1], header = F, stringsAsFactors = F)
SNIFFLESPAV <- read.table(myargs[2], header = F, stringsAsFactors = F)

# SYRIPAV <- read.table("Rabiosa_syri.pav.filtered", header = F, stringsAsFactors = F)
# SNIFFLESPAV <- read.table("sniffles_called_sv.txt", header = F, stringsAsFactors = F)

SYRIPAV_del <- SYRIPAV[SYRIPAV$V5 == "DEL", ]
SYRIPAV_ins <- SYRIPAV[SYRIPAV$V5 == "INS", ]

SNIFFLESPAV_del <- SNIFFLESPAV[SNIFFLESPAV$V5 == "DEL", ]
SNIFFLESPAV_ins <- SNIFFLESPAV[SNIFFLESPAV$V5 == "INS", ]


myfun_ins <- function(x, SNIFFLESPAV_chr = SNIFFLESPAV_chr){
  
  # make the start and end site of INS 10 bp more open for testing 
  start_test <- unique(as.numeric(x[2]) >= (SNIFFLESPAV_chr[, 2] - 10) & as.numeric(x[2]) <= (SNIFFLESPAV_chr[, 3] + 10))
  end_test <- unique(as.numeric(x[3]) >= (SNIFFLESPAV_chr[, 2] - 10) & as.numeric(x[3]) <= (SNIFFLESPAV_chr[, 3] + 10))
  
  if(length(start_test) > 1 | length(end_test) > 1) {
    
    return("common")
    
  }else{
    
    return("unique")
    
  }
  
}

myfun_del <- function(x, SNIFFLESPAV_chr = SNIFFLESPAV_chr){
  
  # one deletion will be considered as a common one only when there is an overlap between sniffles and syri
  overlap1 <- unique(as.numeric(x[2]) < SNIFFLESPAV_chr[, 2] & (as.numeric(x[3]) <= SNIFFLESPAV_chr[, 3] & as.numeric(x[3]) >= SNIFFLESPAV_chr[, 2]))
  overlap2 <- unique(as.numeric(x[2]) < SNIFFLESPAV_chr[, 2] & as.numeric(x[3]) > SNIFFLESPAV_chr[, 3])
  overlap3 <- unique((as.numeric(x[2]) >= SNIFFLESPAV_chr[, 2] & as.numeric(x[2]) <= SNIFFLESPAV_chr[, 3]) & as.numeric(x[3]) > SNIFFLESPAV_chr[, 3])
  overlap4 <- unique((as.numeric(x[2]) >= SNIFFLESPAV_chr[, 2] & as.numeric(x[2]) <= SNIFFLESPAV_chr[, 3]) & (as.numeric(x[3]) <= SNIFFLESPAV_chr[, 3] & as.numeric(x[3]) >= SNIFFLESPAV_chr[, 2]))
  
  if(length(overlap1) > 1 | length(overlap2) > 1 | length(overlap3) > 1 | length(overlap4) > 1) {
    
    return("common")
    
  }else{
    
    return("unique")
    
  }
}


mychr <- sort(unique(SYRIPAV$V1))

myins <- NULL
mydel <- NULL

for(x in mychr){
  
  SYRIPAV_del_chr <- SYRIPAV_del[SYRIPAV_del$V1 == x, ]
  SNIFFLESPAV_del_chr <- SNIFFLESPAV_del[SNIFFLESPAV_del$V1 == x, ]
  
  SYRIPAV_ins_chr <- SYRIPAV_ins[SYRIPAV_ins$V1 == x, ]
  SNIFFLESPAV_ins_chr <- SNIFFLESPAV_ins[SNIFFLESPAV_ins$V1 == x, ]
  
  SYRIPAV_del_chr$V6 <- unlist(apply(SYRIPAV_del_chr, 1, myfun_del, SNIFFLESPAV_chr = SNIFFLESPAV_del_chr))
  
  SYRIPAV_ins_chr$V6 <- unlist(apply(SYRIPAV_ins_chr, 1, myfun_ins, SNIFFLESPAV_chr = SNIFFLESPAV_ins_chr))
  
  mydel <- rbind.data.frame(mydel, SYRIPAV_del_chr)
  myins <- rbind.data.frame(myins, SYRIPAV_ins_chr)
  
}

mypav <- rbind.data.frame(myins, mydel)
mypav <- mypav[order(mypav$V1, mypav$V2), ]

write.table(mypav, "PAV_overlap.txt", quote = F, sep = "\t", row.names = F, col.names = F)
