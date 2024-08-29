#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# args[1] = reference paf file
# args[2] = swap paf file
# args[3] = query prefix
# args[4] = reference prefix
# args[5] = output paf file

mypaf <- read.table(args[1], header = F, stringsAsFactors = F)
myswappaf <- read.table(args[2], header = F, stringsAsFactors = F)

if(nrow(mypaf) == nrow(myswappaf)){
  
  print("Both pafs have same row number, now add query length, start and end to reference paf")
  
  # from myswappaf, add query length, start and end to mypaf
  mypaf$V2 <- myswappaf$V7
  mypaf$V3 <- myswappaf$V8
  mypaf$V4 <- myswappaf$V9
  
  # change query name
  mypaf$V1 <- sub(unique(mypaf$V1), paste(args[3], unique(mypaf$V1), sep = ""), mypaf$V1)
  
  # change reference name
  mypaf$V6 <- sub(unique(mypaf$V6), paste(args[4], unique(mypaf$V6), sep = ""), mypaf$V6)
  
  write.table(mypaf, args[5], quote = F, row.names = F, col.names = F, sep = "\t")
  
}else{
  
  print(paste(args[1], args[2], "don't have same length, check your paf"))
  
}