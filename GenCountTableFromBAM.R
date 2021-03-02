suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(foreach)))

args = commandArgs(trailingOnly=TRUE)

referenceAll = read.csv(args[1], sep = "\t", stringsAsFactors = F, header = F)
colnames(referenceAll) = c("rname","seq")

param = ScanBamParam(what = c("rname", "seq", "qual", "cigar"), tag=c("MD"))
myBam = scanBam(file = args[2], param = param)

outfile = file_path_sans_ext(args[2])


myBam = as.data.frame(myBam)
colnames(myBam) = c("rname", "cigar", "seq", "qual", "MD")
head(myBam)
myBam$seq = substr(myBam$seq,7,26) # Trim to 20 bp protospacer if required
#myBam = myBam[!grepl("\\d", myBam$qual),]
#myBam = myBam[!grepl("[:punct:]", myBam$qual),]
myBam = myBam[!grepl("(G|C|T)", myBam$MD),]
myBam = myBam[!grepl("(I|D)", myBam$cigar),]

readLength = 20

dt = data.table::as.data.table(myBam[,c("rname", "seq")])
dt = dt[, .N, by = c('rname','seq')]
dt = unique(dt)
dt = dt[order(dt$rname, -dt$N),]
head(dt)

dt.mat = cbind(dt, matrix(unlist(strsplit(as.character(dt$seq), "")), ncol = readLength, byrow = T),stringsAsFactors = F)
head(dt.mat)

rnameInBam = as.data.frame(unique(dt.mat$rname))
colnames(rnameInBam) = c("rname")
reference = merge(rnameInBam, referenceAll, by = "rname", all.x = T)
ref.mat = cbind(reference, matrix(unlist(strsplit(as.character(reference$seq), "")), ncol = readLength, byrow = T),stringsAsFactors = F)
head(ref.mat)

system.time({
  foreach (i = as.character(ref.mat$rname)) %do% {
    a = as.matrix(dt.mat[which(dt.mat$rname == i), 3:(readLength+3)])
    b = as.matrix(ref.mat[which(ref.mat$rname == i), 3:(readLength+2)])
    print(i)
    for (j in 1:nrow(a))  {
      
      for (k in 1:(readLength))  {
        if (a[j,k+1] == b[,k]) {
          a[j,k+1] = 0
        } else if (a[j,k+1] == "G"){
          a[j,k+1] = a[j,1]  
        } else {a[j,k+1] = NA}
      }
    }
    mode(a) = "numeric"
    dt.mat[which(dt.mat$rname == i), 3:(readLength+3) ] = as.data.frame(a)
  }
})

head(dt.mat)
dt.filt = as.data.frame(na.omit(dt.mat))
dt.filt[, 4:(readLength+3)] = sapply(dt.filt[, 4:(readLength+3)], as.numeric)

system.time({
foreach (i = as.character(ref.mat$rname)) %do% {
  print(i)
  dt.filt[which(dt.filt$rname == i),"allCounts"] = sum(dt.filt[which(dt.filt$rname == i),"N"])
}
})

head(dt.filt)
dt.filt$prop = dt.filt$N/dt.filt$allCounts
propTable = dt.filt[,c("rname", "seq", "allCounts", "prop")]
head(propTable)
propTable = merge(propTable, reference, by = "rname", all.x = T)
colnames(propTable) = c("rname", "seq", "allCouts", "prop", "refSeq")

count.mat = aggregate(dt.filt[,3:(readLength+3)], list(dt.filt$rname), FUN = sum)
head(count.mat)
colnames(count.mat) = c("rname","allCounts", 1:(readLength) )
count.mat = merge(referenceAll, count.mat, by = "rname", all.x = T)

#write.table(propTable, cat(outfile, "_ProportionTable.csv", sep = ""), sep = "\t", row.names = F)
write.table(count.mat, cat(outfile, "_CountTable.csv", sep = ""), sep = "\t", row.names = F)



