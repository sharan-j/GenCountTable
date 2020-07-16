suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(Biostrings)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(dplyr)))

args = commandArgs(trailingOnly=TRUE)

reference = read.csv(args[1], sep = "\t", stringsAsFactors = F, header = F, row.names = 1)
colnames(reference) = c( "seq")
#head(reference)

param = ScanBamParam(what = c("rname", "seq"))
myBam = scanBam(file = args[2], param = param)
GuideCount = as.matrix((table(myBam[[1]]$rname)))

print("Generating the count table...", quote=FALSE)
count.mat = NULL
for (i in row.names(GuideCount)) {
  count = as.matrix(GuideCount[i,])
  sequences = DNAStringSet(myBam[[1]]$seq[which(myBam[[1]]$rname == i)])
  con.mat = consensusMatrix(sequences, as.prob = T)
  con.mat = round(con.mat, 3)
  withCount = cbind(count,t(con.mat["G",]))
  count.mat = rbind(count.mat,withCount)
}

colnames(count.mat) = c("allCounts", c(1:20))
CountTable = merge(reference,count.mat, by = "row.names", all.y = T)
CountTable = rename(CountTable, ID = Row.names)

CT_splt = colsplit(CountTable$seq, "", names = c(1:20))
CT_bin = as.matrix(CountTable[4:23])
CT_bin[CT_splt != "A"] = 0

CountTable[4:23] = CT_bin
write.table(CountTable, args[3] , sep = "\t", row.names = F)


