
#install the package and call the library
install.packages("Rsubread")
library("Rsubread")

#set input and output directory
dir_input <- "/dfs6/pub/byukimti/alignGenome/finishedAligning/sortedBams/"
dir_output <- "/dfs6/pub/byukimti/alignGenome/finishedAligning/rsubread/"

#set working directory 
setwd(dir_output)

#set the gtf file
gtf.file <- "/dfs6/pub/byukimti/Mus_musculus.GRCm39.108.gtf"

#assign variables to each .bam file
ZO04 <- paste0(dir_input, "304Aligned.sortedByCoord.out.bam")
ZO05 <- paste0(dir_input, "305Aligned.sortedByCoord.out.bam")
ZO06 <- paste0(dir_input, "306Aligned.sortedByCoord.out.bam")
IO07 <- paste0(dir_input, "307Aligned.sortedByCoord.out.bam")
IO08 <- paste0(dir_input, "308Aligned.sortedByCoord.out.bam")
IO09 <- paste0(dir_input, "309Aligned.sortedByCoord.out.bam")
ZI10 <- paste0(dir_input, "310Aligned.sortedByCoord.out.bam")
ZI11 <- paste0(dir_input, "311Aligned.sortedByCoord.out.bam")
ZI12 <- paste0(dir_input, "312Aligned.sortedByCoord.out.bam")
II13 <- paste0(dir_input, "313Aligned.sortedByCoord.out.bam")
II14 <- paste0(dir_input, "314Aligned.sortedByCoord.out.bam")
II15 <- paste0(dir_input, "315Aligned.sortedByCoord.out.bam")
bam.files <- c(ZO04,ZO05,ZO06,IO07,IO08,IO09,ZI10,ZI11,ZI12,II13,II14,II15)
