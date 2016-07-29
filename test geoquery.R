setwd("/Users/gomk/Documents/table")
library(GEOquery)
# import raw data
rawdata= getGEOSuppFiles(GEO='GSE18388', makeDirectory = TRUE, baseDir = getwd())
# import meta data
untar("./GSE18388/GSE18388_RAW.tar")
# gds <- getGEO("GSE18388") # all normalized data + metadata
# gsm <- getGEO(filename=system.file("extdata/GSE18388.txt.gz",package="GEOquery"))

gds <- getGEO("GSE18388",GSEMatrix = F,getGPL=T,AnnotGPL=T)

GSMList(gds)[[1]] # sample 1
GSMList(gds)[[2]] # sample 2
#head(Meta(gds))
#Columns(gds)
Meta(gds)$platform_id

Table(GSMList(gds)[[1]])[1:5,]
Meta(GSMList(gds)[[2]])$geo_accession # for sample 1
Meta(GSMList(gds)[[1]])$title
Meta(GSMList(gds)[[1]])$description
length(GSMList(gds))
mytable=matrix("",length(GSMList(gds)),4)
colnames(mytable)=c("gsm","title","description","group")
for (k in 1:length(GSMList(gds)))
{
mytable[k,] <-c(Meta(GSMList(gds)[[k]])$geo_accession, Meta(GSMList(gds)[[k]])$title, Meta(GSMList(gds)[[k]])$description, "")
  

}

mytable <- data.frame(mytable)

