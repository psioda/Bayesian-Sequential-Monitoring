setwd("/nas/longleaf/home/ekwiatko/FDA/riskdiff/output/Table3/")

file_list <- list.files("/nas/longleaf/home/ekwiatko/FDA/riskdiff/output/Table3/")

dataset <- do.call("cbind",lapply(file_list,
           FUN=function(files){read.csv(files, header=TRUE)[,2]}))

row.names(dataset)<-read.csv("1Table3.csv")[,1]
colnames(dataset)<-gsub("Table3.csv", "", file_list)

write.csv(dataset,file="/nas/longleaf/home/ekwiatko/FDA/riskdiff/latex/3dataset.csv")
