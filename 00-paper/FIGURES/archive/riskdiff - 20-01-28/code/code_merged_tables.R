for (table in c("Table1","Table2","Table3")){
  file_list <- list.files(paste0("../output/",table))
  dataset <- do.call("rbind", lapply(file_list,FUN=function(files){read.csv(paste0("../output/",table,"/",files), header=TRUE)}))
  write.csv(dataset,file=paste0("../output/",table,"_merged.csv"))
}

# for (table in c("Table1","Table2","Table3")){
#   file_list <- list.files(paste0("../output/",table))
#   dataset <- do.call("cbind",lapply(file_list,FUN=function(files){read.csv(paste0("../output/",table,"/",files), header=TRUE)[,2]}))
#   row.names(dataset)<-read.csv(paste0("../output/",table,"/1",table,".csv"))[,1]
#   colnames(dataset)<-gsub(paste0(table,".csv"), "", file_list)
#   write.csv(dataset,file=paste0("../output/",table,"_merged.csv"))
# }