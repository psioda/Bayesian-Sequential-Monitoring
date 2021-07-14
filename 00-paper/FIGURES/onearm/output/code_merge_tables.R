for (table in c("Table1","Table2","Table3")){
  file_list <- list.files(paste0("../output/",table))
  dataset <- do.call("rbind", lapply(file_list,FUN=function(files){read.csv(paste0("../output/",table,"/",files), header=TRUE)}))
  write.csv(dataset,file=paste0(table,"_merged.csv"))
}
