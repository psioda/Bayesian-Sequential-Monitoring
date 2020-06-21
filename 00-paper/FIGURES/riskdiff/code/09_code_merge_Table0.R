read.tcsv = function(file, header=TRUE, sep=",", ...) {
  
  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)
  
  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }
  
  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep) 
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
  
}

table <- "Table0_052020"
file_list <- list.files(paste0("../output/", table))

for (i in 1:200){
  temp <- read.tcsv(paste0("../output/",table,"/", file_list[i]), header = TRUE)
  temp$idx <- temp[26,2]
  temp <- temp[1:25,]
  
  if (i==1) {
    final <- temp
  } else {
    final <- rbind(final,temp)
  }
}


args_simulation <- read.csv(file = "../output/args_simulation_052020.csv", header = TRUE, sep = ",")
final <- final[,-1]
combined <- merge(args_simulation, final, by.x = "X", by.y = "idx")
View(combined)
row1 <- colMeans(combined[combined$p.IP==0.39,])
row2 <- colMeans(combined[combined$p.IP==0.45,])
row3 <- colMeans(combined[combined$p.IP==0.51,])
row4 <- colMeans(combined[combined$p.IP==0.57,])
row5 <- colMeans(combined[combined$p.IP==0.63,])