##################################
# Model parameters
# Evan Kwiatkowski, Jun 22, 2020
##################################

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

file_list <- list.files(paste0("../output/Table0"))

for (i in 1:length(file_list)){
  temp <- read.tcsv(paste0("../output/Table0/", file_list[i]), header = TRUE)
  temp$idx <- temp[nrow(temp),2]
  temp <- temp[1:(nrow(temp)-1),]
  
  if (i==1) {
    final <- temp
  } else {
    final <- rbind(final,temp)
  }
}

args_simulation <- read.csv(file = "args_simulation.csv", header = TRUE, sep = ",")
final           <- final[,-1]
combined        <- merge(args_simulation, final, by.x = "X", by.y = "idx")

write.csv(combined, file = "Table0_merged.csv") # will be used to bring back from longleaf
