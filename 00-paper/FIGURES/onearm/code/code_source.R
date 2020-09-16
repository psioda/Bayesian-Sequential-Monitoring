rm(list = ls())

# library(foreach)
# library(doParallel)
# registerDoParallel(detectCores())
# getDoParWorkers()
# 
# start_time      <- Sys.time()
args_simulation <- read.csv(file="args_simulation.csv",header=TRUE,sep=",")

# foreach (idx = 1:nrow(args_simulation)) %dopar% {
if (Sys.getenv("USER") == "kwiatkoe") {
  root  <-  "/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/onearm/code"
  setwd(root)
  #idx  <-  101
} else {
  args <- commandArgs(trailingOnly = TRUE)
  idx  <- as.numeric(args[1])
}

source("code_functions.R", local = TRUE)
source("args_model.R")
args_simulation <- read.csv(file="args_simulation.csv",header=TRUE,sep=",")

# using idx to identify model
p.range   <- args_simulation$p.range[idx]
freq.mntr <- args_simulation$freq.mntr[idx]
enr.shape <- args_simulation$enr.shape[idx]
out.mean  <- args_simulation$out.mean[idx]
if (args_simulation$skpt_spike[idx]==0){prior.nc.skpt <- skpt_prior_default()}
if (args_simulation$skpt_spike[idx]==1){prior.nc.skpt <- skpt_prior_custom(scale=0.75)}
if (args_simulation$enth_flat[idx]==0) {prior.nc.enth <- enth_prior_default()}
if (args_simulation$enth_flat[idx]==1) {prior.nc.enth <- enth_prior_custom(scale=1.5)}

source("code_main.R", local = TRUE)

Table1     <- data.frame(t(outer[,,]))
Table1$idx <- idx
write.csv(Table1,file=paste0("../output/Table1/",idx,"Table1.csv"))

Table2     <- data.frame(t(outer.p[,,]))
Table2$idx <- idx
write.csv(Table2,file=paste0("../output/Table2/",idx,"Table2.csv"))

Table3     <- data.frame(t(outer.p.agree[,,]))
Table3$idx <- idx
write.csv(Table3,file=paste0("../output/Table3/",idx,"Table3.csv"))
# }

# end_time  <- Sys.time()
# diff_time <- difftime(end_time, start_time, units = "auto")
# cat("Started  ", as.character(start_time), "\n",
#     "Finished ", as.character(end_time), "\n",
#     "Time difference of ", diff_time, " ", attr(diff_time, "units"), "\n",
#     "Used ", foreach::getDoParWorkers(), " cores\n",
#     "Used ", foreach::getDoParName(), " as backend\n",
#     sep = "")

