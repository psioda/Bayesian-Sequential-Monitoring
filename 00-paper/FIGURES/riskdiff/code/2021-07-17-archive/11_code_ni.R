rm(list=ls())

library(foreach)
library(doParallel)
library(pracma)
library(gnorm)

registerDoParallel(detectCores())
getDoParWorkers()
load(file = 'args_model.RData')

Table0 <- read.csv("../output/Table0_merged_071620.csv")
Table0 <- data.frame(Table0[Table0$eff.mix.prob==10,])

Table0.initial <- Table0[!duplicated(Table0[c("y1.PC.initial","y0.PC.initial","y1.IP.initial","y0.IP.initial")]),
                         c("y1.PC.initial","y0.PC.initial","y1.IP.initial","y0.IP.initial")]
names(Table0.initial) <- c("y1.PC","y0.PC","y1.IP","y0.IP")
Table0.final   <- Table0[!duplicated(Table0[c("y1.PC.final","y0.PC.final","y1.IP.final","y0.IP.final")]),
                         c("y1.PC.final","y0.PC.final","y1.IP.final","y0.IP.final")]
names(Table0.final) <- c("y1.PC","y0.PC","y1.IP","y0.IP")

Table0 <- rbind(Table0.initial, Table0.final)
Table0 <- Table0[!duplicated(Table0),c("y1.PC","y0.PC","y1.IP","y0.IP")]

Table0$box.ni.final <- NA

prior_dat_conflict <- function(y1.IP, y0.IP, y1.PC, y0.PC){
  ni.post.nc         <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  
  for (i in 1:(y1.IP + y0.IP + 1)){
    for (j in 1:(y1.PC + y0.PC + 1)){
      
      ni.post.1 <- function(x, y){ # for x > 0 (theta > 0)
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
            log(pgnorm(q = 1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
                  pgnorm(q = -1,  delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) +
            dgnorm(y,          mu,            ni.alpha0,    ni.beta0, log = TRUE) -
            log(pgnorm(q = 1-x,mu,            ni.alpha0,    ni.beta0) -
                  pgnorm(q = 0,   mu,            ni.alpha0,    ni.beta0))
        )
      }
      ni.post.2 <- function(x, y){ # for x < 0 (theta < 0)
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,           delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0, log = TRUE) -
            log(pgnorm(q = 1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0) -
                  pgnorm(q = -1, delta.ni.skpt, ni.rd.alpha0, ni.rd.beta0)) + 
            dgnorm(y,         mu,            ni.alpha0,    ni.beta0, log = TRUE) -
            log(pgnorm(q = 1, mu,            ni.alpha0,    ni.beta0) -
                  pgnorm(q = -x, mu,            ni.alpha0,    ni.beta0))
        )
      }
      ni.post.nc[i, j] <- integrate_debug(ni.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
        integrate_debug(ni.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
    }
  }
 
  # prior data conflict for non-informative prior
  box.ni <- sum(ni.post.nc[ni.post.nc <= ni.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with ni prior (should be 1): ", sum(ni.post.nc)))
  print(paste0("Non-informative prior compatibility: ", box.ni))
  return(box.ni)
}

start_time <- Sys.time()
summary(Table0$box.ni.final)

x <- foreach (i = 1:nrow(Table0), .combine='c') %dopar% {
  prior_dat_conflict(Table0$y1.IP[i],
                     Table0$y0.IP[i],
                     Table0$y1.PC[i],
                     Table0$y0.PC[i])
}

x.t         <- df(t(matrix(data = t(x), nrow = 4)))
names(x.t)  <- 
Table0$box.ni.final <- x

summary(Table0$box.ni.final)
end_time  <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
cat("Started  ", as.character(start_time), "\n",
    "Finished ", as.character(end_time), "\n",
    "Time difference of ", diff_time, " ", attr(diff_time, "units"), "\n",
    sep = "")