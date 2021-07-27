rm(list=ls())
library(pracma)
library(gnorm)
library(foreach)
library(doParallel)

registerDoParallel(detectCores())
getDoParWorkers()

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/")
load(file = '../args_model.RData') # loads all model information include prior parameters AND SETS SEED

# create data matrix
y1.IP        <- seq(0, 58)
y0.IP        <- 58 - y1.IP
Table0       <- data.frame(y1.IP, y0.IP)
Table0$y1.PC <- 16
Table0$y0.PC <- 26

prior_dat_conflict <- function(y1.IP, y0.IP, y1.PC, y0.PC){
  
  skpt.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  enth.post.nc       <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)
  ni.post.nc         <- matrix(NA, 
                               nrow = y1.IP + y0.IP + 1, 
                               ncol = y1.PC + y0.PC + 1)

  for (i in 1:(y1.IP + y0.IP + 1)){
    for (j in 1:(y1.PC + y0.PC + 1)){
      skpt.post.1   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
                  pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
            dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
            log(pgnorm(q = 1-x,  mu, skpt.alpha0, skpt.beta0) - 
                  pgnorm(q = 0,  mu, skpt.alpha0, skpt.beta0))
        )
      }
      skpt.post.2   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.skpt, skpt.rd.alpha0, skpt.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.skpt, skpt.rd.alpha0, skpt.rd.beta0) - 
                  pgnorm(q = -1, delta.skpt, skpt.rd.alpha0, skpt.rd.beta0)) +
            dgnorm(y,            mu, skpt.alpha0, skpt.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    mu, skpt.alpha0, skpt.beta0) - 
                  pgnorm(q = -x, mu, skpt.alpha0, skpt.beta0))
        )
      }
      skpt.post.nc[i, j] <- integrate_debug(skpt.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
        integrate_debug(skpt.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
      enth.post.1   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
                  pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
            dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
            log(pgnorm(q = 1-x,  mu, enth.alpha0, enth.beta0) - 
                  pgnorm(q = 0,  mu, enth.alpha0, enth.beta0))
        )
      }
      enth.post.2   <- function(x, y){
        exp(
          dbinom(i-1,   y1.IP + y0.IP, x + y, log = TRUE) +
            dbinom(j-1, y1.PC + y0.PC, y, log = TRUE) + 
            dgnorm(x,            delta.enth, enth.rd.alpha0, enth.rd.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    delta.enth, enth.rd.alpha0, enth.rd.beta0) - 
                  pgnorm(q = -1, delta.enth, enth.rd.alpha0, enth.rd.beta0)) +
            dgnorm(y,            mu, enth.alpha0, enth.beta0, log = TRUE) - 
            log(pgnorm(q = 1,    mu, enth.alpha0, enth.beta0) - 
                  pgnorm(q = -x, mu, enth.alpha0, enth.beta0))
        )
      }
      enth.post.nc[i, j] <- integrate_debug(enth.post.1, xmin = 0,  xmax = 1, ymin = 0, ymax = function(x) 1 - x) +
        integrate_debug(enth.post.2, xmin = -1, xmax = 0, ymin = function(x) -x, ymax = 1)
      
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
  # prior data conflict for skeptical prior
  box.skpt <- sum(skpt.post.nc[skpt.post.nc <= skpt.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(skpt.post.nc)))
  print(paste0("Skeptical prior compatibility: ", box.skpt))
  
  # prior data conflict for enthusiastic prior
  box.enth <- sum(enth.post.nc[enth.post.nc <= enth.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with enth prior (should be 1): ", sum(enth.post.nc)))
  print(paste0("Enthuastic prior compatibility: ", box.enth))
  
  # prior data conflict for non-informative prior
  box.ni <- sum(ni.post.nc[ni.post.nc <= ni.post.nc[y1.IP + 1, y1.PC + 1]])
  print(paste0("Sum marginal prob of data with ni prior (should be 1): ", sum(ni.post.nc)))
  print(paste0("Non-informative prior compatibility: ", box.ni))
  
  # compute SKEPTICAL COMPONENT mixing weight
  eff.mix.prob <- 1 - max(box.enth - box.skpt, 0)
  return(cbind(eff.mix.prob, box.skpt, box.enth, box.ni))
}

start_time <- Sys.time()

x <- foreach (i = 1:nrow(Table0), .combine='c') %dopar% {
  prior_dat_conflict(Table0$y1.IP[i],
                     Table0$y0.IP[i],
                     Table0$y1.PC[i],
                     Table0$y0.PC[i])
}

x.t         <- data.frame(t(matrix(data = x, nrow = 4)))
names(x.t)  <- c("eff.mix.prob", "box.skpt", "box.enth", "box.ni")

end_time  <- Sys.time()
diff_time <- difftime(end_time, start_time, units = "auto")
cat("Started  ", as.character(start_time), "\n",
    "Finished ", as.character(end_time), "\n",
    "Time difference of ", diff_time, " ", attr(diff_time, "units"), "\n",
    sep = "")

Table0 <- cbind(Table0, x.t)
Table0$risk.diff <- Table0$y1.IP/(Table0$y1.IP + Table0$y0.IP) - Table0$y1.PC/(Table0$y1.PC + Table0$y0.PC)

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
# 
# output_png <- TRUE
# width.scale <- 7
# 
# if(output_png){png('3 part compatibility/3-part-compatibility-1.png',
#                    width = 450*width.scale, 
#                    height = 300*width.scale,
#                    pointsize=16,
#                    res=300)}
# 
# x <- c(Table0$risk.diff)
# y <- c(Table0$box.skpt)
# plot(x, y, xlab = "Observed Response Difference", ylab = "Box's p-value", pch = 19, cex = 0.25, type = 'l',
#      lty='longdash', yaxt='n')
# # model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
# # myPredict <- predict( model ) 
# # ix <- sort(x,index.return=T)$ix
# # lines(x[ix], myPredict[ix], lwd=2 )  
# 
# y <- c(Table0$box.enth)
# lines(x,y, lty='dotted')
# # model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4) + I(x^5) + I(x^6) + I(x^7) + I(x^8) + I(x^9) + I(x^10))
# # myPredict <- predict( model ) 
# # ix <- sort(x,index.return=T)$ix
# # lines(x[ix], myPredict[ix], col=2, lwd=2 )  
# 
# # y <- c(Table0$box.ni)
# # lines(x,y, lty="solid")
# # model <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
# # myPredict <- predict( model ) 
# # ix <- sort(x,index.return=T)$ix
# # lines(x[ix], myPredict[ix], col="blue", lwd=2 )  
# 
# abline(v=0, col="grey")
# abline(v=0.12, col="grey")
# axis(2,
#      las = 2,
#      at = seq(0, 1, by = 0.1),
#      labels = format(seq(0, 1, by = 0.1), nsmall = 1))
# abline(h = seq(0, 1, by = 0.1),
#        col = 'grey')
# 
# legend("topright",#text.width = 0.05,
#        legend = c("Skeptical",
#                   "Enthusiastic"),
#        lty = c('longdash','dotted'))
# 
# if(output_png){dev.off()}

############################################################################################################################################
############################################################################################################################################
############################################################################################################################################

# delta.list <- seq(0, 0.25, by = 0.05)
# beta.list <- c(1, 1.3220, 1.7370, 2.3220, 3.3220)
# adapt.mat <- expand.grid(delta.list, beta.list)
# adapt.mat$id <- 101:130
# 
# if (eff.mix.prob > 100){
#   delta.a <- adapt.mat[adapt.mat$id == eff.mix.prob, "Var1"]
#   beta.a <- adapt.mat[adapt.mat$id == eff.mix.prob, "Var2"]
#   
#   prior_data_conflict_result <- prior_dat_conflict(y1.IP, y0.IP, y1.PC, y0.PC)
#   
#   box.skpt                   <- prior_data_conflict_result[,"box.skpt"]
#   box.enth                   <- prior_data_conflict_result[,"box.enth"]
#   box.ni                     <- prior_data_conflict_result[,"box.ni"]
#   BoxPE                      <- box.enth
#   eff.mix.prob               <- 1 - (1 - delta.a)*pbeta(BoxPE, 1, beta.a) # amount assigned to skeptical prior
# }

output_png <- TRUE
width.scale <- 7
# default
delta.a <- 0
beta.a <- 1
Table0$eff.mix.prob1               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)
delta.a <- 0.05
Table0$eff.mix.prob2               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)
delta.a <- 0.1
Table0$eff.mix.prob3               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)
delta.a <- 0.15
Table0$eff.mix.prob4               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)
delta.a <- 0.2
Table0$eff.mix.prob5               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)
delta.a <- 0.25
Table0$eff.mix.prob6               <- (1 - delta.a)*pbeta(Table0$box.enth, 1, beta.a)

if(output_png){png('3 part compatibility/3-part-compatibility-2.png',
                   width = 450*width.scale, 
                   height = 300*width.scale,
                   pointsize = 16,
                   res = 300)}

# omega.skpt <- Table0$box.skpt
# omega.enth <- Table0$box.enth
# omega.ni   <- Table0$box.ni
# omega.ni   <- pmax(omega.ni-pmax(omega.skpt,omega.enth),0)
# omega.sum  <- omega.skpt+omega.enth+omega.ni
# omega.skpt <- omega.skpt/omega.sum
# omega.enth <- omega.enth/omega.sum
# omega.ni   <- omega.ni/omega.sum

x <- c(Table0$risk.diff)
y <- Table0$eff.mix.prob1
plot(x, y, 
     # pch = 19, 
     # cex = 0.25, 
     ylim = c(0,1), 
     xlab = "Observed Response Difference", 
     ylab = "Enthusiastic Mixture Weight", type = 'l',
     lty=1, 
     xaxt = 'n',
     yaxt='n')

axis(1,
     las = 0,
     at = seq(-4, 0.6, by = 0.05),
     labels = format(seq(-4, 0.6, by = 0.05), nsmall = 2),
     line = 0)
axis(2,
     las = 2,
     at = seq(0, 1, by = 0.1),
     labels = format(seq(0, 1, by = 0.1), nsmall = 1))
abline(h = seq(0, 1, by = 0.1),
       col = 'grey')
# abline(v=0, col="grey")
# abline(v=0.12, col="grey")


lines(x, Table0$eff.mix.prob1)
text(0.125, 1, labels = 1)
lines(x, Table0$eff.mix.prob2)
text(0.125, 0.95, labels = 2)
lines(x, Table0$eff.mix.prob3)
text(0.125, 0.9, labels = 3)
lines(x, Table0$eff.mix.prob4)
text(0.125, 0.85, labels = 4)
lines(x, Table0$eff.mix.prob5)
text(0.125, 0.8, labels = 5)
lines(x, Table0$eff.mix.prob6)
text(0.125, 0.75, labels = 6)





legend('topleft',
       title = "Skeptical Prior Weight Minimum",
       legend= c(as.expression(bquote("1: "*delta*"=0.00")),
                 as.expression(bquote("2: "*delta*"=0.05")),
                 as.expression(bquote("3: "*delta*"=0.10")),
                 as.expression(bquote("4: "*delta*"=0.15")),
                 as.expression(bquote("5: "*delta*"=0.20")),
                 as.expression(bquote("6: "*delta*"=0.25"))))

# legend("top",#text.width = 0.05,
#        legend = c("Skeptical",
#                   "Enthusiastic"),
#        lty = c('solid','longdash','dotted'))
mtext("(A)", side = 2, line = 3, at = 1, las = 1)

if(output_png){dev.off()}