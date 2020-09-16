output_png <- TRUE

enth_prior_custom <- function(scale){
  
  mu0.enth    <- p.intr
  sigma0.seq  <- seq(.01, 0.5, by=0.01)
  lambda0.seq <- seq(1, 7, by=0.1)
  result1     <- matrix(NA, nrow=length(sigma0.seq), ncol=length(lambda0.seq))
  result2     <- matrix(NA, nrow=length(sigma0.seq), ncol=length(lambda0.seq))
  
  for (i in 1:length(sigma0.seq)){
    for (j in 1:length(lambda0.seq)){
      
      sigma0.enth  <- sigma0.seq[i]
      lambda0.enth <- lambda0.seq[j]
      
      prior.enth <- function(x){
        exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
      }
      
      nc.enth <- tryCatch(integrate(prior.enth, lower=-Inf, upper=Inf)[[1]], error = function(e) NA)
      
      if (nc.enth != 0 & !is.na(nc.enth)){
        prior.nc.enth <- function(x){
          exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
        }
        result1[i, j] <- integrate(prior.nc.enth, lower=-Inf,                     upper=mu0.enth-(p.enth-p.skpt))[[1]]    # change 7/16/20
        result2[i, j] <- integrate(prior.nc.enth, lower=mu0.enth-(p.enth-p.skpt), upper=p.skpt)[[1]]                      # change 7/16/20
      }
      
    }
  }
  result3 <- abs(result1-tail.enth)+abs(result2-(pnorm(qnorm(tail.enth)/2)-tail.enth)*scale)
  index   <- which(result3 == min(result3, na.rm = T),  arr.ind = TRUE)
  
  i <- index[1]
  j <- index[2]
  
  sigma0.enth  <- sigma0.seq[i]
  lambda0.enth <- lambda0.seq[j]
  
  prior.enth <- function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)
  }
  nc.enth <- integrate(prior.enth, lower=0+epsilon, upper=1-epsilon)[[1]]
  prior.nc.enth <- function(x){
    exp(-(abs(x-mu0.enth)/sigma0.enth)^lambda0.enth)/nc.enth
  }
  
  print(paste0("mu: ", mu0.enth, ",  sigma: ", sigma0.enth, ",  lambda: ", lambda0.enth))
  print(paste0("Tail area: ", result1[i, j]))
  print(paste0("Half-width area: ", result2[i, j]))
  
  assign("mu0.enth", mu0.enth, envir = .GlobalEnv)
  assign("sigma0.enth", sigma0.enth, envir = .GlobalEnv)
  assign("lambda0.enth", lambda0.enth, envir = .GlobalEnv)
  assign("tail.enth.actual", result1[i, j], envir = .GlobalEnv)
  
  return(prior.nc.enth)
}

p.skpt    <- 0.40     # response rate for skeptic,  enthusiast,  futility
p.enth    <- 0.67
p.intr    <- (p.skpt+p.enth)/2
tail.skpt <- 0.025  # tail probabilities for priors
tail.enth <- 0.025
cred.tail <- 0.05
sig.fut   <- 0.975
sig.eff   <- 0.975
epsilon   <- 0 # used to stop numerical error from integration
max.ss    <- 112
reps      <- 100000
mu0.skpt  <- p.skpt
mu0.enth  <- p.intr

width.scale <- 5
if(output_png){png('figure1d.png', width = 300*width.scale,  height = 300*width.scale, pointsize=12, res=300)}
scale         <-  1.5
prior.nc.enth <- enth_prior_custom(scale=scale)

xmin  <-  p.enth - 0.5
xmax  <-  p.enth + 0.5
ymax  <-  6
x <- seq(xmin, 
         xmax, 
         by=0.005)
plot(x, prior.nc.enth(x), type="l", 
     xlab="", 
     ylab="", 
     main="", 
     xaxt="n", 
     yaxt="n", 
     xlim=c(xmin, xmax), 
     ylim=c(0, ymax)) # 20-01-02
#axis(2, at=c(0, 1, 2, 3), labels=c(0, 1, 2, 3))
axis(1, at=c(p.enth, p.skpt#, (p.enth+p.skpt)/2, (3*p.skpt-p.enth)/2
             ), 
     labels=c(as.expression(bquote(theta[1])), 
              as.expression(bquote(theta[0])) 
              #as.expression(bquote((theta[0]+theta[1])/2)),
              #as.expression(bquote((3*theta[0]-theta[1])/2))
              ))
title(ylab="Density Value",  line=1)
title(xlab="Response Probability", line=2)
#title(xlab="Density Value", line=2)

# polygon(c(x[x<=mu0.enth-(p.enth-p.skpt)], mu0.enth-(p.enth-p.skpt)), 
#         c(prior.nc.enth(x)[x<=mu0.enth-(p.enth-p.skpt)], 0), col="black")

# polygon(c(mu0.enth-(p.enth-p.skpt), x[x>=mu0.enth-(p.enth-p.skpt) & x<=p.skpt], p.skpt), 
#         c(0, prior.nc.enth(x)[x>=mu0.enth-(p.enth-p.skpt) & x<=p.skpt], 0), col="lightgrey")

segments(x0=p.enth, y0=0, y1=prior.nc.enth(p.enth))
segments(x0=p.skpt, y0=0, y1=prior.nc.enth(p.skpt))

#segments(x0=p.intr, y0=0, y1=prior.nc.enth(p.intr))

legend("top", 
       legend= c(
         #as.expression(bquote(mode(theta) == (theta[0]+theta[1])/2)), 
         as.expression(bquote(pi[E](theta[0])==pi[E](theta[1]))),
         as.expression(bquote(pi[E](theta)%~~%"c for "*theta %in% (theta[0]*", "*theta[1])))
         #as.expression(bquote(P(theta< (3*theta[0]-theta[1])/2)==.(tail.enth))), 
         #as.expression(bquote(P(theta %in% ((3*theta[0]-theta[1])/2*", "*theta[0])==.(round((pnorm(qnorm(tail.enth)/2)-tail.enth)*scale, 3)))))#, 
         #as.expression(bquote(GN(mu==theta[1], alpha==.(sigma0.enth), beta==.(lambda0.enth))))
       ))

mtext("(D)", side=2, line=1, at=ymax, las=1)
if(output_png){dev.off()}