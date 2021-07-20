setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")
dat <- read.csv("Table0_merged-2021-07-19-vFixed.csv")
head(dat)
table(dat$p.IP)
dat$success <- dat$eff.prob.initial >= 0.975
dat$successfinal <- dat$eff.prob.final >= 0.975
dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
dat <- dat[dat$p.IP %in% c(0.39, 0.42, 0.45, 0.48, 0.51) & dat$eff.mix.prob == 0, ]
table(dat$p.IP)
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots/evidence-decrease")

output_png <- TRUE
sig.eff    <- 0.975

width.scale <- 7
if(output_png){
  png('figure6a.png',
      width = 450*width.scale,
      height = 300*width.scale,
      pointsize=16,
      res=300)
}

table(dat$success, dat$successfinal)
reduced_mat <- dat[, c("p.IP", "eff.prob.final")][dat[, "eff.prob.initial"] > sig.eff & dat[, "eff.prob.final"] < sig.eff, ]

# quantiles
q <- c(0.1, 0.25, 0.5, 0.75, 0.9)
quants <- reduced_mat %>%
  group_by(p.IP) %>%
  summarize(quant10 = quantile(eff.prob.final, probs = q[1]), 
            quant25 = quantile(eff.prob.final, probs = q[2]),
            quant50 = quantile(eff.prob.final, probs = q[3]), 
            quant75 = quantile(eff.prob.final, probs = q[4]),
            quant90 = quantile(eff.prob.final, probs = q[5]))
quants <- data.frame(quants)

# % AGR
for(i in 1:nrow(dat)){
  if(dat$success[i] == 1 & dat$successfinal[i] == 1){
    dat$agree[i] <- 1
  } else if (dat$success[i] == 1 & dat$successfinal[i] == 0){
    dat$agree[i] <- 0
  } else {
    dat$agree[i] <- NA
  }
}

agr <- aggregate(dat$agree, by=list(dat$p.IP), FUN=mean, na.rm = TRUE)

table(dat$successfinal[dat$success == 1])

plot(quants$p.IP, quants$quant50,
     type='l',
     ylim=c(0.96, 0.975),
     ylab="",
     xlab="",
     main="",
     axes=FALSE)
box()

axis(2,las=2,at=seq(0.96,1,by=0.005),labels=format(seq(0.96,1,by=0.005),nsmall=1))
title(ylab="Final Efficacy Criteria", line=3, #cex.lab=1.2
)
abline(h=seq(0.96,1,by=0.005),col='grey')
lines(quants$p.IP,quants$quant10)
lines(quants$p.IP,quants$quant25)
lines(quants$p.IP,quants$quant75)
lines(quants$p.IP,quants$quant90)
lines(c(p.skpt,p.enth),rep(0.975,2))

polygon(x=c(quants$p.IP,rev(quants$p.IP),quants$p.IP[1]),
        y=c(quants$quant10,rev(quants$quant25),quants$quant10[1]),
        col='lightgrey')

polygon(x=c(quants$p.IP,rev(quants$p.IP),quants$p.IP[1]),
        y=c(quants$quant90,rev(quants$quant75),quants$quant90[1]),
        col='lightgrey')

polygon(x=c(quants$p.IP,rev(quants$p.IP),quants$p.IP[1]),
        y=c(quants$quant75,rev(quants$quant25),quants$quant75[1]),
        col='darkgrey')

polygon(x=c(quants$p.IP,rev(quants$p.IP),quants$p.IP[1]),
        y=c(quants$quant90,rep(0.975,length(quants$p.IP)),quants$quant90[1]),
        col='black',density=15,angle=45,lty=1)

lines(quants$p.IP,quants$quant50)
## SUBSET DATA ##
axis(1,las=0,at=quants$p.IP[seq(1,length(quants$p.IP))],
     labels=format(quants$p.IP[seq(1,length(quants$p.IP))] - 0.39,nsmall=2))
text(quants$p.IP,
     quants$quant50,
     labels=format(round(quants$quant50,digits=3),nsmall=3),
     pos=1)
points(quants$p.IP,quants$quant50,pch=20)

colnames(agr) <- c("p.IP", "agr")
for (j in seq(1,length(agr$p.IP))){

  mtext(text=paste0(format(round(agr$agr[j],digits=3),nsmall=3)),
        side=1,
        line=2,
        at=agr$p.IP[j])
}

lines(x=c(p.skpt,p.enth),rep(0.975,2))

stretch<-0.375
mtext(bquote(theta),side=1,line=1,at=stretch)
mtext(text="% AGR",side=1,line=2,at=stretch)
leg <-legend("bottomright",
       c("Median","P10-P90","P25-75","P90-P100"), 
       fill=c("black","lightgray","darkgray","black"),
       density=c(0,NA,NA,10),
       pch=c(20,NA,NA,NA),
       lty=c(1,NA,NA,NA),
       border=c(NA,'black','black','black'),
       x.intersp=c(0,0,0,0))

if(output_png){dev.off()}
