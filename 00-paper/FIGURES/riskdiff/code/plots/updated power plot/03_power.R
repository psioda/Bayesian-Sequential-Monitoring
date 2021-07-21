#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################

rm(list = ls())
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")

output_png <- FALSE
sig.fut    <- 0.975
sig.eff    <- 0.975

width.scale <- 7
if(output_png){
  png('figure6.png',
      width = 450*width.scale,
      height = 300*width.scale,
      pointsize=16,
      res=300)
}
par(mar=c(5.1 + 2, 4.1 + 0.5, 2.1, 2.1 + 0.5)) #c(bottom, left, top, right)

#par(mar=c(5.1+2, 4.1, 4.1, 2.1)) #c(bottom, left, top, right)
stretch <- -0.05# to add x-axis table under graph

# set initial plotting area
plot(NULL,
     type = 'l',
     xlim = c(0,0.12),
     ylim = c(0,1),
     lwd  = 1,
     ylab = "Probability",
     xlab = "", #Treatment Response Probability",
     main = "",
     xaxt = "n",
     yaxt = "n")
abline(h = seq(0, 1, by = 0.1),
       col = 'grey')
axis(2,
     las = 2,
     at = seq(0, 1, by = 0.1),
     labels = format(seq(0, 1, by = 0.1), nsmall = 1))
mtext(text=c(as.expression(bquote(theta))),side=1,line=1,at=stretch,adj=0)
legend('topleft',
       legend= c("1: 100% Enthusiastic",
                 as.expression(bquote("2: Adaptive "*delta*"=0.1, "*beta*"=0.7")),
                 as.expression(bquote("3: Adaptive "*delta*"=0.1, "*beta*"=0.5")),
                 "4: 50% Enthusiastic",
                 "5: 100% Skeptical"))
# PLOT
x <- seq(0, 0.12, by = 0.03)

axis(1,
     las = 0,
     at = x,
     labels = format(x, nsmall = 2))


################################################################################################################################################

# LOAD & RECODE v101-106
dat <- read.csv("Table0_merged_2021-07-20-v101-106.csv")
# head(dat)
# table(dat$p.IP)
dat$success <- dat$eff.prob.initial >= 0.975
dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
final <- aggregate(dat, list(dat$p.IP), mean)
library(dplyr)
final <- dat %>%
  group_by(p.IP, eff.mix.prob) %>% 
  summarise_each(funs(mean))
# final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
final <- data.frame(final)

temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51 & final$eff.mix.prob == 101, ]
row <- 2
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51 & final$eff.mix.prob == 101, ]
row <- 2
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

# make loop
loop.list <- 101:106
for (i in 1:length(loop.list)){
  row <- i + 1
  temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51 & final$eff.mix.prob == loop.list[i], ]
  y <- temp$success
  lines(x,y)
  for (j in seq(1,length(temp$p.IP)#, by=6
  )){
    mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                      " + ",
                      format(round(temp$ss.final[j]-
                                     temp$ss.initial[j],digits=1),nsmall=1),
                      " = ",
                      format(round(temp$ss.final[j],digits=1),nsmall=1)),
          side=1,line=row,at=x[j])
  }
  text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
  mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)
}




# # plot all fixed weight priors
# Table1        <- read.csv(file = "../../output/table1031920.csv", header = T)
# #Table1_ref    <- read.csv(file = "../../output/output062120/Table1/1Table1.csv")
# #names(Table1) <- names(Table1_ref)
# args          <- read.csv(file = "../../output/args_simulation031920.csv", header = F)
# args_ref      <- read.csv(file = "../args_simulation.csv")
# names(args)   <- names(args_ref)
# combined1     <- merge(args, Table1, by.x = "X", by.y = "idx")
# figure3       <- combined1

# Table0                 <- read.csv("../../output/Table0_merged_071620.csv")
# Table0$ss.final        <- Table0$y0.IP.final + 
#   Table0$y0.PC.final +
#   Table0$y1.IP.final +
#   Table0$y1.PC.final
# Table0$ss.initial        <- Table0$y0.IP.initial + 
#   Table0$y0.PC.initial +
#   Table0$y1.IP.initial +
#   Table0$y1.PC.initial
# Table0$eff.mon.initial <- Table0$eff.prob.initial >= sig.eff
# figure3         <- aggregate(x   = Table0,
#                              by  = list(Table0$p.IP, Table0$p.PC, Table0$eff.mix.prob),
#                              FUN = mean)

## may need to change eff.mix.prob to eff.mix.prob.x

# LOAD & RECODE FIXED WEIGHTS
dat <- read.csv("Table0_merged-2021-07-19-vFixed.csv")
# head(dat)
# table(dat$p.IP)
dat$success <- dat$eff.prob.initial >= 0.975
dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
final <- aggregate(dat, list(dat$p.IP), mean)
library(dplyr)
final <- dat %>%
  group_by(p.IP, eff.mix.prob) %>% 
  summarise_each(funs(mean))
# final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
final <- data.frame(final)



#FIRST LINE

temp <- final[final$eff.mix.prob == 1 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
y <- temp$success
row = 6
lines(x,y)

for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)


## SECOND LINE
row <- 5
temp <- final[final$eff.mix.prob == 0.5 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

## THIRD LINE
row <- 2
temp <- final[final$eff.mix.prob == 0 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

# FORTH LINE
# LOAD & RECODE v103
dat <- read.csv("Table0_merged-2021-07-19-v103.csv")
# head(dat)
# table(dat$p.IP)
dat$success <- dat$eff.prob.initial >= 0.975
dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
final <- aggregate(dat, list(dat$p.IP), mean)
library(dplyr)
final <- dat %>%
  group_by(p.IP, eff.mix.prob) %>% 
  summarise_each(funs(mean))
# final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
final <- data.frame(final)

temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51, ]

row <- 4
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

# FIFTH LINE
# LOAD & RECODE v115
dat <- read.csv("Table0_merged_2021-07-19-v115.csv")
# head(dat)
# table(dat$p.IP)
dat$success <- dat$eff.prob.initial >= 0.975
dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
final <- aggregate(dat, list(dat$p.IP), mean)
library(dplyr)
final <- dat %>%
  group_by(p.IP, eff.mix.prob) %>% 
  summarise_each(funs(mean))
# final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
final <- data.frame(final)

temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51, ]

row <- 3
y <- temp$success
lines(x,y)
for (j in seq(1,length(temp$p.IP)#, by=6
)){
  mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(temp$ss.final[j]-
                                   temp$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(temp$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=x[j])
}
text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
mtext(text=paste0(row - 1),side=1,line=row,at=stretch,adj=0)

if(output_png){dev.off()}

