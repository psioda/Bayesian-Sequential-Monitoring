#############################################
#### Figure 6, Risk Diff Inference Plots ####
#############################################

rm(list = ls())
setwd("/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")

output_png <- TRUE
sig.fut    <- 0.975
sig.eff    <- 0.975

width.scale <- 9
if(output_png){
  png('../code/plots/updated power plot/figure6.png',
      width = 450*width.scale,
      height = 300*width.scale,
      pointsize=16,
      res=300)
}
par(mar=c(5.1 + 5, 4.1 + 2, 2.1, 2.1 + 0.5)) #c(bottom, left, top, right)

#par(mar=c(5.1+2, 4.1, 4.1, 2.1)) #c(bottom, left, top, right)
stretch1 <- -0.03# to add x-axis table under graph
stretch2 <- -0.03
# set initial plotting area
plot(NULL,
     type = 'l',
     xlim = c(-0.005,0.125),
     ylim = c(0,0.8),
     lwd  = 1,
     ylab = "Probability of Interim Stoppage for Efficacy",
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
# mtext(text=c(as.expression(bquote(theta))),side=1,line=1,at=stretch1,adj=0)
mtext(text="Risk Difference",side=1,line=2,at=0.05,adj=0)

# PLOT
x <- seq(0, 0.12, by = 0.03)

axis(1,
     las = 0,
     at = x,
     labels = format(x, nsmall = 2),
     line = 0)


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


# make loop
loop.list <- 101:106
delta.list <- seq(0, 0.25, by = 0.05)
for (i in 1:length(loop.list)){
  row <- i + 3
  temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51 & final$eff.mix.prob == loop.list[i], ]
  y <- temp$success
  lines(x,y)
  for (j in seq(1,length(temp$p.IP)#, by=6
  )){
    mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
                      " / ",
                      format(round(temp$ss.final[j],digits=1),nsmall=1)),
          side=1,line=row,at=x[j])
    # mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
    #                   " + ",
    #                   format(round(temp$ss.final[j]-
    #                                  temp$ss.initial[j],digits=1),nsmall=1),
    #                   " = ",
    #                   format(round(temp$ss.final[j],digits=1),nsmall=1)),
    #       side=1,line=row,at=x[j])
  }
  # text(temp$p.IP[j], temp$eff.prob.initial[j], row - 1)
  text(x = temp$p.IP[j] - 0.39 + 0.0025, 
       y = temp$success[j] , 
       labels = format(round(temp$success[j],digits=2),nsmall=2))
  if (i %in% c(1,2,3)){
  text(x = temp$p.IP[1] - 0.39 - 0.0025, 
       y = temp$success[1], 
       labels = format(round(temp$success[1],digits=2),nsmall=2))
  }
  if (i == 4){
    text(x = temp$p.IP[1] - 0.39 - 0.0025, 
         y = 0.04, 
         labels = format(round(temp$success[1],digits=2),nsmall=2))
  }
  if (i == 5){
    text(x = temp$p.IP[1] - 0.39 - 0.0025, 
         y = 0.015, 
         labels = format(round(temp$success[1],digits=2),nsmall=2))
  }
  if (i == 6){
    text(x = temp$p.IP[1] - 0.39 - 0.0025, 
         y = -0.01, 
         labels = format(round(temp$success[1],digits=2),nsmall=2))
  }
  text(x = 0.06, 
       y = temp$success[3], 
       labels = row - 3)
}

legend('topleft',
       title = "Skeptical Prior Weight Minimum",
       legend= c(as.expression(bquote("1: "*delta*"=0.00")),
                 as.expression(bquote("2: "*delta*"=0.05")),
                 as.expression(bquote("3: "*delta*"=0.10")),
                 as.expression(bquote("4: "*delta*"=0.15")),
                 as.expression(bquote("5: "*delta*"=0.20")),
                 as.expression(bquote("6: "*delta*"=0.25"))))

mtext(text="Expected Sample Size (Interim/Final)",
      side=1,line=3,at=stretch2,adj=0)

# mtext(text=as.expression(bquote(delta*"=0.00")),
#       side=1,line=2,at=stretch1,adj=0)
# mtext(text=as.expression(bquote(delta*"=0.05")),
#       side=1,line=3,at=stretch1,adj=0)
# mtext(text=as.expression(bquote(delta*"=0.10")),
#       side=1,line=4,at=stretch1,adj=0)
# mtext(text=as.expression(bquote(delta*"=0.15")),
#       side=1,line=5,at=stretch1,adj=0)
# mtext(text=as.expression(bquote(delta*"=0.20")),
#       side=1,line=6,at=stretch1,adj=0)
# mtext(text=as.expression(bquote(delta*"=0.25")),
#       side=1,line=7,at=stretch1,adj=0)

mtext(text=as.expression(bquote("1: "*delta*"=0.00")),
      side=1,line=4,at=stretch1,adj=0)
mtext(text=as.expression(bquote("2: "*delta*"=0.05")),
      side=1,line=5,at=stretch1,adj=0)
mtext(text=as.expression(bquote("3: "*delta*"=0.10")),
      side=1,line=6,at=stretch1,adj=0)
mtext(text=as.expression(bquote("4: "*delta*"=0.15")),
      side=1,line=7,at=stretch1,adj=0)
mtext(text=as.expression(bquote("5: "*delta*"=0.20")),
      side=1,line=8,at=stretch1,adj=0)
mtext(text=as.expression(bquote("6: "*delta*"=0.25")),
      side=1,line=9,at=stretch1,adj=0)

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
# 
# # LOAD & RECODE FIXED WEIGHTS
# dat <- read.csv("Table0_merged-2021-07-19-vFixed.csv")
# # head(dat)
# # table(dat$p.IP)
# dat$success <- dat$eff.prob.initial >= 0.975
# dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
# final <- aggregate(dat, list(dat$p.IP), mean)
# library(dplyr)
# final <- dat %>%
#   group_by(p.IP, eff.mix.prob) %>% 
#   summarise_each(funs(mean))
# # final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
# final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
# final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
# final <- data.frame(final)
# 
# 
# 
# #FIRST LINE
# 
# temp <- final[final$eff.mix.prob == 1 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
# y <- temp$success
# row = 6
# lines(x,y)
# 
# for (j in seq(1,length(temp$p.IP)#, by=6
# )){
#   mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
#                     " + ",
#                     format(round(temp$ss.final[j]-
#                                    temp$ss.initial[j],digits=1),nsmall=1),
#                     " = ",
#                     format(round(temp$ss.final[j],digits=1),nsmall=1)),
#         side=1,line=row,at=x[j])
# }
# text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
# mtext(text=paste0(row - 1),side=1,line=row,at=stretch1,adj=0)
# 
# 
# ## SECOND LINE
# row <- 5
# temp <- final[final$eff.mix.prob == 0.5 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
# y <- temp$success
# lines(x,y)
# for (j in seq(1,length(temp$p.IP)#, by=6
# )){
#   mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
#                     " + ",
#                     format(round(temp$ss.final[j]-
#                                    temp$ss.initial[j],digits=1),nsmall=1),
#                     " = ",
#                     format(round(temp$ss.final[j],digits=1),nsmall=1)),
#         side=1,line=row,at=x[j])
# }
# text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
# mtext(text=paste0(row - 1),side=1,line=row,at=stretch1,adj=0)
# 
# ## THIRD LINE
# row <- 2
# temp <- final[final$eff.mix.prob == 0 & final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
# y <- temp$success
# lines(x,y)
# for (j in seq(1,length(temp$p.IP)#, by=6
# )){
#   mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
#                     " + ",
#                     format(round(temp$ss.final[j]-
#                                    temp$ss.initial[j],digits=1),nsmall=1),
#                     " = ",
#                     format(round(temp$ss.final[j],digits=1),nsmall=1)),
#         side=1,line=row,at=x[j])
# }
# text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
# mtext(text=paste0(row - 1),side=1,line=row,at=stretch1,adj=0)
# 
# # FORTH LINE
# # LOAD & RECODE v103
# dat <- read.csv("Table0_merged-2021-07-19-v103.csv")
# # head(dat)
# # table(dat$p.IP)
# dat$success <- dat$eff.prob.initial >= 0.975
# dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
# final <- aggregate(dat, list(dat$p.IP), mean)
# library(dplyr)
# final <- dat %>%
#   group_by(p.IP, eff.mix.prob) %>% 
#   summarise_each(funs(mean))
# # final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
# final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
# final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
# final <- data.frame(final)
# 
# temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
# 
# row <- 4
# y <- temp$success
# lines(x,y)
# for (j in seq(1,length(temp$p.IP)#, by=6
# )){
#   mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
#                     " + ",
#                     format(round(temp$ss.final[j]-
#                                    temp$ss.initial[j],digits=1),nsmall=1),
#                     " = ",
#                     format(round(temp$ss.final[j],digits=1),nsmall=1)),
#         side=1,line=row,at=x[j])
# }
# text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
# mtext(text=paste0(row - 1),side=1,line=row,at=stretch1,adj=0)
# 
# # FIFTH LINE
# # LOAD & RECODE v115
# dat <- read.csv("Table0_merged_2021-07-19-v115.csv")
# # head(dat)
# # table(dat$p.IP)
# dat$success <- dat$eff.prob.initial >= 0.975
# dat$pm.rd.f2 <- (dat$y1.IP.final / (dat$y1.IP.final + dat$y0.IP.final)) - (dat$y1.PC.final / (dat$y1.PC.final + dat$y0.PC.final))
# final <- aggregate(dat, list(dat$p.IP), mean)
# library(dplyr)
# final <- dat %>%
#   group_by(p.IP, eff.mix.prob) %>% 
#   summarise_each(funs(mean))
# # final$mean.ss <- final$y1.IP.f + final$y0.IP.f + final$y1.PC.f + final$y0.PC.f # avg ss
# final$ss.initial <- final$y1.IP.initial + final$y0.IP.initial + final$y1.PC.initial + final$y0.PC.initial # avg ss
# final$ss.final <- final$y1.IP.final + final$y0.IP.final + final$y1.PC.final + final$y0.PC.final # avg ss
# final <- data.frame(final)
# 
# temp <- final[final$p.IP >= 0.39 & final$p.IP <= 0.51, ]
# 
# row <- 3
# y <- temp$success
# lines(x,y)
# for (j in seq(1,length(temp$p.IP)#, by=6
# )){
#   mtext(text=paste0(format(round(temp$ss.initial[j],digits=1),nsmall=1),
#                     " + ",
#                     format(round(temp$ss.final[j]-
#                                    temp$ss.initial[j],digits=1),nsmall=1),
#                     " = ",
#                     format(round(temp$ss.final[j],digits=1),nsmall=1)),
#         side=1,line=row,at=x[j])
# }
# text(temp$p.IP[j], temp$eff.mon.initial[j], row - 1)
# mtext(text=paste0(row - 1),side=1,line=row,at=stretch1,adj=0)
# 

mtext("(B)", side = 2, line = 3, at = 0.8, las = 1)

if(output_png){dev.off()}