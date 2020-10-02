###########################################################
#### Figure XXXX, Two Arm Sequential Design Properties ####
###########################################################
rm(list = ls())

output_png <- TRUE
root       <- "/Users/kwiatkoe/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/code/plots"
setwd(root)

########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

Table0    <- read.csv("../Table0_merged.csv") # CHANGE HERE
Table0    <- within(Table0, rm(X,X.1,
                               box.ni.initial,box.ni.final,
                               box.skpt.initial,box.skpt.final,
                               box.enth.initial,box.enth.final)) # get rid of non-numeric columns
#Table0    <- read.csv("../../output/Table0_merged_071620.csv") # CHANGE HERE
#Table0    <- Table0[Table0$eff.mix.prob==10,]                  # CHANGE HERE
#Table0    <- within(Table0, rm(X,X.1,box.ni.initial,box.ni.final))


## RUN TO INCLUDE PM & CP
# ni_unique <- read.csv("../../output/Table0_ni_unique_pm_cp_072320.csv")
# ni_unique <- within(ni_unique, rm(X))
# names(ni_unique) <- paste(names(ni_unique),".final", sep="")
# Table0           <- merge(Table0, 
#                           ni_unique, 
#                           by=c("y1.IP.final","y0.IP.final","y1.PC.final","y0.PC.final",
#                                "box.skpt.final","box.enth.final"))
# 
# 
# names(ni_unique) <- gsub("final", "initial", names(ni_unique))
# Table0           <- merge(Table0, 
#                           ni_unique, 
#                           by=c("y1.IP.initial","y0.IP.initial","y1.PC.initial","y0.PC.initial",
#                                "box.skpt.initial","box.enth.initial"))
# # coverage probability
# Table0$cp_rd.initial <- (Table0$lower_rd.initial < Table0$pm_rd.initial & Table0$pm_rd.initial < Table0$upper_rd.initial)
# Table0$cp_rd.final   <- (Table0$lower_rd.final   < Table0$pm_rd.final   & Table0$pm_rd.final   < Table0$upper_rd.final)
# 
# Table0$cp_pi0.initial <- (Table0$lower_pi0.initial < Table0$pm_pi0.initial & Table0$pm_pi0.initial < Table0$upper_pi0.initial)
# Table0$cp_pi0.final   <- (Table0$lower_pi0.final   < Table0$pm_pi0.final   & Table0$pm_pi0.final   < Table0$upper_pi0.final)
# 
# Table0$ss.initial <- Table0$y1.IP.initial + Table0$y0.IP.initial + Table0$y1.PC.initial + Table0$y0.PC.initial
# Table0$ss.final   <- Table0$y1.IP.final   + Table0$y0.IP.final   + Table0$y1.PC.final   + Table0$y0.PC.final

Table0$efficacy.initial <- (Table0$eff.prob.initial >= 0.975)
Table0$efficacy.final   <- (Table0$eff.prob.final   >= 0.975)

Table0$futility.initial <- (Table0$fut.prob.initial >= 0.975)
Table0$futility.final   <- (Table0$fut.prob.final   >= 0.975)

Table0$inconcl.initial  <- 1 - Table0$efficacy.initial - Table0$futility.initial
Table0$inconcl.final    <- 1 - Table0$efficacy.final   - Table0$futility.final


figure9 <- aggregate(. ~ eff.mix.prob + p.IP, data = Table0, mean)
########################################################################################################################################################
########################################################################################################################################################
########################################################################################################################################################

label_main =  ""
stretch     <- Table0$p.PC[1]*0.9
width.scale <- 7

if(output_png){png('../../../figure9.png',
    width = 450*width.scale, 
    height = 300*width.scale,
    pointsize=16,
    res=300)}

par(mar=c(5.1,4.1 + 0.5, 2.1, 2.1 + 0.5)) #c(bottom, left, top, right)
plot(figure9$p.IP, 
     figure9$efficacy.initial,
     type = 'l',
     ylim = c(0,1),
     lwd = 2,
     ylab = "Probability",
     xlab = "",
     #main = label_main,
     axes = FALSE)
box()

lines(figure9$p.IP,
      figure9$futility.initial,
      lwd = 2,
      lty = 'longdash')
lines(figure9$p.IP,
      figure9$inconcl.initial,
      lwd = 2,
      lty = 'dotted')

axis(2,
     las = 2,
     at = seq(0, 1, by = 0.1),
     labels = format(seq(0, 1, by = 0.1), nsmall = 1))
abline(h = seq(0, 1, by = 0.1),
       col = 'grey')

axis(1,
     las = 0,
     at = figure9$p.IP,
     labels = format(figure9$p.IP-figure9$p.PC, nsmall = 2))

mtext(bquote(theta),
      side=1,
      line=1,
      at=stretch)

row <- 2
for (j in seq(1,length(figure9$p.IP))){
  mtext(text = paste0(format(round(figure9$ss.initial[j],digits = 1),nsmall = 1),
                    " + ",
                    format(round(figure9$ss.final[j]-
                                 figure9$ss.initial[j],digits = 1),nsmall = 1),
                    " = ",
                    format(round(figure9$ss.final[j],digits = 1),nsmall = 1)),
        side = 1,line = row,at = figure9$p.IP[j])
}
mtext(text = "SS",side = 1,line = row,at = stretch)
row <- 3
for (j in seq(1,length(figure9$p.IP))){
  mtext(text = paste0("(I) ",
                    format(round(figure9$pm_rd.initial[j],digits = 3),nsmall = 3),
                    " (F) ",
                    format(round(figure9$pm_rd.final[j],digits = 3),nsmall = 3)),
        side = 1,line = row,at = figure9$p.IP[j])
}
mtext(text = "PM",
      side = 1,
      line = row,
      at = stretch)
row <- 4
for (j in seq(1,length(figure9$p.IP))){
  mtext(text = paste0("(I) ",
                    format(round(figure9$cp_rd.initial[j],digits = 3),nsmall = 3),
                    " (F) ",
                    format(round(figure9$cp_rd.final[j],digits = 3),nsmall = 3)),
        side = 1,line = row,at = figure9$p.IP[j])
}
mtext(text = "CP",
      side = 1,
      line = row,
      at = stretch)
row <- 5
for (j in seq(1,length(figure9$p.IP))){
  mtext(text = paste0("(I) ",
                      format(round(figure9$pm_pi0.initial[j],digits = 3),nsmall = 3),
                      " (F) ",
                      format(round(figure9$pm_pi0.final[j],digits = 3),nsmall = 3)),
        side = 1,line = row,at = figure9$p.IP[j])
}
mtext(text = "PM",side = 1,line = row,at = stretch)

row <- 6
for (j in seq(1,length(figure9$p.IP))){
  mtext(text = paste0("(I) ",
                      format(round(figure9$cp_pi0.initial[j],digits = 3),nsmall = 3),
                      " (F) ",
                      format(round(figure9$cp_pi0.final[j],digits = 3),nsmall = 3)),
        side = 1,line = row,at = figure9$p.IP[j])
}
mtext(text = "CP",side = 1,line = row,at = stretch)

points(figure9$p.IP[seq(1,length(figure9$p.IP))],
       figure9$eff.mon.initial[seq(1,length(figure9$p.IP))],
       pch = 20)

text(figure9$p.IP[seq(1,length(figure9$p.IP))],
     figure9$efficacy.initial[seq(1,length(figure9$p.IP))],
     labels = format(round(
       figure9$efficacy.initial[seq(1,length(figure9$p.IP))],
       digits = 3),nsmall = 3),
     pos = 3)

legend("topright",#text.width = 0.05,
       legend = c("Stop Early for Efficacy",
                "Stop Early for Futility",
                "Inconclusive Findings"),
       lty = c('solid','longdash','dotted'))

if(output_png){dev.off()}