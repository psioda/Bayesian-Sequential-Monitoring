rm(list = ls())

for (idx in 1:4){
if (.Platform$OS.type == "windows") {
  root<-"D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES"
  setwd(root)
}

### Batch arguments   ####################################
prior_indicator_skpt<-c(0,0,1,1)
prior_indicator_enth<-c(0,1,0,1)

### Model arguments  ####################################
p.skpt<-0.40     # response rate for skeptic, enthusiast, futility
p.enth<-0.67
p.intr<-(p.skpt+p.enth)/2
tail.skpt<-0.025  # tail probabilities for priors
tail.enth<-0.025
cred.tail<-0.05
sig.fut<-0.975
sig.eff<-0.975
lambda0.skpt<-2
lambda0.enth<-2
epsilon<-0 # used to stop numerical error from integration

par(mfrow=c(2,2))
### Parameterize generalized normal priors  #############
if (prior_indicator_skpt[idx]==0){
  source("scripts/code_skpt_prior_default.R")
  source("plots/plots_priors.R")
}
if (prior_indicator_skpt[idx]==1){
  scale<-1.15
  source("scripts/code_skpt_prior_custom.R")
  source("plots/plots_priors.R")
}
if (prior_indicator_enth[idx]==0){
  source("scripts/code_enth_prior_default.R")
  source("plots/plots_priors.R")
}
if (prior_indicator_enth[idx]==1){
  scale<-0.85
  source("scripts/code_enth_prior_custom.R")
  source("plots/plots_priors.R")
}
# print the overall labels
mtext('Response Probability', side = 1, outer = TRUE, line = .3)
mtext('Density Value', side = 2, outer = TRUE, line = .3)

### Prior plots  ########################################
png(filename=paste0("normal_priors",prior_indicator_skpt[idx],prior_indicator_enth[idx],".png"))
source("plots/plots_priors.R")
dev.off()

### Functions  ##########################################
source("functions/functions.R")

# ### Violin plots - Efficacy #############################
# n<-c(0,6,12,18) # 
# y1<-c(0,4,8,12) 
# miss<-c(0,0,0,0)
# y0=n-y1
# p.range<-1 # temp
# j<-1       # temp
# stretch<-1.5
# 
# source("plots/plots_violin.R")
# plot1 <- ggplot_gtable(ggplot_build(p))
# plot1$layout$clip[plot1$layout$name == "panel"] <- "off"
# #png(filename="violin_efficacy_955.png")
# grid.arrange(plot1)
# #dev.off()
# 
# ### Violin plots - Efficacy #############################
# n<-   c(0,6,12,18,24,30) 
# y1<-  c(0,2,4, 6, 8,10)
# miss<-c(0,0,0, 0,0,0)
# y0=n-y1
# stretch<-2.5 # "Sample Size" labels xlim
# source("plots/plots_violin.R")
# plot2 <- ggplot_gtable(ggplot_build(p))
# plot2$layout$clip[plot2$layout$name == "panel"] <- "off"
# 
# # plot side-by-side
# grid.arrange(plot1,plot2,ncol=2)

### Simulations ########################################
max.ss<-112
reps<-1000
p.range<-seq(p.skpt,p.intr,by=(p.intr-p.skpt)/2)
#p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05)
freq.mntr<-2
#freq.mntr<-rep(c(1,2,4,8,16,76),2)     # frequency of monitoring
enr.shape<-1
#enr.shape<-c(rep(1,6),rep(0.25,6))      # shape gamma dist enrollment
out.mean<-4
#out.mean<-c(rep(4,12))       # mean normal dist outcome
source("scripts/code_main.R")

### Sequential design properties ########################
label_main=""
stretch<-p.skpt*0.8
pdf(file=paste0("seq_dsn",prior_indicator_skpt[idx],prior_indicator_enth[idx],".pdf"),height=6,width=15)
source("plots/plots_seq_design_prop.R")
dev.off()
}
