rm(list = ls())
setwd("D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES")

### Model arguments  ####################################
p.skpt<-0.40     # response rate for skeptic, enthusiast, futility
p.enth<-0.67
p.intr<-(p.skpt+p.enth)/2
tail.skpt<-0.05  # tail probabilities for priors
tail.enth<-0.05
cred.tail<-0.05
sig.fut<-0.95
sig.eff<-0.95
lambda0.skpt<-2
lambda0.enth<-2
epsilon<-0 # used to stop numerical error from integration

### Parameterize generalized normal priors  #############
source("scripts/code_headers_gen_nrml.R")

### Prior plots  ########################################
#png(filename="normal_priors.png")
source("plots/plots_priors.R")
#dev.off()

### Functions  ##########################################
source("functions/functions.R")

### Violin plots ########################################
n<-c(0,6,12,18) # Efficacy example
y1<-c(0,4,8,12) 
miss<-c(0,0,0,0)
y0=n-y1
p.range<-1 # temp
j<-1       # temp
stretch<-1.5

source("plots/plots_violin.R")
plot1 <- ggplot_gtable(ggplot_build(p))
plot1$layout$clip[plot1$layout$name == "panel"] <- "off"
png(filename="violin_efficacy_955.png")
grid.arrange(plot1)
dev.off()

n<-   c(0,6,12,18,24,30) # Futility example
y1<-  c(0,2,4, 6, 8,10)
miss<-c(0,0,0, 0,0,0)
y0=n-y1
stretch<-2.5 # "Sample Size" labels xlim
source("plots/plots_violin.R")
plot2 <- ggplot_gtable(ggplot_build(p))
plot2$layout$clip[plot2$layout$name == "panel"] <- "off"

# plot side-by-side
grid.arrange(plot1,plot2,ncol=2)

### Simulations ########################################
max.ss<-80
reps<-2000
p.range<-seq(p.skpt,p.enth,by=(p.enth-p.skpt)/5)
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
#pdf('spike_2019-09-18.pdf',height=6,width=15)
source("plots/plots_seq_design_prop.R")
#dev.off()

