rm(list = ls())

setwd("D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/")

source("code_functions.R")

p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors (low, high)
tail.enth<-0.05
source("code_headers.R")

sig.fut<-0.85
sig.eff<-0.95
cred.tail<-0.05
max.ss<-76
reps<-50000
p.range<-0.20
#p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05)
freq.mntr<-c(1,2,4,8,16,76)
enr.shape<-1
out.mean<-8
source("code_main_nested.R")

spike<-0 # spike/slab version or regular version
sig.fut<-0.85
sig.eff<-0.95
cred.tail<-0.05
max.ss<-76
reps<-10000
p.range=0.2
#p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05) # range of response proportion
freq.mntr<-rep(c(1,2,4,8,16,76),2)     # frequency of monitoring
enr.shape<-c(rep(1,6),rep(0.25,6))      # shape gamma dist enrollment
out.mean<-c(rep(4,12))       # mean normal dist outcome
source("code_main.R")
label_main=""
stretch<-0.125
source("plots_seq_design_prop.R")

# p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
# p.enth<-0.40
# p.intr<-0.30
# tail.skpt<-0.045  # tail probabilities for priors (low, high)
# tail.enth<-0.05
# source("code_headers.R")
# 
# spike<-0 # spike/slab version or regular version
# sig.fut<-0.85
# sig.eff<-0.95
# cred.tail<-0.05
# max.ss<-76
# reps<-60000
# #p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05) # range of response proportion
# p.range<-0.2
# #freq.mntr<-2      # frequency of monitoring
# freq.mntr<-c(1,2,4,8,16,76)
# #freq.mntr<-rep(c(1,2,4,8,16,80),4)
# #enr.shape<-1      # shape gamma dist enrollment
# #enr.shape<-c(rep(1,12),rep(0.25,12))
# enr.shape<-rep(1,6)
# #out.mean<-4       # mean normal dist outcome
# out.mean<-rep(4,6)
# #out.mean<-rep(c(rep(4,6),rep(8,6)),2)
# source("code_main.R")
# source("plots_t1e.R")

# spike<-0 # spike/slab version or regular version
# sig.fut<-0.85
# sig.eff<-0.95
# cred.tail<-0.05  
# max.ss<-76  
# reps<-10000      
# p.range<-seq(p.skpt-0.05,p.enth+0.05,by=0.05) # range of response proportion
# #p.range<-0.15                             
# freq.mntr<-2      # frequency of monitoring
# #freq.mntr<-c(1,2,4,8,16,76)
# #freq.mntr<-rep(c(1,2,4,8,16,80),4) 
# enr.shape<-1      # shape gamma dist enrollment
# #enr.shape<-c(rep(1,12),rep(0.25,12))
# #enr.shape<-rep(1,6)
# out.mean<-4       # mean normal dist outcome
# #out.mean<-rep(4,6)
# #out.mean<-rep(c(rep(4,6),rep(8,6)),2)
# source("code_main.R")

## Sequential design properties
# label_main=""
# stretch<-0.125
# source("plots_seq_design_prop.R")

## Futility example
#n<-   c(0,8,16,24,32,40)
#y1<-  c(0,1,2, 5, 5,8)
#miss<-c(0,6,4, 13,8,0)
#stretch<-2.5 # "Sample Size" labels xlim
## Efficacy example
# n<-c(0,8,16,19)
# y1<-c(0,4,9,11)
# miss<-c(0,11,3,0)
# y0=n-y1
# stretch<-1.5
# source("plots_violin.R")
# gt <- ggplot_gtable(ggplot_build(p))
# gt$layout$clip[gt$layout$name == "panel"] <- "off"
# grid.draw(gt)
