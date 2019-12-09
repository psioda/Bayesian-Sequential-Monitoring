rm(list = ls())

setwd("D:/Users/ekwiatko/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/")

source("code_functions.R")

p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors
tail.enth<-0.05
source("code_headers_beta.R")

# Efficacy example
n<-c(0,8,16,19)
y1<-c(0,4,9,11)
miss<-c(0,11,3,0)
y0=n-y1
stretch<-1.5
source("plots_violin.R")
plot1 <- ggplot_gtable(ggplot_build(p))
plot1$layout$clip[plot1$layout$name == "panel"] <- "off"
grid.arrange(plot1,ncol=1)


# Futility example
n<-   c(0,8,16,24,32,40)
y1<-  c(0,1,2, 5, 5,8)
miss<-c(0,6,4, 13,8,0)
y0=n-y1
stretch<-2.5 # "Sample Size" labels xlim
source("plots_violin.R")
plot2 <- ggplot_gtable(ggplot_build(p))
plot2$layout$clip[plot2$layout$name == "panel"] <- "off"
grid.arrange(plot2,ncol=1)


# plot side-by-side
grid.arrange(plot1,plot2,ncol=2)
