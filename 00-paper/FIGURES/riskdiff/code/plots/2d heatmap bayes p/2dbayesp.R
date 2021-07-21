rm(list = ls())
setwd("/Users/evankwiatkowski/Documents/GitHub/Bayesian-Sequential-Monitoring/00-paper/FIGURES/riskdiff/output")

dat <- read.csv("Table0_merged_2021-07-21-bayes-p-grid.csv")

head(dat)
summary(dat$eff.mix.prob.initial)
summary(dat$eff.mix.prob.final)
summary(dat$box.enth.initial)
summary(dat$box.skpt.initial)
