## The purpose of this program is to produce violin plots that show sample
## trajectories of sample data

rm(list = ls())

library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)
library(latex2exp)

## Design parameters ##
p.skpt<-0.20      # response rate for skeptic, enthusiast, futility
p.enth<-0.40
p.intr<-0.30
tail.skpt<-0.045  # tail probabilities for priors (low, high)
tail.enth<-0.05
sig.fut<-0.85     # significant trial result threshold
sig.eff<-0.95
cred.tail<-0.05   # credible interval is 1-cred.tail
max.ss<-76        # maximum sample size

par(mfrow=c(1,2))
#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################

# Step 1: Create grid for possible values of phi
phi.range<-seq(0,100,by=0.001)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to tail.skpt
quantiles.skpt<-qbeta(tail.skpt,(p.skpt)*phi.range,(1-(p.skpt))*phi.range,
                      lower.tail=FALSE)
# lower tail probability equal to tail.enth
quantiles.enth<-qbeta(tail.enth,(p.enth)*phi.range,(1-(p.enth))*phi.range,
                      lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi.range[which.min(abs(p.enth-quantiles.skpt))] # fixed 5/13/19
phi_H<-phi.range[which.min(abs(p.skpt-quantiles.enth))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha.skpt<-(p.skpt)*phi_L
beta.skpt<-(1-(p.skpt))*phi_L
alpha.enth<-(p.enth)*phi_H
beta.enth<-(1-(p.enth))*phi_H

#################################################################################################
## FUNCTIONS NEEDED #############################################################################
#################################################################################################

## posterior quantile
posterior.quantile<-function(a1,b1,a2,b2,y0,y1,q){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
       beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*qbeta(q,a1+y1,b1+y0)+
    (1-c)*qbeta(q,a2+y1,b2+y0)
  
  return(result)}

## posterior mean
posterior.mean<-function(a1,b1,a2,b2,y0,y1){
  
  c<-beta(a1+y1,b1+y0)/beta(a1,b1)/
    (beta(a1+y1,b1+y0)/beta(a1,b1)+
       beta(a2+y1,b2+y0)/beta(a2,b2))
  
  result<-c*(a1+y1)/(a1+b1+y1+y0)+
    (1-c)*(a2+y1)/(a2+b2+y1+y0)
  
  return(result)}

#################################################################################################
## GENERATE DATA ###############################################################################
#################################################################################################

n<-c(0,6,12,18)
y1<-c(0,0,1,1)
y0=n-y1

futility<-pbeta(p.intr,alpha.enth+y1,beta.enth+y0,lower.tail=TRUE)
efficacy<-pbeta(p.skpt,alpha.skpt+y1,beta.skpt+y0,lower.tail=TRUE)
mean<-posterior.mean(alpha.skpt,beta.skpt,alpha.enth,beta.enth,y0,y1)
lower<-posterior.quantile(alpha.skpt,beta.skpt,alpha.enth,beta.enth,y0,y1,0.025)
upper<-posterior.quantile(alpha.skpt,beta.skpt,alpha.enth,beta.enth,y0,y1,0.975)

#################################################################################################
## CREATE PLOTTING DATAFRAMES ###################################################################
#################################################################################################

grid<-seq(0,.8,length=200)

dens<-array(NA,dim=c(length(n),2,length(grid)),
                dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
loc<-array(NA,dim=c(length(n),2,length(grid)),
           dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
x<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
m<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))

for (i in 1:length(n)){
  dens[i,"skpt",]<-dbeta(grid,shape1=alpha.skpt+y1[i],shape2=beta.skpt+n[i]-y1[i])*(-1)+(i-1)*10
  dens[i,"enth",]<-dbeta(grid,shape1=alpha.enth+y1[i],shape2=beta.enth+n[i]-y1[i])+(i-1)*10
  loc[i,"skpt",]<-grid
  loc[i,"enth",]<-grid
  x[i,"skpt",]<-rep(as.character(i),length(grid))
  x[i,"enth",]<-rep(as.character(i),length(grid))
  m[i,"skpt",]<-rep("s",length(grid))
  m[i,"enth",]<-rep("e",length(grid))
}
pdat=data.frame(dens=as.vector(dens),
                loc=as.vector(loc),
                x=as.factor(as.vector(x)),
                m=as.factor(as.vector(m)))

pdat2<-pdat[order(pdat$x,pdat$m,-pdat$loc),]

#################################################################################################
## BEGIN GGPLOT #################################################################################
#################################################################################################
# initial plot
df<-data.frame(x=c(1,2),y=c(.2,.3))
p<- ggplot(pdat2, aes(dens, loc, fill = m, group = interaction(m, x))) + 
  geom_hline(yintercept=seq(0,.8,by=0.2), colour="black")+
  geom_polygon() + theme_classic() +
  theme(plot.margin = unit(c(0.5,0.5,3,0.5), "cm"))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y=element_blank())+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(0,0.8,by=0.2)) +
  geom_point(x=0,y=mean[1])+
  geom_point(x=10,y=mean[2])+
  geom_point(x=20,y=mean[3])+
  geom_point(x=30,y=mean[4])+
  geom_point(x=0,y=lower[1])+
  geom_point(x=10,y=lower[2])+
  geom_point(x=20,y=lower[3])+
  geom_point(x=30,y=lower[4])+
  geom_point(x=0,y=upper[1])+
  geom_point(x=10,y=upper[2])+
  geom_point(x=20,y=upper[3])+
  geom_point(x=30,y=upper[4])+
  geom_segment(aes(x = 0, y = mean[1], xend = 10, yend = mean[2]))+
  geom_segment(aes(x = 10, y = mean[2], xend = 20, yend = mean[3]))+
  geom_segment(aes(x = 20, y = mean[3], xend = 30, yend = mean[4]))+
  geom_segment(aes(x = 0, y = lower[1], xend = 0, yend = upper[1]))+
  geom_segment(aes(x = 10, y = lower[2], xend = 10, yend = upper[2]))+
  geom_segment(aes(x = 20, y = lower[3], xend = 20, yend = upper[3]))+
  geom_segment(aes(x = 30, y = lower[4], xend = 30, yend = upper[4]))
p


for (i in 1:length(n)){
p <- p + annotation_custom(
  grob = textGrob(label = 
  paste0(format(round(futility[i]*100,digits=1),nsmall=1),"%"),
  hjust = 0, gp = gpar(cex = 1)),
  ymin = -.15,ymax = -.15,xmin = (i-1)*10,xmax = (i-1)*10)+
  annotation_custom(
  grob = textGrob(label =
 paste0(format(round((1-efficacy[i])*100,digits=1),nsmall=1),"%"),
  hjust = 0, gp = gpar(cex = 1)),
  ymin = -.12,ymax = -.12,xmin = (i-1)*10,xmax = (i-1)*10)+
 annotation_custom(
  grob = textGrob(label = y1[i], hjust = 0, gp = gpar(cex = 1)),
  ymin = -.09,ymax = -.09,xmin = (i-1)*10,xmax = (i-1)*10)+
 annotation_custom(
  grob = textGrob(label = n[i], hjust = 0, gp = gpar(cex = 1)),
  ymin = -.06,ymax = -.06,xmin = (i-1)*10,xmax = (i-1)*10)
}
p
p <- p + annotation_custom(
  grob = textGrob(label =TeX('$P(\\theta\\leq 0.30|\\mathbf{D},\\pi_E)$'), hjust = 0, gp = gpar(cex = 1)),
  ymin = -.15,ymax = -.15,xmin = -5,xmax = -5)+
  annotation_custom(
    grob = textGrob(label = TeX('$P(\\theta>0.20|\\mathbf{D},\\pi_S)$'), hjust = 0, gp = gpar(cex = 1)),
    ymin = -.12,ymax = -.12,xmin = -5,xmax = -5)+
  annotation_custom(
    grob = textGrob(label = "# Successes", hjust = 0, gp = gpar(cex = 1)),
    ymin = -.09,ymax = -.09,xmin = -5,xmax = -5)+
  annotation_custom(
    grob = textGrob(label = "Sample Size", hjust = 0, gp = gpar(cex = 1)),
    ymin = -.06,ymax = -.06,xmin = -5,xmax = -5)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)

#################################################################################################
## DO NOT USE BELOW #############################################################################
#################################################################################################


# p<- ggplot(pdat2, aes(dens, loc, fill = m, group = interaction(m, x))) + 
#   geom_polygon() + theme_classic() +
#   theme(panel.border = element_rect(colour = "black", fill=NA, size=0.5)) +
#   geom_point(aes(x=10,y=0))
# scale_x_continuous(breaks = c(-10,0,10,20,30)) +
#   ylab('Probability of Response') +
#   #ylim(-.2,0.8) +
#   scale_y_continuous(breaks=seq(0,0.6,by=0.2)) +
#   theme_minimal() +
#   theme(axis.title.x = element_blank())

