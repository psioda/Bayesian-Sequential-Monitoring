library(ggplot2)
library(dplyr)
library(ggpubr)
library(grid)
library(latex2exp)
library(gridExtra)

#################################################################################################
## GENERATE DATA ###############################################################################
#################################################################################################

efficacy<-NA
futility<-NA
mean<-NA
lower<-NA
upper<-NA
for (i in 1:length(n)){
efficacy[i]<-eff_fut(i)[1]
futility[i]<-eff_fut(i)[2]
mean[i]<-pm_cp(i)[1]
lower[i]<-pm_cp(i)[4]
upper[i]<-pm_cp(i)[5]
}

#################################################################################################
## CREATE PLOTTING DATAFRAMES ###################################################################
#################################################################################################
grid<-seq(0+1E-4,1-1E-4,length=1000)
spacing<-(seq(1:length(n))-1)*16
dens<-array(NA,dim=c(length(n),2,length(grid)),
                dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
loc<-array(NA,dim=c(length(n),2,length(grid)),
           dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
x<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))
m<-array(NA,dim=c(length(n),2,length(grid)),
         dimnames=list(seq_len(length(n)),c("skpt","enth"),seq_len(length(grid))))

for (i in 1:length(n)){
  posterior.skpt<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))
  }
  skpt.nc<-integrate(posterior.skpt,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.skpt<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.skpt)/sigma0.skpt)^lambda0.skpt))/skpt.nc
  }
  posterior.enth<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))
  }
  enth.nc<-integrate(posterior.enth,lower=0+epsilon,upper=1-epsilon)[[1]]
  posterior.nc.enth<-function(x){
    exp(y1[i]*log(x)+y0[i]*log(1-x)-((abs(x-mu0.enth)/sigma0.enth)^lambda0.enth))/enth.nc
  }

  dens[i,"skpt",]<-posterior.nc.skpt(grid)*(-1)+spacing[i]
  dens[i,"enth",]<-posterior.nc.enth(grid)+spacing[i]
  loc[i,"skpt",]<-grid
  loc[i,"enth",]<-grid
  x[i,"skpt",]<-rep(as.character(i),length(grid))
  x[i,"enth",]<-rep(as.character(i),length(grid))
  m[i,"skpt",]<-rep("s",length(grid))
  m[i,"enth",]<-rep("e",length(grid))
}

dens[1,"skpt",1]<-0 # 12/8/19
dens[1,"enth",length(grid)]<-0 # 12/8/19

pdat=data.frame(dens=as.vector(dens),
                loc=as.vector(loc),
                x=as.factor(as.vector(x)),
                m=as.factor(as.vector(m)))

pdat2<-pdat[order(pdat$x,pdat$m,-pdat$loc),]

#################################################################################################
## BEGIN GGPLOT #################################################################################
#################################################################################################
par(mfrow=c(1,1))
# initial plot
df<-data.frame(x=c(1,2),y=c(.2,.3))
p<- ggplot(pdat2, aes(dens, loc, fill = m, group = interaction(m, x))) + 
  geom_hline(yintercept=seq(0,1,by=0.2), colour="black")+
  geom_polygon() + theme_classic() +
  scale_fill_manual(values = c("darkgrey", "lightgrey"))+
  # c(bottom, left, top, right)
  theme(plot.margin = unit(c(0,0,2,4), "cm"))+ # TRouBLe 12/8/19 change Bottom to 2cm to include table underneath
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  theme(axis.title.x = element_blank())+
  theme(axis.text.x=element_blank())+
  theme(axis.title.y = element_blank())+
  theme(axis.text.y=element_blank())+
  #ylab("Probability of Response")+
  #theme(axis.title.y=element_blank())+
  #theme(axis.text.y="Probability of Success")+
  theme(legend.position="none")+
  scale_y_continuous(breaks=seq(0,1,by=0.2))

for (i in 1:length(n)){
p<-p+ geom_point(x=spacing[i],y=mean[i])+
      geom_point(x=spacing[i],y=lower[i],shape="-")+
      geom_point(x=spacing[i],y=upper[i],shape="-")+
      geom_segment(x=spacing[i],y=lower[i],xend=spacing[i],yend=upper[i])
}
  
for (i in 2:length(n)){
p<-p+geom_segment(x = spacing[i-1], 
                      y = mean[i-1],
                      xend = spacing[i], 
                      yend = mean[i])
}
p <- p +  geom_segment(x = spacing[2]/2, y = 0,xend = spacing[2]/2, yend = 1,
                 linetype=2)+
          geom_segment(x = spacing[length(n)]-spacing[2]/2, y = 0, 
                     xend = spacing[length(n)]-spacing[2]/2, yend = 1,
                 linetype=2)

for (i in 1:length(n)){
p <- p + annotation_custom(
  grob = textGrob(label =
  paste0(format(round(futility[i],digits=3),nsmall=3)),
  hjust = 0, gp = gpar(cex = 1)),
  ymin = -.18,ymax = -.18,xmin = spacing[i],xmax = spacing[i])+
  annotation_custom(
  grob = textGrob(label =
 paste0(format(round((1-efficacy[i]),digits=3),nsmall=3)),
  hjust = 0, gp = gpar(cex = 1)),
  ymin = -.15,ymax = -.15,xmin = spacing[i],xmax = spacing[i])+
  annotation_custom(
    grob = textGrob(label = miss[i], hjust = 0, gp = gpar(cex = 1)),
    ymin = -.12,ymax = -.12,xmin = spacing[i],xmax = spacing[i])+
 annotation_custom(
  grob = textGrob(label = y1[i], hjust = 0, gp = gpar(cex = 1)),
  ymin = -.09,ymax = -.09,xmin = spacing[i],xmax = spacing[i])+
 annotation_custom(
  grob = textGrob(label = n[i], hjust = 0, gp = gpar(cex = 1)),
  ymin = -.06,ymax = -.06,xmin = spacing[i],xmax = spacing[i])
}


p <- p +
  annotation_custom(
     grob = textGrob(label =TeX('$P(\\theta\\leq \\theta_0|\\mathbf{D},\\pi_E)$'), hjust = 0, gp = gpar(cex = 1)),
     ymin = -.18,ymax = -.18,xmin = -spacing[2]*stretch,xmax = -spacing[2]*stretch)+
   annotation_custom(
     grob = textGrob(label = TeX('$P(\\theta > (\\theta_1+\\theta_0)/2|\\mathbf{D},\\pi_S)$'), hjust = 0, gp = gpar(cex = 1)),
     ymin = -.15,ymax = -.15,xmin = -spacing[2]*stretch,xmax = -spacing[2]*stretch)+
  annotation_custom(
    grob = textGrob(label = "# Missing", hjust = 0, gp = gpar(cex = 1)),
    ymin = -.12,ymax = -.12,xmin = -spacing[2]*stretch,xmax = -spacing[2]*stretch)+
  annotation_custom(
    grob = textGrob(label = "# Successes", hjust = 0, gp = gpar(cex = 1)),
    ymin = -.09,ymax = -.09,xmin = -spacing[2]*stretch,xmax = -spacing[2]*stretch)+
  annotation_custom(
    grob = textGrob(label = "Sample Size", hjust = 0, gp = gpar(cex = 1)),
    ymin = -.06,ymax = -.06,xmin = -spacing[2]*stretch,xmax = -spacing[2]*stretch)


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

