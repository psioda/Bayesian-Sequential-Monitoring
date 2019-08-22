library(ggplot2)
library(dplyr)
library(ggpubr)
n<-c(0,6,12,18)
y1<-c(0,1,2,3)

grid<-seq(0,1,length=50)

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

ggplot(pdat2, aes(dens, loc, fill = m, group = interaction(m, x))) + 
  geom_polygon() +
  scale_x_continuous(breaks = c(0,10,20,30)) +
  ylab('Probability of Response') +
  scale_y_continuous(breaks=seq(0,0.8,by=0.2)) +
  theme_minimal() +
  theme(axis.title.x = element_blank()) +
  geom_text(size = 3.5)
