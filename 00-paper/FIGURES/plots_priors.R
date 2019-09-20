# Step 5: Plot Skeptical and enthuastic priors separately
## 7/19/19 match slide 14 from FDA presentation, shade in gray scale, set 350 DPI
par(ask=FALSE)
par(oma = c(1.25, 1.25, 0, 0)) # make room (i.e. the 4's) for the overall x and y axis titles
par(mar=c(5.1-1,4.1,4.1-4,2.1),mai = c(0.3, 0.3, 0.1, 0.1)) #bottom, left, top, and right.
par(mfrow = c(1, 3)) 
x<-seq(0,1,by=0.01)
# Make 2 boxplots
#low (skeptical)
plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",
     xlab="",ylab="Density Value",
     main="",
     #xaxt="n",
     ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt))))

polygon(c(x[x<=0.4],0.4),c(dbeta(x,alpha.skpt,beta.skpt)[x<=0.4],0),col="lightgrey")
polygon(c(x[x>=0.4],0.4),c(dbeta(x,alpha.skpt,beta.skpt)[x>=0.4],0),col="black")
#text(0.5, 4,expression(P(theta>0.40)==0.045),pos=4)
#text(0.5,3.8,"Approximately N=14.0 Subjects",pos=4)
#text(0.5,3.6,expression(E(theta)==0.20),pos=4)
#axis(1,at=p.skpt,labels=expression(theta[0]))
#axis(1,at=seq(0,1,by=0.2))
#high (enthuastic)
plot(x,dbeta(x,alpha.enth,beta.enth),type="l",
     xlab="Response Probability",
     ylab="",
     main="",
     #xaxt="n",
     yaxt="n",
     ylim=c(0,max(dbeta(x,alpha.enth,beta.enth))))
axis(2,at=c(0,1,2,3),labels=c(0,1,2,3))
polygon(c(x[x<=0.2],0.2),c(dbeta(x,alpha.enth,beta.enth)[x<=0.2],0),col="black")
polygon(c(x[x>=0.2],0.2),c(dbeta(x,alpha.enth,beta.enth)[x>=0.2],0),col="lightgrey")
#text(0.5, 3,expression(P(theta<0.20)==0.05),pos=4)
#text(0.5,2.85,"Approximately N=14.0 Subjects",pos=4)
#text(0.5,2.7,expression(E(theta)==0.40),pos=4)



# Step 6: Plot inference priors (mixtures)
## 7/19/19 make 50:50 solid black line, 25:75 dark grey with different dash patterns,
# make skeptical/enthuastic light grey with smaller thinkness
# no legend, put omega on actual lines


plot(x,1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth),
     type="l",xlab="",ylab="",yaxt="n",
     main="",
     ylim=c(0,max(1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth))))
y<-1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth)
axis(2, at = c(0,1,2), labels=c("0","1","2"))
polygon(c(x,1),c(y,0),col="lightgrey")
polygon(c(0.2,x[x>=0.2 & x<=0.4],0.4),c(0,y[x>=0.2 & x<=0.4],0),col="darkgrey")

# print the overall labels
mtext('Response Probability', side = 1, outer = TRUE, line = 0)
mtext('Density Value', side = 2, outer = TRUE, line = 0)
# par(mfrow = c(1, 1)) 
# plot(x,dbeta(x,alpha.skpt,beta.skpt),type="l",col="black",xlab="Response Probability",ylab="Density Value",main="Inference Priors",
#      xaxt="n",
#      ylim=c(0,max(dbeta(x,alpha.skpt,beta.skpt),dbeta(x,alpha.enth,beta.enth))))
# axis(1,at=p.skpt,labels=expression(theta[0]))
# lines(x,dbeta(x,alpha.enth,beta.enth),type="l",col="lightgrey")
# axis(1,at=p.enth,labels=expression(theta[A]))
# axis(1,at=p.skpt,labels=expression(theta[0]))
# points(x,1/4*dbeta(x,alpha.skpt,beta.skpt)+3/4*dbeta(x,alpha.enth,beta.enth),type='l',lty=3)
# points(x,1/2*dbeta(x,alpha.skpt,beta.skpt)+1/2*dbeta(x,alpha.enth,beta.enth),type='l',lty=2)
# points(x,3/4*dbeta(x,alpha.skpt,beta.skpt)+1/4*dbeta(x,alpha.enth,beta.enth),type='l',lty=1)
# legend(x=0.55,y=4.5,
#        legend=c("Skeptical","75:25","50:50","25:75","Enthuastic"),
#        col=c("black","black","black","black","lightgrey"),lty=c(1,1,2,3,1),
#        cex=0.8,text.font=3)