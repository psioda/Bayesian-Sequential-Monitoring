model1<-expand.grid(1,seq(p.skpt,p.enth,by=0.01),1,0,2,1,4)
model2<-expand.grid(2,p.skpt,1,0,seq(1,max.ss,by=1),1,4)
model3<-expand.grid(3,seq(p.skpt,p.enth,length=5),1,0,2,c(1,0.25),4)
#model3<-expand.grid(3,seq(p.skpt,p.enth,length=3),1,0,c(1,2,3,5,10,20,61,112),c(1,0.25),4)
model4<-expand.grid(4,seq(p.skpt,p.enth,length=5),c(0,1),c(0,1),2,1,4)
model<-rbind(model1,model2,model3,model4)
names(model)<-c("model","p.range","skpt_spike","enth_flat",
                "freq.mntr","enr.shape","out.mean")
write.csv(model,file="args_simulation.csv")

#freq.mnt=frequency of monitoring
#enr.shape=shape gamma dist enrollment
#out.mean=mean normal dist outcome