model<-expand.grid(seq(0,0.5,by=0.1),seq(0,0.5,by=0.1),10,1,4)
names(model)<-c("p.IP","p.PC","freq.mntr","enr.shape","out.mean")
model<-model[(
  (round(model$p.IP-model$p.PC+.1,1)==0) |
    (round(model$p.IP-model$p.PC,1)==0)  |
    (round(model$p.IP-model$p.PC-.1,1)==0)  |
    (round(model$p.IP-model$p.PC-.2,1)==0) ),]

write.csv(model,file="model.csv",row.names=FALSE)