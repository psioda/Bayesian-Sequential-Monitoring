model<-expand.grid(seq(0.39,0.51,by=0.005),.39,10,1,4)
names(model)<-c("p.IP","p.PC","freq.mntr","enr.shape","out.mean")
write.csv(model,file="args_simulation.csv")