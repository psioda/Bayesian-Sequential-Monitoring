model<-expand.grid(seq(0.39,0.51,by=0.01),.39,2,1,4,seq(0,1,by=0.25))
names(model)<-c("p.IP","p.PC","freq.mntr","enr.shape","out.mean","mix.prob")
write.csv(model,file="args_simulation.csv")