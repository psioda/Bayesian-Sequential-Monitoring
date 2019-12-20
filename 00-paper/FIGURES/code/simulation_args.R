model<-expand.grid(seq(p.skpt,p.enth,by=0.01),0,0,2,1,4)
names(model)<-c("p.range","prior_indicator_skpt","prior_indicator_enth",
                "freq.mntr","enr.shape","out.mean")
write.csv(model,file="simulation_args.csv")

#freq.mntr<-rep(c(1,2,4,8,16,76),2)     # frequency of monitoring
#enr.shape<-c(rep(1,6),rep(0.25,6))      # shape gamma dist enrollment
#out.mean<-c(rep(4,12))       # mean normal dist outcome