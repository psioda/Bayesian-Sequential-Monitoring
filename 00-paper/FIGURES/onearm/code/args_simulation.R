model1 <- expand.grid(1,
                    seq(p.skpt, p.enth, length = 25),
                    c(0, 1),
                    c(0, 1),
                    10,
                    1,
                    4)
# model2 <- expand.grid(2,
#                     p.skpt,
#                     1,
#                     0,
#                     c(1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60),
#                     c(1, 0.25),
#                     4)
#  
# model <- rbind(model1, model2)

model <- model1
names(model) <- c("model", 
                "p.range", 
                "skpt_spike", 
                "enth_flat", 
                "freq.mntr", 
                "enr.shape", 
                "out.mean")
write.csv(model, file="args_simulation.csv")

#freq.mnt=frequency of monitoring
#enr.shape=shape gamma dist enrollment
#out.mean=mean normal dist outcome