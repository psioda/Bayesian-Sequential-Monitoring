for (table in c("Table1", "Table2", "Table3")){
  file_list <- list.files(paste0("../output/output062120/", table))
  dataset <- do.call("rbind", 
                     lapply(file_list, 
                            FUN = function(files){read.csv(paste0("../output/output062120/",table,"/", files), header = TRUE)}))
  write.csv(dataset, file=paste0("../output/output062120/", table, "_merged.csv"))
}

probs.p <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1) # posterior probability
inner.p <- array(NA, 
                 dim = c(reps, 2),
                 dimnames = list(seq_len(reps), c("initial.p", "final.p")))

outer.p.agree <- rep(NA, 3)
names(outer.p.agree) <- c("p.agree", "efficacy", "conditional")

outer                        <- apply(inner, MARGIN = 2, FUN=mean)
outer.p                      <- quantile(inner.p[,"final.p"][inner.p[,"initial.p"] > sig.eff & inner.p[,"final.p"] < sig.eff], probs = probs.p)
outer.p.agree["p.agree"]     <- sum((inner.p[,"initial.p"] > sig.eff) == (inner.p[,"final.p"] > sig.eff)) / reps
outer.p.agree["efficacy"]    <- sum((inner.p[,"initial.p"] > sig.eff) & (inner.p[,"final.p"] > sig.eff)) / reps
outer.p.agree["conditional"] <- sum((inner.p[,"initial.p"] > sig.eff) & (inner.p[,"final.p"] > sig.eff)) / sum(inner.p[,"initial.p"] > sig.eff)

Table1     <- data.frame(t(outer))
Table1$idx <- idx
write.csv(Table1, file = paste0("../output/Table1/", idx, "Table1.csv"))

Table2     <- data.frame(t(outer.p))
Table2$idx <- idx
write.csv(Table2, file = paste0("../output/Table2/", idx, "Table2.csv"))

Table3     <-data.frame(t(outer.p.agree))
Table3$idx <- idx
write.csv(Table3, file = paste0("../output/Table3/", idx, "Table3.csv"))