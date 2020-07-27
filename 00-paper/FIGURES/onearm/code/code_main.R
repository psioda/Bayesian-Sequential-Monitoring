matrix.names<-c("eff.mon.initial", "eff.mon.final", 
                "eff.inf.initial", "eff.inf.final", 
                "fut.mon.initial", "fut.mon.final", 
                "fut.inf.initial", "fut.inf.final", 
                "phat.initial", "phat.final", 
                "ss.initial", "ss.final", 
                "post.mean.initial", "post.mean.final", 
                "cov.initial", "cov.final")
outer<-array(NA, 
             dim=c(length(freq.mntr), length(p.range), length(matrix.names)), 
             dimnames=list(seq_len(length(freq.mntr)), p.range, matrix.names))
inner<-array(NA, 
             dim=c(reps, length(matrix.names)), 
             dimnames=list(seq_len(reps), matrix.names))

# probs.p <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
probs.p <- seq(0,  1,  by = 0.01)

outer.p<-array(NA, 
               dim=c(length(freq.mntr), length(p.range), length(probs.p)), 
               dimnames=list(seq_len(length(freq.mntr)), p.range, probs.p))
inner.p<-array(NA, 
               dim=c(reps, 2), 
               dimnames=list(seq_len(reps), c("initial.p", "final.p")))
outer.p.agree<-array(NA, 
                     dim=c(length(freq.mntr), length(p.range), length(c("p.agree", "efficacy", "conditional"))), 
                     dimnames=list(seq_len(length(freq.mntr)), p.range, c("p.agree", "efficacy", "conditional")))

for (i in 1:length(freq.mntr)){
  
  for (j in 1:length(p.range)){
    
    for (k in 1:reps){
      
      if (k%%100==0){print(paste0("Model ", idx, ",  Response p=", p.range[j], ",  Simulation ", k))}
      
      enr.times     <- cumsum(rgamma(n=max.ss, shape=enr.shape[i], scale=0.5))
      outcome.times <- sort(enr.times+rnorm(n=max.ss, mean=out.mean[i], sd=0.25))
      
      responses <- rbinom(n=max.ss, size=1, prob=p.range[j])
      y1        <- cumsum(responses)
      y0        <- seq(1:length(responses))-y1
      
      futility <- rep(NA, length=max.ss)
      efficacy <- rep(NA, length=max.ss)
      
      for (h in 1:max.ss){
        
        results     <- eff_fut(h)
        efficacy[h] <- results[1]
        futility[h] <- results[2]
      }
      
      n.initial<-min(which(((futility>sig.fut) | ((1-efficacy)>sig.eff)) & (seq(1:max.ss)%%freq.mntr[i]==0)), 
                     max.ss, 
                     na.rm=TRUE)
      
      cutoff.time <- outcome.times[n.initial]
      n.final     <- length(responses[enr.times<=cutoff.time])
      
      time <- c("initial", "final")
      n    <- c(n.initial, n.final)
      for (l in 1:2){
        inner[k, paste0("fut.mon.", time[l])]   <- (futility[n[l]]>sig.fut)
        inner[k, paste0("eff.mon.", time[l])]   <- (1-efficacy[n[l]]>sig.eff)
        inner[k, paste0("ss.", time[l])]        <- n[l]
        inner[k, paste0("phat.", time[l])]      <- y1[n[l]]/n[l]
        
        results                                 <- pm_cp(n[l])
        inner[k, paste0("post.mean.", time[l])] <- results[1]
        inner[k, paste0("eff.inf.", time[l])]   <- results[2]
        inner[k, paste0("cov.", time[l])]       <- results[3]
        inner.p[k, paste0(time[l], ".p")]       <- (1-efficacy[n[l]])
      }
    }
    
    outer[i, j, ]   <- apply(inner, MARGIN=2, FUN=mean)
    outer.p[i, j, ] <- quantile(inner.p[, "final.p"][inner.p[, "initial.p"]>sig.eff & inner.p[, "final.p"]<sig.eff], 
                                probs=probs.p)
    
    outer.p.agree[i, j, "p.agree"]     <- sum((inner.p[, "initial.p"]>sig.eff)==(inner.p[, "final.p"]>sig.eff))/reps
    outer.p.agree[i, j, "efficacy"]    <- sum((inner.p[, "initial.p"]>sig.eff) & (inner.p[, "final.p"]>sig.eff))/reps
    outer.p.agree[i, j, "conditional"] <- sum((inner.p[, "initial.p"]>sig.eff) & (inner.p[, "final.p"]>sig.eff))/sum(inner.p[, "initial.p"]>sig.eff)
  }
}