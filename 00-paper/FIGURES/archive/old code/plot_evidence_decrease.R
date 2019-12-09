outer.p[i,j,]<-quantile(inner.p[,"final.p"][inner.p[,"initial.p"]>sig.eff],probs=probs.p)
outer.p.agree[i,j,"p.agree"]<-sum((inner.p[,"initial.p"]>sig.eff)==(inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"efficacy"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/reps
outer.p.agree[i,j,"conditional"]<-sum((inner.p[,"initial.p"]>sig.eff) & (inner.p[,"final.p"]>sig.eff))/
  
  
outer.p
outer.p.agree

outer.p
