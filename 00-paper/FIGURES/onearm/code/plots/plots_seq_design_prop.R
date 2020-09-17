plot(figure3$p.range,figure3$eff.mon.initial,
     type='l',ylim=c(0,1),lwd=2,
     ylab="Probability",xlab="",
     main=label_main,
     axes=FALSE)
box()

lines(figure3$p.range,figure3$fut.mon.initial,lwd=2,lty='longdash')
lines(figure3$p.range,1-figure3$fut.mon.initial-figure3$eff.mon.initial,lwd=2,lty='dotted')


axis(2,las=2,at=seq(0,1,by=0.1),labels=format(seq(0,1,by=0.1),nsmall=1))
abline(h=seq(0,1,by=0.1),col='grey')

axis(1,las=0,at=figure3_table$p.range[seq(1,length(figure3_table$p.range))],
     labels=format(figure3_table$p.range[seq(1,length(figure3_table$p.range))],nsmall=2))

mtext(bquote(theta),side=1,line=1,at=stretch)
row<-2
for (j in seq(1,length(figure3_table$p.range))){
  mtext(text=paste0(format(round(figure3_table$ss.initial[j],digits=1),nsmall=1),
                    " + ",
                    format(round(figure3_table$ss.final[j]-
                                 figure3_table$ss.initial[j],digits=1),nsmall=1),
                    " = ",
                    format(round(figure3_table$ss.final[j],digits=1),nsmall=1)),
        side=1,line=row,at=figure3_table$p.range[j])
}
mtext(text="SS",side=1,line=row,at=stretch)
row<-3
for (j in seq(1,length(figure3_table$p.range))){
  mtext(text=paste0("(I) ",
                    format(round(figure3_table$post.mean.initial[j],digits=3),nsmall=3),
                    " (F) ",
                    format(round(figure3_table$post.mean.final[j],digits=3),nsmall=3)),
        side=1,line=row,at=figure3_table$p.range[j])
}
mtext(text="PM",side=1,line=row,at=stretch)
row<-4
for (j in seq(1,length(figure3_table$p.range))){
  mtext(text=paste0("(I) ",
                    format(round(figure3_table$cov.initial[j],digits=3),nsmall=3),
                    " (F) ",
                    format(round(figure3_table$cov.final[j],digits=3),nsmall=3)),
        side=1,line=row,at=figure3_table$p.range[j])
}
mtext(text="CP",side=1,line=row,at=stretch)

points(figure3_table$p.range[seq(1,length(figure3_table$p.range))],
       figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))]
       ,pch=20)

text(figure3_table$p.range[seq(1,length(figure3_table$p.range))],
     figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))],
     labels=format(round(
       figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))],
       digits=3),nsmall=3),
     pos=3)

legend("right",#text.width=0.05,
       legend=c("Stop Early for Efficacy",
                "Stop Early for Futility"),
       lty=c('solid','longdash'))

# row<-2
# for (j in seq(1,length(figure3$p.range),by=9)){
#   mtext(text=paste0(format(round(figure3$eff.mon.initial[j],digits=3),nsmall=3)),
#         side=1,line=row,at=figure3$p.range[j])
# }
# mtext(text="EFF",side=1,line=row,at=stretch)
# row<-3
# for (j in seq(1,length(figure3$p.range),by=9)){
#   mtext(text=paste0(format(round(figure3$fut.mon.initial[j],digits=3),nsmall=3)),
#         side=1,line=row,at=figure3$p.range[j])
# }
# mtext(text="FUT",side=1,line=row,at=stretch)
# row<-4
# for (j in seq(1,length(figure3$p.range),by=9)){
#   mtext(text=paste0(format(round(1-figure3$fut.mon.initial[j]-figure3$eff.mon.initial[j],
#                                  digits=3),nsmall=3)),
#         side=1,line=row,at=figure3$p.range[j])
# }
# mtext(text="INC",side=1,line=row,at=stretch)