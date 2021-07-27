lines(figure3$p.range,
      figure3$eff.mon.initial,
     lwd=2)

lines(figure3$p.range,
      figure3$fut.mon.initial,
      lwd=2,
      lty='longdash')

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

points(figure3_table$p.range[seq(1,length(figure3_table$p.range))],
       figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))]
       ,pch=20)
mtext(text="SS",side=1,line=row,at=stretch)

text(figure3_table$p.range[seq(1,length(figure3_table$p.range))],
     figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))],
     labels=format(round(
       figure3_table$eff.mon.initial[seq(1,length(figure3_table$p.range))],
       digits=3),nsmall=3),
     pos=3)