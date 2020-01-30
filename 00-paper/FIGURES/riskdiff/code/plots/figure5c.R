plot(grid1d,
    skpt.x,
    type='l',
    xlab="",
    ylab="",
    main="",
    xaxt="n",
    yaxt="n")

axis(1,
    at=c(mu),
    labels=c(as.expression(bquote(mu))))

segments(x0=mu,y0=0,y1=skpt.x[mu*(x.len-1)+1])

title(ylab="Density Value", line=1)
title(xlab="Response Probability",line=2)