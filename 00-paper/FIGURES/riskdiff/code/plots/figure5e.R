plot(c(rd.grid.lower, rd.grid.upper),
     c(enth.lower, enth.upper),
     type = 'l',
     xlim = c(-0.4, 0.4),
     ylim = c(0, 9.25),
     xlab = "",
     ylab = "",
     main = "",
     xaxt = "n",
     yaxt = "n")

title(ylab = "Density Value", line = 1)
title(xlab = "Response Probability", line = 2)

axis(1,
     at = c(delta.skpt, delta.enth),
     labels = c(as.expression(bquote(delta[S])),
              as.expression(bquote(delta[E]))))

polygon(c(rd.grid.lower, 0), c(enth.lower, 0), col = "black")
polygon(c(rd.grid.upper, 0), c(enth.upper, 0), col = "lightgrey")

segments(x0=delta.skpt, y0 = 0, y1 = enth.upper[which(rd.grid.upper==delta.skpt)])
segments(x0=delta.enth, y0 = 0, y1 = enth.upper[which(rd.grid.upper==delta.enth)])

legend('topleft',
       legend= c(as.expression(bquote(P(theta[1]-theta[0]<delta[S])==.(1-sig.eff)))))