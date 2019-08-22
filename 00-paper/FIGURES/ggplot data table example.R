rm(list = ls())
library (ggplot2)
library(grid)

df=data.frame(y=c("cat1","cat2","cat3"),x=c(12,10,14),n=c(5,15,20))

p <- ggplot(df, aes(x,y)) + geom_point() +            # Base plot
  theme(plot.margin = unit(c(1,3,1,1), "lines"))   # Make room for the grob
p
i<-1

  p <- p + annotation_custom(
    grob = textGrob(label = df$n[i], hjust = 0, gp = gpar(cex = 1.5)),
    ymin = df$y[i],      # Vertical position of the textGrob
    ymax = df$y[i],
    xmin = 14.3,         # Note: The grobs are positioned outside the plot area
    xmax = 14.3)
  
# Code to override clipping
gt <- ggplot_gtable(ggplot_build(p))
gt$layout$clip[gt$layout$name == "panel"] <- "off"
grid.draw(gt)



for (i in 1:length(df$n))  {
  p <- p + annotation_custom(
    grob = textGrob(label = df$n[i], hjust = 0, gp = gpar(cex = 1.5)),
    ymin = df$y[i],      # Vertical position of the textGrob
    ymax = df$y[i],
    xmin = 14.3,         # Note: The grobs are positioned outside the plot area
    xmax = 14.3)
}    

