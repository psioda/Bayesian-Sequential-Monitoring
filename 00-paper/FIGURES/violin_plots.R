rm(list = ls())

set.seed(20160229)

my_data = data.frame(
  y=c(rnorm(1000), rnorm(1000, 0.5), rnorm(1000, 1), rnorm(1000, 1.5)),
  x=c(rep('a', 2000), rep('b', 2000)),
  m=c(rep('i', 1000), rep('j', 2000), rep('i', 1000))
)

GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}

ggplot(my_data, aes(x, y, fill = m)) + geom_split_violin()

#################################################################################################
## INITIAL CONDITIONS ###########################################################################
#################################################################################################
# lower target
theta_L<-0.15

# upper target
theta_H<-0.45

# tail probabilities for priors (low, high)
alpha_L<-0.05
alpha_H<-0.05
#################################################################################################
## PRIOR SPECIFICATION ##########################################################################
#################################################################################################
## Parameterize prior models
# Step 1: Create grid for possible values of phi
# E.g. beta parameterization for middle prior: beta(pi_0*phi,(1-pi_0)*phi)
phi_seq<-seq(0,100,by=0.01)

# Step 2: Compute tail probabilities for every possible choice of phi
# upper tail probability equal to alpha_L/2
test_L<-qbeta(alpha_L/2,(theta_L)*phi_seq,(1-(theta_L))*phi_seq,lower.tail=FALSE)
# lower tail probability equal to alpha_H/2
test_H<-qbeta(alpha_H/2,(theta_H)*phi_seq,(1-(theta_H))*phi_seq,lower.tail=TRUE)

# Step 3: Grid search to find value of phi with the desired tail probability for the priors
phi_L<-phi_seq[which.min(abs(theta_H-test_L))] # fixed 5/13/19
phi_H<-phi_seq[which.min(abs(theta_L-test_H))] # fixed 5/13/19

# Step 4: Find parameters for the priors
alpha_L<-(theta_L)*phi_L
beta_L<-(1-(theta_L))*phi_L

alpha_H<-(theta_H)*phi_H
beta_H<-(1-(theta_H))*phi_H
#################################################################################################
## SIMULATIONS ##################################################################################
#################################################################################################
x<-seq(0,1,length.out=1000)
my_data2 = data.frame(
  y=c(qbeta(x,alpha_L,beta_L),qbeta(x,alpha_L+5,beta_L+15),qbeta(x,alpha_L+10,beta_L+30),
      qbeta(x,alpha_H,beta_H),qbeta(x,alpha_H+5,beta_H+15),qbeta(x,alpha_H+10,beta_H+30)),
  m=c(rep('l', 3000), rep('h', 3000)),
  x=c(rep('1', 1000), rep('2', 1000),rep('3',1000),
      rep('1', 1000), rep('2', 1000),rep('3',1000))
)
#################################################################################################
#################################################################################################
#################################################################################################

ggplot(my_data2, aes(x, y, fill = m)) + geom_split_violin()

