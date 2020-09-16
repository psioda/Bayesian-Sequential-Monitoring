binWidth(60, .733, method = "CP")
binWidth(60, .733, method = "Blaker") # WINNER
binWidth(60, .733, method = "AC")
binWidth(60, .733, method = "Score")
binWidth(60, .733, method = "SOC")
binWidth(60, .733, method = "Wald")

binWidth(n = 60, p = .733)

p  <- 44/60
n  <- 60
se <- sqrt(p*(1-p)/n)
z  <- qnorm(0.975)
p + z*se
p - z*se

p  <- 40/60
n  <- 60
se <- sqrt(p*(1-p)/n)
z  <- qnorm(0.975)
p + z*se
p - z*se
z*se # less than 12% precision
