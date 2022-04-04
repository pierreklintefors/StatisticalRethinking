

## 1 Construct a linear regression of weight as predicted by height, using the
# adults (age 18 or greater) from the Howell1 dataset. The heights listed below
# were recorded in the !Kung census, but weights were not recorded for these
# individuals. Provide predicted weights and 89% compatibility intervals for
# each of these individuals. That is, fill in the table below, using model-based
# predictions.
# Individual  height    expected weight     89% interval
#   1           140
#   2           160
#   3           175

library(rethinking)
data("Howell1")
d = Howell1[Howell1$age>=18, ]

height_mean = mean(d$height)

model <- quap(
  alist(
    weight ~dnorm(mu, sigma),
    mu <- a + b*(height-height_mean),
    a ~ dnorm(50, 80),
    b ~ dlnorm(0, 1),
    sigma ~ dunif(0, 30)
  ), data = d
)

post <- extract.samples(model)

individual_heights = c(140, 160, 175) 

for (i in 1:length(individual_heights)) {
  assign(paste("mu_at_", individual_heights[i], sep = ""), post$a + post$b*(individual_heights[i]- height_mean) )
  assign(paste("PI_for_", individual_heights[i], sep = ""), PI(paste("mu_at_", individual_heights[i],sep = "")))
  }

