

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
    a ~ dnorm(100, 30),
    b ~ dlnorm(0, 1),
    sigma ~ dnorm(0, 30)
  ), data = d
)

post <- extract.samples(model)

individual_heights = c(140, 160, 175) 

for (i in 1:length(individual_heights)) {
  assign(paste("mu_at_", individual_heights[i], sep = ""), post$a + post$b*(individual_heights[i]- height_mean) )
  assign(paste("PI_for_", individual_heights[i], sep = ""), PI(paste("mu_at_", individual_heights[i],sep = "")))
  }

PI_140 = round(PI(mu_at_140), 2)
PI_160 = round(PI(mu_at_160), 2)
PI_175 = round(PI(mu_at_175),2)

tab =matrix(nrow = 3, ncol = 4)

predicted_heights = c(mean(mu_at_140), mean(mu_at_160), mean(mu_at_175))

tab[,1] = c(1:3)
tab[,2] = individual_heights
tab[,3] = round(predicted_heights, 2) 
tab[,4] = c(paste(PI_140[1],PI_140[2], sep = "-"),
            paste(PI_160[1],PI_160[2], sep = "-"), 
            paste(PI_175[1],PI_175[2], sep = "-"))


colnames(tab) = c("ID", "Height","Pred wieght", "89% PI")

precis(model)
tab


################################
# 2. From the Howell1 dataset, consider only the people younger than 13 years
# old. Estimate the causal association between age and weight. Assume that
# age influences weight through two paths. First, age influences height, and
# height influences weight. Second, age directly influences weight through age-
#   related changes in muscle growth and body proportions. All of this implies
# this causal model (DAG):
# w <- A -> H -> w
# Use a linear regression to estimate the total (not just direct) causal effect of
# each year of growth on weight. Be sure to carefully consider the priors. Try
# using prior predictive simulation to assess what they imply.

d2 = Howell1[Howell1$age<13,]

plot(weight~age+height, d2)


height_mean13 = mean(d2$height)
weight_mean13 = mean(d2$weight)
model2 <- quap(
  alist(
    weight ~ dnorm(mu, sigma),
    mu <- a+ b1*age + b2*(height- height_mean13) + b3*(age*height),
    a ~ dnorm(weight_mean13, 20),
    b1 ~ dlnorm(0, 1),
    b2 ~ dlnorm(0, 1),
    b3 ~ dnorm(0, 2),
    sigma ~ dexp(1)
  ), data = d2
)

age_seq = 0:13
height_seq = seq(50, 170, 13)
mu = link(model2, data = data.frame(age= age_seq, height = height_seq) )

precis(model2)
