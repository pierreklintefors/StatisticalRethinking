#R code 3.2

p_grid <- seq( from=0 , to=1 , length.out=1000 )

prob_p <- rep( 1 , 1000 )

prob_data <- dbinom( 6 , size=9 , prob=p_grid )

posterior <- prob_data * prob_p

posterior <- posterior / sum(posterior)

#R code 3.3

samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )



#R code 3.4

plot( samples )


#R code 3.5

library(rethinking)

dens( samples )


#R code 3.6
# add up posterior probability where p < 0.5
sum( posterior[ p_grid < 0.5 ] )

#R code 3.7
sum( samples < 0.5 ) / 1e4

#R code 3.8
sum( samples > 0.5 & samples < 0.75 ) / 1e4

#R code 3.9
quantile( samples , 0.8 )

# R code 3.10
quantile( samples ,c( 0.1 , 0.9 ) )
 
# R code 3.11
p_grid <- seq( from=0 , to=1 , length.out=1000 )

prior <- rep(1,1000)

likelihood <- dbinom( 3 , size=3 , prob=p_grid )

posterior <- likelihood * prior

posterior <- posterior / sum(posterior)

samples <- sample( p_grid , size=1e4 , replace=TRUE , prob=posterior )
 
#  R code 3.12
PI( samples , prob=0.8 )


#R code 3.13
HPDI( samples , prob=0.5 )

#R code 3.14
p_grid[ which.max(posterior) ]

#R code 3.15
chainmode( samples , adj=0.01 )


#R code 3.16
mean( samples )
median( samples )

#R code 3.17
sum( posterior*abs( 0.5 - p_grid ) )

#R code 3.18
loss <- sapply( p_grid , function(d) sum( posterior*abs( d - p_grid ) ) )

#R code 3.19
p_grid[ which.min(loss) ]

#R code 3.20
dbinom( 0:2 , size=2 , prob=0.7 )


#R code 3.21
rbinom( 1 , size=2 , prob=0.7 )

# R code 3.22
rbinom( 10 , size=2 , prob=0.7 )

#R code 3.23
dummy_w <- rbinom( 1e5 , size=2 , prob=0.7 )
table(dummy_w)/1e5
dummy_w

#R code 3.24
dummy_w <- rbinom( 1e5 , size=9 , prob=0.7 )
simplehist( dummy_w , xlab="dummy water count" )

#R code 3.25
w <- rbinom( 1e4 , size=9 , prob=0.6 )


#R code 3.26
w <- rbinom(1e4 , size=9 , prob=samples)

simplehist(w)

####################################################################
######################## Practice ##################################
####################################################################

#R code 3.27

p_grid <- seq( from=0 , to=1 , length.out=1000 )

prior <- rep( 1 , 1000 )

likelihood <- dbinom( 6 , size=9 , prob=p_grid )

posterior <- likelihood * prior

posterior <- posterior / sum(posterior)

set.seed(100)

samples <- sample( p_grid , prob=posterior , size=1e4 , replace=TRUE )

######### EASY ##################
e1 = sum(samples < 0.2)/length(samples)
cat("Posterior probability below 0.2 is: " , e1)
e2 = sum(samples> 0.8)/length(samples)
cat("Posterior probability over 0.8 is: " , e2)
e3 = sum(samples >.2 & samples < .8) / length(samples)
cat("The proportion of the posterior probability between 0.2 and 0.8 is: ", e3)
e4 = quantile(posterior, 0.2)
cat("20% of the posterior is below: ", e4)
e5 = quantile(posterior, 0.8)
cat("20% of the posterior is above: ", e5)

e6 = HPDI(samples, prob=0.6)
e6
e7 = PI(samples, prob =  0.66)

####### Medium  ####################

#M1
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- rep( 1 , 1000 )
likelihood <- dbinom( 8 , size=15 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(101)
plot(p_grid,posterior)

#M2
samples = sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)
HPDI(samples, prob = .9)

#M3
pred_sim = rbinom(n=1e4, size = 15, prob =samples )
simplehist(pred_sim, xlab= "Number of water tosses")
cat("The probability of getting 8 water out of 15 tosses is: ", sum(pred_sim==8) /length(pred_sim))

#M4

pred_sim2 = rbinom(n=1e4, size = 9, prob = samples)
cat("The probability of getting 6 water in 9 tosses is: ", sum(pred_sim2==6)/length(pred_sim2))
simplehist(pred_sim2)

#M5 - The above assignments with new prior
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- ifelse(p_grid< 0.5 , 0, 1 )
likelihood <- dbinom( 8 , size=15 , prob=p_grid )
posterior <- likelihood * prior
posterior <- posterior / sum(posterior)
set.seed(101)
plot(p_grid,posterior)
  #M5.2
  samples = sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)
  HPDI(samples, prob = .9)
  #M5.3
  pred_sim = rbinom(n=1e4, size = 15, prob =samples )
  simplehist(pred_sim, xlab= "Number of water out of 15 tosses ", main = "New prior")
  cat("(The probability of getting 8 water out of 15 tosses is: ", sum(pred_sim==8) /length(pred_sim))
  #M5.4
  
  pred_sim2 = rbinom(n=1e4, size = 9, prob = samples)
  cat("The probability of getting 6 water in 9 tosses is: ", sum(pred_sim2==6)/length(pred_sim2))
  simplehist(pred_sim2, xlab = "Number of water out of 9 tosses", main = "New prior")

#M6
tosses = 15
precision = 1
p_grid <- seq( from=0 , to=1 , length.out=1000 )
prior <- ifelse(p_grid< 0.5 , 0, 1 )
while (precision > 0.05) {
  tosses = tosses +1
  likelihood <- dbinom( 8 , size=tosses , prob=p_grid )
  posterior <- likelihood * prior
  posterior <- posterior / sum(posterior)
  samples = sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)
  precision = HPDI(samples, prob = .99)[2] - HPDI(samples, prob = .99)[1]
  
}

cat("To get a width of .05 of the 99th percentile you ned to toss the globe", tosses, "times", sep = " ")

########## HARD ###############
library(rethinking)
data(homeworkch3)
total_boys = sum(birth1) + sum(birth2)
num_births = length(birth1) + length(birth2)
#H1
p_grid = seq(from=0, to= 1, length.out=1000)
prior = rep(1, length(p_grid))
likelihood = dbinom(total_boys, size = num_births, prob = p_grid)
posterior = prior * likelihood
posterior = posterior / sum(posterior)

p_grid[which.max(posterior)]


#H2
birth_samples = sample(p_grid, size = 1e4, replace = TRUE, prob = posterior)
HPDI(birth_samples, prob = c(.5, .89, .97))

#H3
birth_sim = rbinom(1e4, size = 200, prob = birth_samples)
dens(birth_sim, adj = 0.1, main = "Boy birth simulation")
abline(v= total_boys, col = "red")


#H4
first_born_sim = rbinom(1e4, size = 100, prob = birth_samples)
dens(first_born_sim, adj = .1,  main = "First born boys simulation")
abline(v=sum(birth1), col="red")


#H5
num_girl_first = sum(birth1==0)
boy_after_girl = sum(birth2[birth1==0])
boy_after_girl_sim = rbinom(1e4, size = num_girl_first, prob = birth_samples)
dens(boy_after_girl_sim, adj = 0.1, main = "Boys after girls simulation")
abline(v=boy_after_girl, col="red")
