######### Example code from the book

library(rethinking)

# R code 4.1
# Normal distribution
pos <- replicate( 1000 , sum( runif(16,-1,1) ) )


hist(pos)

plot(density(pos))


#R code 4.2
#Random growth rate
prod( 1 + runif(12,0,0.1) ) #product of vector elements


#R code 4.3
#Distributions of growth rate

# Generate 10000 of these products of vector elements
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )
dens( growth , norm.comp=TRUE )


#Multiplying small effects can be approximated by addition

#R code 4.4
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )

dens(big, norm.comp = TRUE) # non-Gaussian
dens(small, norm.comp = TRUE) #Gaussian

#R code 4.5
#Large deviates produces Gaussian if converted to log-scale
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )
dens(log.big, norm.comp = TRUE)



# R code 4.6
# Bayes theorem with grid approximation with binomial distribution
w <- 6; n <- 9; # Water outcomes and number of tosses 
P_grid <- seq(from=0,to=1,length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)



#### Making a model of height with Bayesian approxiamtion

#R code 4.7
library(rethinking)
data(Howell1)
d <- Howell1

#4.8
str(d)
#4.9 summary of the data 
precis(d)



#R code 4.11
#Sub-setting data of adults
d2 <- d[ d$age >= 18 , ]

dens(d2$height)


#R code 4.12
# Prior of the mean height 
curve( dnorm( x , 178 , 20 ) , from=100 , to=250 )

#R code 4.13
#Flat Uniform prior of the variance
curve( dunif( x , 0 , 50 ) , from=-10 , to=60 )


#R code 4.14

#Prior predictive simulation to test different priors
sample_mu <- rnorm( 1e4 , 178 , 20 )
sample_sigma <- runif( 1e4 , 0 , 50 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )
 
abline(v = mean(prior_h))

# code 4.15
# larger spread
sample_mu <- rnorm( 1e4 , 178 , 100 )
prior_h <- rnorm( 1e4 , sample_mu , sample_sigma )
dens( prior_h )
abline(v=mean(prior_h))
