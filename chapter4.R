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


#R code 4.16
# Grid approximation of a more complex model with several priors
# Some of this code has not been explained yet
mu.list <- seq( from=150, to=160 , length.out=100 )
sigma.list <- seq( from=7 , to=9 , length.out=100 )
post <- expand.grid( mu=mu.list , sigma=sigma.list )
post$LL <- sapply( 1:nrow(post) , function(i) sum(
  dnorm( d2$height , post$mu[i] , post$sigma[i] , log=TRUE ) ) )
post$prod <- post$LL + dnorm( post$mu , 178 , 20 , TRUE ) +
  dunif( post$sigma , 0 , 50 , TRUE )
post$prob <- exp( post$prod - max(post$prod) )

#R code 4.17
contour_xyz( post$mu , post$sigma , post$prob )

#Or you can plot a simple heat map with:
 #R code 4.18
image_xyz( post$mu , post$sigma , post$prob )

#R code 4.19
#Randomly samples rows from the posterior and pull out parameter values
sample.rows <- sample( 1:nrow(post) , size=1e4 , replace=TRUE ,
                       
                       prob=post$prob )
sample.mu <- post$mu[ sample.rows ]
sample.sigma <- post$sigma[ sample.rows ]


#R code 4.20
# Plot the samples, col.alpha makes the colour transparent
plot( sample.mu , sample.sigma , cex=0.6 , pch=20 , col=col.alpha(rangi2,0.3) )

#R code 4.21
#The shapes of the marginal distribution
dens( sample.mu )
dens( sample.sigma )


#To summarize the widths of these densities with posterior compatibility intervals:
#R code 4.22
PI( sample.mu )
PI( sample.sigma )

#Same analysis with grid approximation but with a small subset of the data
#R code 4.23
d3 <- sample( d2$height , size=20 )

#R code 4.24
mu.list <- seq( from=150, to=170 , length.out=200 )
sigma.list <- seq( from=4 , to=20 , length.out=200 )
post2 <- expand.grid( mu=mu.list , sigma=sigma.list )
post2$LL <- sapply( 1:nrow(post2) , function(i)
  sum( dnorm( d3 , mean=post2$mu[i] , sd=post2$sigma[i] ,
              log=TRUE ) ) )
post2$prod <- post2$LL + dnorm( post2$mu , 178 , 20 , TRUE ) +
  dunif( post2$sigma , 0 , 50 , TRUE )
post2$prob <- exp( post2$prod - max(post2$prod) )
sample2.rows <- sample( 1:nrow(post2) , size=1e4 , replace=TRUE ,
                        prob=post2$prob )
sample2.mu <- post2$mu[ sample2.rows ]
sample2.sigma <- post2$sigma[ sample2.rows ]
plot( sample2.mu , sample2.sigma , cex=0.5 ,
      col=col.alpha(rangi2,0.1) ,
      xlab="mu" , ylab="sigma" , pch=16 )

#R code 4.25
dens( sample2.sigma , norm.comp=TRUE )


### Using quadratic approximation ######

#R code 4.26
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[ d$age >= 18 , ]


#R code 4.27
# List the parameters of the model
flist <- alist(
  height ~ dnorm( mu , sigma ) ,
  mu ~ dnorm( 178 , 20 ) ,
  sigma ~ dunif( 0 , 50 )
)

#R code 4.28
m4.1 <- quap( flist , data=d2 )

#R code 4.29
precis( m4.1 )


#The starting point of quad is random by default but can be specified
#R code 4.30
# object in list gets evaluated/executed while objects in alist do not
start <- list(
  mu=mean(d2$height),
  sigma=sd(d2$height))
m4.1 <- quap( flist , data=d2 , start=start )

#R code 4.31
m4.2 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu ~ dnorm( 178 , 0.1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=d2 )
precis( m4.2 )


#The quad approximation is a multidimensional Gaussian
#A list of means together with matrices of variances and covariances can be used
# to describe the multidimensional Gaussian
#R code 4.32
vcov( m4.1 )

#The variance-covariance matrix above can be factored into two elements
#R code 4.33
diag( vcov( m4.1 ))
cov2cor( vcov( m4.1 ) )


#R code 4.34
# To sample a multidimensional Guassian, vectors of parameter values need to be 
# extracted
library(rethinking)
post <- extract.samples( m4.1 , n=1e4 )
head(post)

# R code 4.35 
precis(post)

plot(post)


#The extract.samples function simulates multivariate distributions 
# This can also be done directly by using mvrnorm()

#R code 4.36
library(MASS)
post <- mvrnorm( n=1e4 , mu=coef(m4.1) , Sigma=vcov(m4.1) )

plot(post)
head(precis(post))


#### Linear prediction ########

#R code 4.37
# Scatter lot of weight and height
library(rethinking)
data(Howell1); d <- Howell1; d2 <- d[ d$age >= 18 , ]
plot( d2$height ~ d2$weight )


#Simulating the priors of intercept and slope
#R code 4.38
set.seed(2971)
N <- 100 # 100 lines
a <- rnorm( N , 178 , 20 )
b <- rnorm( N , 0 , 10 )

#R code 4.39
plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400), xlab="weight" , ylab="height" )
abline( h=0 , lty=2 )
abline( h=272 , lty=1 , lwd=0.5 )
mtext( "b ~ dnorm(0,10)" )
xbar <- mean(d2$weight)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar), 
                        from=min(d2$weight) , to=max(d2$weight) , add=TRUE ,
                        col=col.alpha("black",0.2) )

## The prior can be restricted to positive values to get a more realistic model
# This can be done by sampling them from a normal distribution of log values
#R code 4.40
b <- rlnorm( 1e4 , 0 ,1)
dens( b , xlim=c(0,5) , adj=0.1 )

#Simulate and plot again
#R code 4.41
set.seed(2971)
N <- 100# 100 lines
a <- rnorm( N , 178 , 20 )
b <- rlnorm( N , 0 , 1 )

plot( NULL , xlim=range(d2$weight) , ylim=c(-100,400), xlab="weight" , ylab="height" )
abline( h=0 , lty=2 )
abline( h=272 , lty=1 , lwd=0.5 )
mtext( "b ~ dnorm(0,10)" )
xbar <- mean(d2$weight)
for ( i in 1:N ) curve( a[i] + b[i]*(x - xbar), 
                        from=min(d2$weight) , to=max(d2$weight) , add=TRUE ,
                        col=col.alpha("black",0.2) )


#The new b prior can be put in the quadratic approximation model
# Conventionally <- is used in R instead of = when defining a linear model

#R code 4.42
# load data again, since it's a long way back
library(rethinking)
data(Howell1); d <- Howell1; d2 <- d[ d$age >= 18 , ]

# define the average weight, x-bar
xbar <- mean(d2$weight)

# fit model
m4.3 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*( weight - xbar ) ,
    a ~ dnorm( 178 , 20 ) ,
    b ~ dlnorm( 0 , 1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=d2 )


#R code 4.43
#The log of the b parameter can also be drawn from normal and than used in the linear
# equation with the exp(b) because b = exp(log(b))

m4.3b <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + exp(log_b)*( weight - xbar ),
    a ~ dnorm( 178 , 20 ) ,
    log_b ~ dnorm( 0 , 1),
    sigma ~ dunif( 0 , 50 )
  ) , data=d2 )

precis(m4.3)


#R code 4.45
# Variance-covariance matrix
round( vcov( m4.3 ) , 3 )


#R code 4.46
# Plotting raw data and a single line
plot( height ~ weight , data=d2 , col=rangi2 )
post <- extract.samples( m4.3 )
a_map <- mean(post$a)
b_map <- mean(post$b)
curve( a_map + b_map*(x - xbar) , add=TRUE )

#The single line shows the most probable line but does not entail the uncertainty

#R code 4.47
#Extraction samples of lines
post <- extract.samples( m4.3 )
post[1:5,]

#R code 4.47
post <- extract.samples( m4.3 )
post[1:5,]
#Extracting the first 10 from the data and drawing many lines

#R code 4.48
N <- 200
dN <- d2[ 1:N , ]
mN <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + b*( weight - mean(weight) ) ,
    a ~ dnorm( 178 , 20 ) ,
    b ~ dlnorm( 0 , 1 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=dN )

#Now let’s plot 20 of these lines, to see what the uncertainty looks like.

#R code 4.49

# extract 20 samples from the posterior
post <- extract.samples( mN , n=20 )

# display raw data and sample size
plot( dN$weight , dN$height ,
      xlim=range(d2$weight) , ylim=range(d2$height) ,
      col=rangi2 , xlab="weight" , ylab="height" )

mtext(concat("N = ",N))

# plot the lines, with transparency
for ( i in 1:20 )
  curve( post$a[i] + post$b[i]*(x-mean(dN$weight)) ,
          col=col.alpha("black",0.3) , add=TRUE )


# Showing uncertainty

# Estimating values of length for a 50 kg person
#R code 4.50
post <- extract.samples( m4.3 )
mu_at_50 <- post$a + post$b * ( 50 - xbar )


# R code 4.51
dens( mu_at_50 , col=rangi2 , lwd=2 , xlab="mu|weight=50" )
abline(v=mean(mu_at_50))


#R code 4.52
PI( mu_at_50 , prob=0.89 )


#The link function can be used to this for every case of the data
# and sample from the posterior distribution

# R code 4.53
mu <- link( m4.3 )
str(mu)

# Make a distribution for each weight in sample
#R code 4.54

# define sequence of weights to compute predictions for
# these values will be on the horizontal axis
weight.seq <- seq( from=25 , to=70 , by=1 )

# use link to compute mu
# for each sample from posterior
# and for each weight in weight.seq
mu <- link( m4.3 , data=data.frame(weight=weight.seq) )
str(mu)

# Summarise the distribution for each weight value

#R code 4.56
# summarize the distribution of mu
mu.mean <- apply( mu , 2 , mean )
mu.PI <- apply( mu , 2 , PI , prob=0.89 )

plot(mu.mean)
str(mu.PI)


#Plot the data with line and shade indicating uncertainty for each weight value
#R code 4.57

# plot raw data
# fading out points to make line and interval more visible
plot( height ~ weight , data=d2 , col=col.alpha(rangi2,0.5) )

# plot the MAP line, aka the mean mu for each weight
lines( weight.seq , mu.mean )

# plot a shaded region for 89% PI
shade(mu.PI , weight.seq)


# The link function is just iterating for each case in the data
#The same can be done manually

#R code 4.58

post <- extract.samples(m4.3)
mu.link <- function(weight) post$a + post$b*( weight - xbar )
weight.seq <- seq( from=25 , to=70 , by=1 )
mu <- sapply( weight.seq , mu.link )
mu.mean <- apply( mu , 2 , mean )
mu.CI <- apply( mu , 2 , PI , prob=0.89 )

#R code 4.59
sim.height <- sim( m4.3 , data=list(weight=weight.seq) )
str(sim.height)

# R code 4.60
height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

#R code 4.61

# plot raw data

plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )

# draw MAP line
lines( weight.seq , mu.mean )

# draw HPDI region for line
shade( mu.HPDI , weight.seq )


# draw PI region for simulated heights
shade( height.PI , weight.seq )

#R code 4.62

sim.height <- sim( m4.3 , data=list(weight=weight.seq) , n=1e4 )

height.PI <- apply( sim.height , 2 , PI , prob=0.89 )

#R code 4.61
# plot raw data
plot( height ~ weight , d2 , col=col.alpha(rangi2,0.5) )

# draw MAP line
lines( weight.seq , mu.mean )

# draw HPDI region for line
shade( mu.HPDI , weight.seq )

# draw PI region for simulated heights
shade( height.PI , weight.seq )
