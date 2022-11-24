#Overfitting and undefitting

#R code 7.1
sppnames <- c( "afarensis","africanus","habilis","boisei",
              "rudolfensis","ergaster","sapiens")
brainvolcc <- c( 438 , 452 , 612, 521, 752, 871, 1350 )
masskg <- c( 37.0 , 35.5 , 34.5 , 41.5 , 55.5 , 61.0 , 53.5 )
d <- data.frame( species=sppnames , brain=brainvolcc , mass=masskg )


#R code 7.2
d$mass_std <- (d$mass - mean(d$mass))/sd(d$mass)
d$brain_std <- d$brain / max(d$brain)


#R code 7.3
library(rethinking)
m7.1 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ), #exo log to make it positive
    mu <- a + b*mass_std,
    a ~ dnorm( 0.5 , 1 ), # large prior, including negative values, does not make sense
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d )
precis(m7.1)

#You could use lm function and OLS but it does not provide a posterior for
# the standard deviation (sigma)
#R code 7.4
m7.1_OLS <- lm( brain_std ~ mass_std , data=d )
post <- extract.samples( m7.1_OLS )


#Computing R2 with variance the old fashion way (sqaured dev. from the mean)
#R code 7.5
set.seed(12)
s <- sim( m7.1 )
r <- apply(s,2,mean) - d$brain_std
resid_var <- var2(r) # var2 = squared dev. from the mean
outcome_var <- var2( d$brain_std )
1 - resid_var/outcome_var

#function to do the same as code above
#R code 7.6
R2_is_bad <- function( quap_fit ) {
  s <- sim( quap_fit , refresh=0 )
  r <- apply(s,2,mean) - d$brain_std
  1 - var2(r)/var2(d$brain_std)
}

#New model with a secon quadratic predictor
#R code 7.7
m7.2 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ), #b as an vector, same for both predictors´
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,2)) ) # start is needed to tell how long the vector is


#Model with tird degree polynominal 
#R code 7.8

m7.3 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ), # vectors for all predictors
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,3)) )


#Fouth degree polynominal
m7.4 <- quap(
alist(
  brain_std ~ dnorm( mu , exp(log_sigma) ),
  mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
    b[3]*mass_std^3 + b[4]*mass_std^4,
  a ~ dnorm( 0.5 , 1 ),
  b ~ dnorm( 0 , 10 ),
  log_sigma ~ dnorm( 0 , 1 )
), data=d , start=list(b=rep(0,4)) )


#Fith degree
m7.5 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4 +
      b[5]*mass_std^5,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,5)) )

#6th degree with fixes standard deviation
#R code 7.9
m7.6 <- quap(
  alist(
    brain_std ~ dnorm( mu , 0.001 ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4 +
      b[5]*mass_std^5 + b[6]*mass_std^6,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 )
  ), data=d , start=list(b=rep(0,6)) )

#Now to plot each model. We’ll follow the steps from earlier chapters: extract 
#samples from the posterior, compute the posterior predictivedistribution at 
#each of several locations on the horizontal axis, summarize, and plot. For m7.1:
  
#R code 7.10
post <- extract.samples(m7.1)
mass_seq <- seq( from=min(d$mass_std) , to=max(d$mass_std) , length.out=100 )
l <- link( m7.1 , data=list( mass_std=mass_seq ) )
mu <- apply( l , 2 , mean )
ci <- apply( l , 2 , PI )
plot( brain_std ~ mass_std , data=d )
lines( mass_seq , mu )
shade( ci , mass_seq )

#Plot all the models
models = c(m7.1, m7.2, m7.3, m7.4, m7.5, m7.6 )
par(mfrow=c(3,2))
for ( model in models){
  brain_plot(model)
}
View(brain_loo_plot)

#Remove a data point and see how the model's fit is affected
par(mfrow=c(3,2))
for ( model in models){
  brain_loo_plot(model)
}

###########Entropy and accuracy##################
 
#Entropy = H(p) = -E(log(p))

#Entropy for two events with .3 and 7 probability
#H(p) = − (p1log(p1) + p2log(p2))≈0.61

#R code 7.12
p <- c( 0.3 , 0.7 )
-sum( p*log(p) )

#Computing the log of the average probability for each observation
# using log-pointwise-predictive-density (lppd)

#R code 7.13
set.seed(1)
lppd( m7.1 , n=1e4 )

#The coded needed to replicate the code above
#R code 7.14

set.seed(1)
logprob <- sim( m7.1 , ll=TRUE , n=1e4 )#ll=True returns the log-probability
n <- ncol(logprob)
ns <- nrow(logprob)
f <- function( i ) log_sum_exp( logprob[,i] ) - log(ns)
( lppd <- sapply( 1:n , f ) )

#Log score of each of the previous models, same as R2 it is increased with complexity
#R code 7.15
set.seed(1) 
sapply( list(m7.1,m7.2,m7.3,m7.4,m7.5,m7.6) , function(m) sum(lppd(m)) )

