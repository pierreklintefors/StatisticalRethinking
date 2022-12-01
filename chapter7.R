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


#R code 7.16
N <- 20
kseq <- 1:5
dev <- sapply( kseq , function(k) {
  print(k);
  r <- mcreplicate( 1e2 , sim_train_test( N=N, k=k ) , mc.cores=10 ) # For mac and Linux
  c( mean(r[1,]) , mean(r[2,]) , sd(r[1,]) , sd(r[2,]) )
} )


#R code 7.18
plot( 1:5 , dev[1,] , ylim=c( min(dev[1:2,])-5 , max(dev[1:2,])+10 ) ,
      xlim=c(1,5.1) , xlab="number of parameters" , ylab="deviance" ,
      pch=16 , col=rangi2 )
mtext( concat( "N = ",N ) )
points( (1:5)+0.1 , dev[2,] )
for ( i in kseq ) {
  pts_in <- dev[1,i] + c(-1,+1)*dev[3,i]
  pts_out <- dev[2,i] + c(-1,+1)*dev[4,i]
  lines( c(i,i) , pts_in , col=rangi2 )
  lines( c(i,i)+0.1 , pts_out )
}

#### WAIC calcullation ###########
#R code 7.19

data(cars)
m <- quap(
  alist(
    dist ~ dnorm(mu,sigma),
    mu <- a + b*speed,
    a ~ dnorm(0,100),
    b ~ dnorm(0,10),
    sigma ~ dexp(1)
  ) , data=cars )
set.seed(94)
post <- extract.samples(m,n=1000)

#The likelihood of each observation i at each sample s from the posterior
#R code 7.20

n_samples <- 1000
logprob <- sapply( 1:n_samples ,
         function(s) {
           mu <- post$a[s] + post$b[s]*cars$speed   
           dnorm( cars$dist , mu , post$sigma[s] , log=TRUE )
         } )
#Log of the average by using the log_sum_exp function
#R code 7.21
n_cases <- nrow(cars)
lppd <- sapply( 1:n_cases , function(i) log_sum_exp(logprob[i,]) - log(n_samples) )

#Sum(lppd) will give lppd as defined in text

#The penalty term pWAIC (effective number of parameters/ penalty of the log probability)
#R code 7.22
pWAIC <- sapply( 1:n_cases , function(i) var(logprob[i,]) )

#WAIC = -2(lppd - pWAIC)
##R code 7.23
-2*(sum(lppd) - sum(pWAIC))
WAIC(m, n = n_samples)
#R code 7.24
waic_vec <- -2*( lppd - pWAIC )
sqrt( n_cases*var(waic_vec) )


# Run m6.6, m6.7, m6.8 from chapter 6
#R code 7.25
set.seed(11)
WAIC( m6.7 )

#Compare the models with convient functino
#R code 7.26
set.seed(77)
compare( m6.6 , m6.7 ,m6.8 , func=WAIC )
compare( m6.6 , m6.7 ,m6.8 , func=PSIS )
View(compare)

#Comparing the standard error of differences using pointwise breakdown
#R code 7.27
set.seed(91)
waic_m6.7 <- WAIC( m6.7 , pointwise=TRUE )$WAIC
waic_m6.8 <- WAIC( m6.8 , pointwise=TRUE )$WAIC
n <- length(waic_m6.7)
diff_m6.7_m6.8 <- waic_m6.7 - waic_m6.8
sqrt( n*var( diff_m6.7_m6.8 ) )
    
# a 99 % interval based on this standard errr of differnce
#R code 7.28 
40.0 + c(-1,1)*10.4*2.6


#R code 7.29
plot( compare( m6.6 , m6.7 , m6.8 ) )


#Standard error of difference between the intercept model and model with treatment
# bu not fungus

#R code 7.30
set.seed(92)
waic_m6.6 <- WAIC( m6.6 , pointwise=TRUE )$WAIC
diff_m6.6_m6.8 <- waic_m6.6 - waic_m6.8
sqrt( n*var( diff_m6.6_m6.8 ) )

#R code 7.31
set.seed(93)
compare( m6.6 , m6.7 , m6.8 )@dSE


#How PSIS and WAIC represent imortance and is affected by outliers
#R code 7.32
#Data and models from chapter 5
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )


m5.1 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

m5.2 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )
m5.3 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

precis(m5.1)
precis(m5.2)
precis(m5.3)

#R code 7.33
set.seed(24071847)
compare( m5.1 , m5.2 , m5.3 , func=PSIS )


#Looking at pointwise to see how the outlisers affect PSIS
#R code 7.34
set.seed(24071847)
PSIS_m5.3 <- PSIS(m5.3,pointwise=TRUE)
set.seed(24071847)
WAIC_m5.3 <- WAIC(m5.3,pointwise=TRUE)
plot( PSIS_m5.3$k , WAIC_m5.3$penalty , xlab="“PSIS Pareto k”" ,
      ylab="“WAIC penalty”" , col=rangi2 , lwd=2 )

WAIC(m5.3)

#Estimate age and marrige rates inluence on divorce using a studen T distribution
# The Student T distribution as heavier/thicker tails 

#R code 7.35
m5.3t <- quap(
  alist(
    D ~ dstudent( 2 , mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

PSIS(m5.3t)
PSIS(m5.3)

#comparing their posteriors
precis(m5.3)
precis(m5.3t)


############PRACTICE #############################
e1 = "1. Measure uncertainty continously to avoid that a samll change in 
      probalities resukt in a big change in uncertainty
      2. The uncertainty should increase with the number of possible events.
      3. The uncertainties should be additive where the total uncertainty is 
        the sum of all uncertainties"

e2 = -log(.70, 2) # H(p) = -E(log(p))

#e3
#Entropy = H(p) = -E(log(p))
p = c(.2, .25, .25, .30)
H_p = -sum(p*log(p))

#e4
p = c(.33, .33, .33, 0)
H_p = -sum(p*log(p), na.rm = TRUE)

#m1 Write down and compare the definitions of AIC and WAIC. Which of these criteria is most general? 
#Which assumptions are required to transform the more general criterion into a less general one?

DKL_info = " The Kullback-Leibler divergence is the average difference in log probability
between the target (p) and the model (q).
DKL(p,q)=∑pi(log(pi)−log(qi))=∑(pi*log(pi/qi))
The divergence is the additional entropy,  difference between the cross 
  entropy and the entropy of the target disrtibution.
  DKL(p,q) = H(p,q)-H(p)
  H(p,q) = -∑(pi*log(qi)) 

In the average log probability using a whole distribution can be acheived 
with log pointwise predictive density. The AIC can be used as crude estmimate 
for DKL"

m1 ="The Akaike information criteron:
AIC = -2log(L)+ 2*parameters = -2lppd + 2p
where L is the likelihood
Works when;
(1)  The priors are flat or overwhelmed by the likelihood.
(2)  The posterior distribution is approximately multivariate Gaussian.
(3)  The sample size N is much greater than the number of parameters k.

The Widely applicable information criterion:
WAIC(y|dist) = -2(lppd - ∑( var(θ) * log(p(yi|θ))

y is the observation, θ the parameters and ∑( var(θ) * log(p(yi|θ)) is the penalty term and
is somtimes referred to as the effective number of parameters (PWAIC).


Of the AIC and WAIC, the WAIC is the more general one. Makes no assumption
about the shape of the posterior. It provides an approximation of the 
out-of-sample deviance that converges to the cross-validation approximation 
in a large sample. But in a finite sample, it can disagree. It can disagree 
because it has a different target—it isn’t trying to approximate the 
cross-validation score, but rather guess the out-of-sample KL divergence. 
In the large-sample limit, these tend to be the same."


#m2 Explain the difference between model selection and model comparison. 
#What information is lost under model selection?
m2 = "Model comparison is more general than model selection. In selection
the model with the lowest information criteriion is selected. In comparison
several model are used to understand how introducing new varibles influence
predicability as well as using a causal model to infer causal relationship.
Just using the lowest information criterion may lead to the best prediction but
it does not provide information about the causal relationship. It also discards
the information about the relative accurace wich is contained in the possible
differnce among different selction values such as cross-validation, ´
PARETO-SMOOTHED IMPORTANCE SAMPLING CROSS-VALIDATION (PSIS) and WAIC"

#m3
data(cars)
library(rethinking)
m <- quap(
  alist(
    dist ~ dnorm(mu,sigma),
    mu <- a + b*speed,
    a ~ dnorm(0,100),
    b ~ dnorm(0,10),
    sigma ~ dexp(1)
  ) , data=cars )

WAIC(m, 1000)
WAIC(m, 1500)

m3 = "More observations will provide more information. The information
      criterion is used as approximation of the DKL. In large samples it will 
      converge to a better approximation and differ more for samller samples"

#m4 What happens to the effective number of parameters, as measured by PSIS or 
#WAIC, as a prior becomes more concentrated? Why? Perform some experiments, 
#if you are not sure.

m4 = "It should get smaller because the effective number of parameters is dependent
      on the variance of the log of the probablity"

m_b <- quap(
  alist(
    dist ~ dnorm(mu,sigma),
    mu <- a + b*speed,
    a ~ dnorm(0,10),
    b ~ dnorm(0,1),
    sigma ~ dexp(1)
  ) , data=cars )
set.seed(87)

WAIC(m_b)
WAIC(m)
PSIS(m_b)
PSIS(m)

#m5.  Provide an informal explanation of why informative priors reduce overfitting.
m5 = "When the priors are informative the model is less influenced of the
      by the dat in the sample. If the priors are flat with higher uncertainty, 
      a bigger span of paramter values seems possible and the posterior will
      encode more oof the traing sample."

#m6 Provide an informal explanation of why overly informative priors result in underfitting. 
m6 = "If the priors entail to much information the model will be influenced by
      very narrow span of values from the sample and its features will  be missed. "

# H1 
#In 2007, The Wall Street Journal published an editorial (“We’re Number One, Alas”)
#with a graph of corporate tax rates in 29 countries plotted against tax revenue.
# A badly fit curve was drawn in (reconstructed at right), seemingly by hand, to 
# make the argument that the relationship between tax rate and tax revenue increases 
# and then declines, such that higher tax rates can actually produce less tax revenue. 
# I want you to actually fit a curve to these data, found in data(Laffer). Consider 
# models that use tax rate to predict tax revenue. Compare, using WAIC or PSIS, a 
# straight-line model to any curved models you like. What do you conclude about the
# relationship between tax rate and tax revenue?

data(Laffer)

laff

#Linear regression
h1_a = quap(
  alist(
    tax_revenue ~ dnorm(mu, sigma),
    mu <- a + bt * tax_rate,
    a ~ dlnorm(0),
    bt ~ dnorm(0,2),
    sigma~dnorm(0,1)
    
  ), data = Laffer
)

#"nd degree polynominal
h1_b = quap(
  alist(
    tax_revenue ~ dnorm(mu, sigma),
    mu <- a + bt[1] * tax_rate + bt[2]*tax_rate^2 ,
    a ~ dlnorm(0),
    bt ~ dnorm(0,1),
    sigma~dexp(1)
  ), data = Laffer, start=list(bt=rep(0,2))
)

h1_c = quap(
  alist(
    tax_revenue ~ dnorm(mu, sigma),
    mu <- a + bt[1] * tax_rate + bt[2]*tax_rate^2 +  bt[3]*tax_rate^3 ,
    a ~ dlnorm(0),
    bt ~ dnorm(0,1),
    sigma~dexp(1)
  ), data = Laffer, start=list(bt=rep(0,3))
)


models = c(h1_a, h1_b, h1_c)
par(mfrow=c(1,length(models)))
for (model in models)
{
  tax_seq = seq(from= min(Laffer$tax_rate), to= max(Laffer$tax_rate), length.out =100)
  l = link(model, data = list(tax_rate = tax_seq))
  mu <- apply( l , 2 , mean )
  ci <- apply( l , 2 , PI )
  with(Laffer, plot(tax_rate, tax_revenue))
  lines(tax_seq, mu)
  shade(ci, tax_seq)
  
}

par(mfrow=c(1,1))
h1_a_waic = WAIC(h1_a, pointwise = TRUE)
h1_a_waic = WAIC(h1_b, pointwise = T)

plot(compare(h1_a, h1_b, h1_c))
precis(h1_b, depth = 2)

#H2
set.seed(24071847)
PSIS_h1_b <- PSIS(h1_b,pointwise=TRUE)
set.seed(24071847)
WAIC_h1_b <- WAIC(h1_b,pointwise=TRUE)

plot( PSIS_h1_a$k , WAIC_h1_a$penalty , xlab="“PSIS Pareto k”" ,
      ylab="“WAIC penalty”" , col=rangi2 , lwd=2 )


h2 = quap(
  alist(
    tax_revenue ~ dstudent(2, mu, sigma),
    mu <- a + bt[1] * tax_rate + bt[2]*tax_rate^2 ,
    a ~ dlnorm(0),
    bt ~ dnorm(0,1),
    sigma~ dexp(1)
  ), data = Laffer, start=list(bt=rep(0,2))
)

tax_seq = seq(from= min(Laffer$tax_rate), to= max(Laffer$tax_rate), length.out =100)
l = link(h2, data = list(tax_rate = tax_seq))
mu <- apply( l ,2, FUN = mean )
ci <- apply( l , 2 , PI )
with(Laffer, plot(tax_rate, tax_revenue))
lines(tax_seq, mu)
shade(ci, tax_seq)


#H3
p_1 = c(.2,.2,.2,.2,.2)
ent_1 = -sum(p_1*log(p_1))

p_2 = c(.8, .1, .05, .025 ,.025)
ent_2 =-sum(p_2*log(p_2))

p_3 = c(0.05, 0.15, 0.7,  0.05, 0.05)
ent_3 = -sum(p_3*log(p_3))


# DKL(p,q) = H(p,q)-H(p)
# H(p,q) = -∑(pi*log(qi)) 

#Island 1 to predict  Island2
H_1_2 = -sum(p_1*log(p_2))
DKL_1_2 = H_1_2 - ent_1
#Island 1 to predict  Island 3
H_1_3 = -sum(p_1*log(p_3))
DKL_1_3 = H_1_3 - ent_1

#Island 2 to predict  Island 1
H_2_1 = -sum(p_2*log(p_1))
DKL_2_1 = H_2_1 - ent_2
#Island 2 to predict  Island 3
H_2_3 = -sum(p_2*log(p_3))
DKL_2_3 = H_2_3 - ent_2

#Island 3 to predict Island 1
H_3_1 = -sum(p_3*log(p_1))
DKL_3_1 = H_2_1 - ent_3
#Island 3 to predict Island 2
H_3_2 = -sum(p_3*log(p_2))
DKL_3_2 = H_3_2 - ent_3

I1_score =DKL_1_2 + DKL_1_3

I2_score =DKL_2_1 + DKL_2_3

I3_score = DKL_3_1 + DKL_3_2

min(I1_score,I2_score,I3_score)

#Island 1 is the best to predict the other because it has the closest approximation
# to the others distribution of birds. It has an equal distribution of the species
# and gives equal probability in its predictions. In contrast Island 2 and 3 has 
# a clear dominating species and few of the others. This make them bad models in
# predicting other Island populations that does not share their skewness. 