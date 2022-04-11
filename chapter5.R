rm(list = ls())
graphics.off()

#R code 5.1
# load data and copy
library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# standardize variables
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )
d$A <- standardize( d$MedianAgeMarriage )

#R code 5.2
sd( d$MedianAgeMarriage )

#R code 5.3
m5.1 <- quap(
  alist(
    D ~ dnorm( mu , sigma ),
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )) , data = d )

#R code 5.4
set.seed(10)
prior <- extract.prior( m5.1 )
mu <- link( m5.1 , post=prior , data=list( A=c(-2,2) ) )
plot( NULL , xlim=c(-2,2) , ylim=c(-2,2) )
for ( i in 1:50 ) lines( c(-2,2) , mu[i,] , col=col.alpha("black",0.4) )

#R code 5.5
# compute percentile interval of mean
A_seq <- seq( from=-3 , to=3.2 , length.out=30 )
mu <- link( m5.1 , data=list(A=A_seq) )
mu.mean <- apply( mu , 2, mean )
mu.PI <- apply( mu , 2 , PI )

# plot it all
plot( D ~ A , data=d , col=rangi2 )
lines( A_seq , mu.mean , lwd=2 )
shade( mu.PI , A_seq )
precis(m5.1)


#R code 5.6
# Model of marriage rate and divorce
m5.2 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

# compute percentile interval of mean
M_seq <- seq( from=-3 , to=60 , length.out=30 )
mu <- link( m5.2 , data=list(M=M_seq) )
mu.mean <- apply( mu , 2, mean )
mu.PI <- apply( mu , 2 , PI )

# plot it all
plot( D ~ M , data=d , col=rangi2 )
lines( M_seq , mu.mean , lwd=2 )
shade( mu.PI , M_seq )
precis(m5.2)


### DAG ############
#R code 5.7
library(dagitty)
dag5.1 <- dagitty( "dag{ A -> D; A -> M; M -> D }" )
coordinates(dag5.1) <- list( x=c(A=0,D=1,M=2) , y=c(A=0,D=1,M=0) )
drawdag( dag5.1 )

#R code 5.8
#Check if there are any conditional implication (confounding backdoors)?
DMA_dag2 <- dagitty('dag{ D <- A -> M }')
impliedConditionalIndependencies( DMA_dag2 )
drawdag(DMA_dag2)

#This shows the Divorce and Marriage rate should be 
#independent after conditioning on Age

#The first DAG has no conditional independences
#R code 5.9
impliedConditionalIndependencies( dag5.1 )

#The dependence of D and M while conditioned on A can be tested with multiple 
#regression

#Quadratic approximation with multiple slope parameters for divorce data
#R code 5.10

m5.3 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )
precis( m5.3 )


#Visualisation of the posterior of the three models
#R code 5.11
plot( coeftab(m5.1,m5.2,m5.3), par=c("bA","bM") )

#The effect of marriage rate (bM) is around zero when age is included in the model

#Simulation of the DAG
#R code 5.12
N <- 50 # number of simulated States
age <- rnorm( N )# sim A
mar <- rnorm( N , -age ) # sim A -> M
div <- rnorm( N , age ) # sim A -> D

# t is transposing so the varibles becomes the columns
sim_dat = as.data.frame(t(rbind(age, mar, div, deparse.level = 2) ))

sim_dat$D <- standardize( d$Divorce )
sim_dat$M <- standardize( d$Marriage )
sim_dat$A <- standardize( d$MedianAgeMarriage )
#Use these variables in the models 

m5.1b <- quap(
  alist(
    D ~ dnorm( mu , sigma ),
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )) , data = sim_dat )

m5.2b <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = sim_dat )

m5.3b <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = sim_dat )

#The same patterns between the posterior means 
plot( coeftab(m5.1b,m5.2b,m5.3b), par=c("bA","bM") )


#### Visualising ###########

##Predictor residual plots

#Marriage rate is approximated from age
#R code 5.13
m5.4 <- quap(
  alist(
    M ~ dnorm( mu , sigma ) ,
    mu <- a + bAM * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bAM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

#Computing the residuals
#R code 5.14
mu <- link(m5.4)
mu_mean <- apply( mu , 2 , mean )
mu_resid <- d$M - mu_mean
mu_PI <- apply(mu,2, PI)

plot(D~mu_resid, data = d)

D_seq <- seq( from=-2 , to=2 , length.out=50 )
lines( D_seq , mu.mean , lwd=2 )
shade(mu_PI , D_seq )


###Posterior prediction plots
#Simulating data

#R code 5.15
# call link without specifying new data
# so it uses original data
mu <- link( m5.3 )

# summarize samples across cases
mu_mean <- apply( mu , 2 , mean )
mu_PI <- apply( mu , 2 , PI )

# simulate observations
# again no new data, so uses original data
D_sim <- sim( m5.3 , n=1e4 )
D_PI <- apply( D_sim , 2 , PI )

#Plot predictions against observed
#R code 5.16
plot( mu_mean ~ d$D , col=rangi2 , ylim=range(mu_PI) ,
      xlab="Observed divorce" , ylab="Predicted divorce" )
abline( a=0 , b=1 , lty=2 )
for ( i in 1:nrow(d) ) lines( rep(d$D[i],2) , mu_PI[,i] , col=rangi2 )

#Tool for finding lables of points 
#R code 5.17
identify( x=d$D , y=mu_mean , labels=d$Loc )


#Simulating spurious (confounding) associations to see how multiple regression 
# can be used to indicate the right predictor
#R code 5.18
N <- 100# number of cases
x_real <- rnorm( N )# x_real as Gaussian with mean 0 and stddev 1
x_spur <- rnorm( N , x_real )# x_spur as Gaussian with mean=x_real
y <- rnorm( N , x_real )# y as Gaussian with mean=x_real
d <- data.frame(y,x_real,x_spur)# bind all together in data frame

pairs(d)


### Counterfactual plots
#R code 5.19
data(WaffleDivorce)
d <- list()
d$A <- standardize( WaffleDivorce$MedianAgeMarriage )
d$D <- standardize( WaffleDivorce$Divorce )
d$M <- standardize( WaffleDivorce$Marriage )
m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 ),
    ## A -> M
    M ~ dnorm( mu_M , sigma_M ),
    mu_M <- aM + bAM*A,
    aM ~ dnorm( 0 , 0.2 ),
    bAM ~ dnorm( 0 , 0.5 ),
    sigma_M ~ dexp( 1 )
  ) , data = d )

precis(m5.3_A)

#Negative correlation between age and marriage rate, what whould happen when age
# manipulated?

#R code 5.20
A_seq <- seq( from=-2 , to=2 , length.out=30 )

#Simulate
#R code 5.21
# prep data
sim_dat <- data.frame( A=A_seq )

# simulate M and then D, using A_seq, vars tells order
s <- sim( m5.3_A , data=sim_dat , vars=c("M","D") )

#R code 5.22
plot( sim_dat$A , colMeans(s$D) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated A" , ylab="counterfactual D" )
shade( apply(s$D,2,PI) , sim_dat$A )
mtext( "Total counterfactual effect of A on D" )

# A -> M
plot( sim_dat$A , colMeans(s$M) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated A" , ylab="counterfactual M" )
shade( apply(s$D,2,PI) , sim_dat$A )
mtext( "Total counterfactual effect of A on M" )
par(xaxt="s")


# The expected effect of increasing median age from 20 to 30
#R code 5.23
# new data frame, standardized to mean 26.1 and std dev 1.24
sim2_dat <- data.frame( A=(c(20,30)-26.1)/1.24 )
s2 <- sim( m5.3_A , data=sim2_dat , vars=c("M","D") )
mean( s2$D[,2] - s2$D[,1]) 

#Manipulating marriage rate when average age is assumed
# R code 5.24
sim_dat <- data.frame( M=seq(from=-2,to=2,length.out=30) , A=0 )
s <- sim( m5.3_A , data=sim_dat , vars="D" )
plot( sim_dat$M , colMeans(s) , ylim=c(-2,2) , type="l" ,
      xlab="manipulated M" , ylab="counterfactual D" )
shade( apply(s,2,PI) , sim_dat$M )
mtext( "Total counterfactual effect of M on D" )

#Overthinking counterfactuals A-> D and A->M->D
#R code 5.25
A_seq <- seq( from=-2 , to=2 , length.out=30 )

#Extracting posterior samples
#R code 5.26
post <- extract.samples( m5.3_A )
M_sim <- with( post , sapply( 1:30 ,
                              function(i) rnorm( 1e3 , aM + bAM*A_seq[i] , sigma_M ) ) )

#R code 5.27
D_sim <- with( post , sapply( 1:30 ,
                              function(i) rnorm( 1e3 , a + bA*A_seq[i] + bM*M_sim[,i] , sigma ) ) )

plot(A_seq, colMeans(D_sim), ylim= c(-2,2), type="l" )

########### Masked relationships ###############
#R code 5.28
library(rethinking)
data(milk)
d <- milk
str(d)

# Standardising variables
#R code 5.29
d$K <- standardize( d$kcal.per.g )
d$N <- standardize( d$neocortex.perc )
d$M <- standardize( log(d$mass) )

# Bivariate regression with vague priors
#R code 5.30
m5.5_draft <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N ,
    a ~ dnorm( 0 , 1 ) ,
    bN ~ dnorm( 0 , 1 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
#The code above creates error becasuse there are NAs in neocortex variable´


#Making a new data set with no NAs
#R code 5.32
dcc <- d[ complete.cases(d$K,d$N,d$M) , ]


#R code 5.33

m5.5_draft <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N ,
    a ~ dnorm( 0 , 1 ) ,
    bN ~ dnorm( 0 , 1 ) ,
    sigma ~ dexp( 1 )
  ) , data=dcc )


#Simulate and plot regression lines of the priors

#R code 5.34
prior <- extract.prior( m5.5_draft )
xseq <- c(-2,2)
mu <- link( m5.5_draft , post=prior , data=list(N=xseq) )
plot( NULL , xlim=xseq , ylim=xseq )
for ( i in 1:50 ) lines( xseq , mu[i,] , col=col.alpha("black",0.3) )

# Poor priors, can be approved
# R code 5.35
m5.5 <- quap(
alist(
  K ~ dnorm( mu , sigma ) ,
  mu <- a + bN*N ,
  a ~ dnorm( 0 , 0.2 ) , # intercept closer to zero
  bN ~ dnorm( 0 , 0.5 ) ,
  sigma ~ dexp( 1 )
) , data=dcc )

prior <- extract.prior( m5.5 )
xseq <- c(-2,2)
mu <- link( m5.5 , post=prior , data=list(N=xseq) )
plot( NULL , xlim=xseq , ylim=xseq )
for ( i in 1:50 ) lines( xseq , mu[i,] , col=col.alpha("black",0.3) )

