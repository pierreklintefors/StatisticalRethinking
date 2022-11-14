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
lines( D_seq , mu_mean , lwd=2 )
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

#R code 5.36
precis( m5.5 )


#R code 5.37
xseq <- seq( from=min(dcc$N)-0.15 , to=max(dcc$N)+0.15 , length.out=30 )
mu <- link( m5.5 , data=list(N=xseq) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( K ~ N , data=dcc )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )

# New bivariate model of calories with the logarithmized mass (M) as predictor.
#The log is used because it translates the relationship to the magnitude of mass
#R code 5.38

m5.6 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=dcc )
precis(m5.6)

xseq <- seq( from=min(dcc$M)-0.15 , to=max(dcc$M)+0.15 , length.out=30 )
mu <- link( m5.6 , data=list(M=xseq) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( K ~ M , data=dcc )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )


#New model with both neocortex proportional size and log of of body mass as predictors
#R code 5.39
m5.7 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bN ~ dnorm( 0 , 0.5 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=dcc )
precis(m5.7)

#Visualize and comparing the posteriors of the models
#R code 5.40
plot( coeftab( m5.5 , m5.6 , m5.7 ) , pars=c("bM","bN") )

#Plot the correlation between kilocalories per gram, log(body mass) and neocortex
pairs( ~K + M + N , dcc )


#R code 5.41

xseq <- seq( from=min(dcc$M)-0.15 , to=max(dcc$M)+0.15 , length.out=30 )
mu <- link( m5.7 , data=data.frame( M=xseq , N=0 ) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( NULL, xlim=range(dcc$M) , ylim=range(dcc$K) )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )


#Predictions of neocortex effect on milk energy
xseq <- seq( from=min(dcc$N)-0.15 , to=max(dcc$N)+0.15 , length.out=30 )
mu <- link( m5.7 , data=data.frame( N=xseq , M=0 ) )
mu_mean <- apply(mu,2,mean)
mu_PI <- apply(mu,2,PI)
plot( NULL , xlim=range(dcc$N) , ylim=range(dcc$K) )
lines( xseq , mu_mean , lwd=2 )
shade( mu_PI , xseq )

# DAG where M->K<-N
#R code 5.42
# M -> K <- N
# M -> N
n <- 100
M <- rnorm( n )
N <- rnorm( n , M )
K <- rnorm( n , N - M )
d_sim <- data.frame(K=K,N=N,M=M)


#Using the simulated data in the models
#

m5.5b <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N ,
    a ~ dnorm( 0 , 0.2 ) , # intercept closer to zero
    bN ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim )

m5.6b <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim )

m5.7b <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bN ~ dnorm( 0 , 0.5 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim )

# The slopes in the m5.7b becomes extreme
precis(m5.5b)
precis(m5.6b)
precis(m5.7b)
plot( coeftab( m5.5b , m5.6b , m5.7b ) , pars=c("bM","bN") )


#Two other DAG simulated data set
#R code 5.43
# M -> K <- N
# N -> M
n <- 100
N <- rnorm( n )
M <- rnorm( n , N )
K <- rnorm( n , N - M )
d_sim2 <- data.frame(K=K,N=N,M=M)


m5.5c <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N ,
    a ~ dnorm( 0 , 0.2 ) , # intercept closer to zero
    bN ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim2 )

m5.6c <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim2 )

m5.7c <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bN ~ dnorm( 0 , 0.5 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim2 )

# The slopes in t
precis(m5.5c)
precis(m5.6c)
precis(m5.7c)
plot( coeftab( m5.5c , m5.6c , m5.7c ) , pars=c("bM","bN") )

# M -> K <- N
# M <- U -> N
n <- 100
U <- rnorm( n )
N <- rnorm( n , U )
M <- rnorm( n , U )
K <- rnorm( n , N - M )
d_sim3 <- data.frame(K=K,N=N,M=M)


m5.5d <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N ,
    a ~ dnorm( 0 , 0.2 ) , # intercept closer to zero
    bN ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim3 )

m5.6d <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim3 )

m5.7d <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bN*N + bM*M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bN ~ dnorm( 0 , 0.5 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d_sim3 )

# The slopes in the m5.7b becomes extreme
precis(m5.5d)
precis(m5.6d)
precis(m5.7d)
plot( coeftab( m5.5d , m5.6d , m5.7d ) , pars=c("bM","bN") )

## The Markov equivalence set can be computed using dagitty
#R code 5.44
dag5.7 <- dagitty( "dag{
M -> K <- N
M -> N }" )
coordinates(dag5.7) <- list( x=c(M=0,K=1,N=2) , y=c(M=0.5,K=1,N=0.5) )
MElist <- equivalentDAGs(dag5.7)

drawdag(MElist)


############# Categories #################
#R code 5.45
data(Howell1)
d <- Howell1
str(d)

#If a prior is assigned to the difference in the criterion variable
#it assumes that one of the categories is more uncertain than the other

#The prior distribution of male and females
#R code 5.46
mu_female <- rnorm(1e4,178,20)
mu_male <- rnorm(1e4,178,20) + rnorm(1e4,0,10)
precis( data.frame( mu_female , mu_male ) ) #This makes the interval for males wider

#R code 5.47
d$sex <- ifelse( d$male==1 , 2 , 1 ) # making female=1 and male=2
str( d$sex )
#Now, the same prior can be used for both categories, by indexing a based on sex
m5.8 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a[sex] ,
    a[sex] ~ dnorm( 178 , 20 ) ,
    sigma ~ dunif( 0 , 50 )
  ) , data=d )
precis( m5.8 , depth=2 )

#Expected difference by extracting samples
#R code 5.49
post <- extract.samples(m5.8)
post$diff_fm <- post$a[,1] - post$a[,2]
precis( post , depth=2 )

#Making an index variable of multiple categories from milk data
#R code 5.50
data(milk)
d <- milk
levels(d$clade)

#R code 5.51
d$clade_id <- as.integer( d$clade )

#R code 5.52
d$K <- standardize( d$kcal.per.g )
m5.9 <- quap(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a[clade_id],
    a[clade_id] ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=d )

labels <- paste( "a[" , 1:4 , "]:" , levels(d$clade) , sep="" )
plot( precis( m5.9 , depth=2 , pars="a" ) , labels=labels ,
      xlab="expected kcal (std)" )

#R code 5.53
set.seed(63)
d$house <- sample( rep(1:4,each=8) , size=nrow(d) )

# R code 5.54
m5.10 <- quap(
  alist(
    K ~ dnorm( mu , sigma ),
    mu <- a[clade_id] + h[house],
    a[clade_id] ~ dnorm( 0 , 0.5 ),
    h[house] ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ) , data=d )

precis(m5.10, depth = 2)
labels <- paste( "a[" , 1:4 , "]:" , c("GRyffindor", "Hufflepuff", "Ravenclaw",
                                       "Slytherin") , sep="" )
plot( precis( m5.10 , depth=2 , pars="a" ) , labels=labels ,
      xlab="expected kcal (std)" )

############# PRACTICE ####################################
#EASY
e1 = "μi=βxxi+βz & μi=α+βxxi+βzzi"
e2 = "AD ~Normal(mu_i,sigma),
      mu_i = a + BL * Li + BP * Pi"
   
e3 = "PhdT ~ Normal(mu_i, sigma), 
      mu_i = a + BS*Si + BF*Fi"
      
e4 = "
(1)  μi = α + βAAi + βBBi + βDDi
(2)  μi = α + βAAi + βBBi + βCCi + βDDi
(3)  μi = α + βBBi + βCCi + βDDi
(4)  μi = αAAi + αBBi + αCCi + αDDi
(5)  μi = αA(1−Bi−Ci−Di) + αBBi+αCCi + αDDi 
Model 1 and 3-5 are inferentially equivalent.
Models 1 and 3 both make use of 3 of the 4 total indicator
variables which means that we can always derive the 
parameter estimates for the 4th indicator variable from 
a combination of the three parameter estimates present 
as well as the intercept. Model 4 is akin to an index 
variable approach which is inferentially the same as 
an indicator approach (Models 1 and 3). Model 5 is the
same as Model 4, so long as we assume that each observation 
has to belong to one of the four indicator variables."

#M1

#Salary and free time has both spurious correlation with happiness  

#Monthly salary in thousand dollars
S = rnorm(1000, 3, 0.5)

#Free time (hours/week)
FT = rnorm(1000, 10-S, 2) 

#Happiness 
H = S*FT

data = data.frame(S, FT, H)

data$S = standardize(S)
data$H = standardize(H)
data$FT = standardize(FT)

cor(FT, H)
cor(S, H)
cor(S,FT)

m1S <- quap(
  alist(
    H ~ dnorm(mu, sigma),
    mu <- a + BS * S*-1,
    a ~ dnorm(0,0.3),
    BS ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = data
)

m1FT <- quap(
  alist(
    H ~ dnorm(mu, sigma),
    mu <- a + BFT * FT,
    a ~ dnorm(0,0.3),
    BFT ~ dnorm(0,1),
    sigma ~ dexp(1)
  ), data = data
)


m1SFT <- quap(
  alist(
    H ~ dnorm(mu, sigma),
    mu <- a + BFT * FT + BS * S*-1,
    a ~ dnorm(0,0.3),
    BFT ~ dnorm(0,0.1),
    BS ~ dnorm(0,0.1),
    sigma ~ dexp(1)
  ), data = data
)

plot(coeftab(m1S, m1FT, m1SFT), par= c("BS", "BFT"))

#H1
cor(d$D, d$M)
cor(d$D, d$A)
cor(d$M, d$A)

#DAG
mad_dag <- dagitty("dag{M -> A -> D}")
impliedConditionalIndependencies(mad_dag)
equivalentDAGs(mad_dag)

