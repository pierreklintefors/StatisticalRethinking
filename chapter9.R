#The metropolis algorithm
#R code 9.1 #The king Markov example
num_weeks <- 1e5
positions <- rep(0,num_weeks)
current <- 10
for ( i in 1:num_weeks ) {
  ## record current position
  positions[i] <- current
  ## flip coin to generate proposal
  proposal <- current + sample( c(-1,1) , size=1 )
  ## now make sure he loops around the archipelago
  if ( proposal < 1 ) proposal <- 10
  if ( proposal > 10 ) proposal <- 1
  ## move?
  prob_move <- proposal/current
  current <- ifelse( runif(1) < prob_move , proposal , current )
}

#R code 9.2
plot( 1:100 , positions[1:100] )

#R code 9.3
plot( table( positions ) )


#Sample form a highdimensional distribution to will not be near the mode
#R code 9.4
library(rethinking)
D <- 10
T <- 1e3
Y <- rmvnorm(T,rep(0,D),diag(D))
rad_dist <- function( Y ) sqrt( sum(Y^2) )
Rd <- sapply( 1:T , function(i) rad_dist( Y[i,] ) )
dens( Rd )

#Both Metropolis and Gibbs tend to get stuck in a local
#minumum without exploring the total posterior for high dimensional distributions

#Hamiltonian Monte Carlo
#Requires continious parameters
#More complex computation but can represent the posterior in fewer steps
#Uses a physics engine 

#The HMC need five things
# 1 a function that returns the negative log-probability of the data at the current position
# 2. a function that returns the gradient of the negative log-probability (dervivatives)
# 3. a step size, epsilion
# 4. a count of steps, L
# 5. the position, current_q

#R code 9.5
# U needs to return neg-log-probability
U <- function( q , a=0 , b=1 , k=0 , d=1 ) {
  muy <- q[1]
  mux <- q[2]
  U <- sum( dnorm(y,muy,1,log=TRUE) ) + sum( dnorm(x,mux,1,log=TRUE) ) +
    dnorm(muy,a,b,log=TRUE) + dnorm(mux,k,d,log=TRUE)
  return( -U )
}

#The derivative of of the logarithm of any univariate Gaussian with mean a and 
#standard deviation b with respect to a is:
# ∂logN(y|a,b) = y−a/b2
#     ∂a

#R code 9.6

# gradient function
# need vector of partial derivatives of U with respect to vector q
U_gradient <- function( q , a=0 , b=1 , k=0 , d=1 ) {
  muy <- q[1]
  mux <- q[2]
  G1 <- sum( y - muy ) + (a - muy)/b^2 #dU/dmuy
  G2 <- sum( x - mux ) + (k - mux)/d^2 #dU/dmux
  return( c( -G1 , -G2 ) ) # negative bc energy is neg-log-prob
}
# test data
set.seed(7)
y <- rnorm(50)
x <- rnorm(50)
x <- as.numeric(scale(x))
y <- as.numeric(scale(y))

#R code 9.7
library(shape) # for fancy arrows
Q <- list()
Q$q <- c(-0.1,0.2)
pr <- 0.3
plot( NULL , ylab="muy" , xlab="mux" , xlim=c(-pr,pr) , ylim=c(-pr,pr) )
step <- 0.03
L <- 11 # 0.03/28 for U-turns — 11 for working example
n_samples <- 4
path_col <- col.alpha("black",0.5)
points( Q$q[1] , Q$q[2] , pch=4 , col="black" )
for ( i in 1:n_samples ) {
  Q <- HMC2( U , U_gradient , step , L , Q$q )
  if ( n_samples < 10 ) {
    for ( j in 1:L ) {
      K0 <- sum(Q$ptraj[j,]^2)/2 # kinetic energy
      lines( Q$traj[j:(j+1),1] , Q$traj[j:(j+1),2] , col=path_col , lwd=1+2*K0 )
    }
    points( Q$traj[1:L+1,] , pch=16 , col="white" , cex=0.35 )
    Arrows( Q$traj[L,1] , Q$traj[L,2] , Q$traj[L+1,1] , Q$traj[L+1,2] ,
            arr.length=0.35 , arr.adj = 0.7 )
    text( Q$traj[L+1,1] , Q$traj[L+1,2] , i , cex=0.8 , pos=4 , offset=0.4 )
  }
  points( Q$traj[L+1,1] , Q$traj[L+1,2] , pch=ifelse( Q$accept==1 , 16 , 1 ) ,
          col=ifelse( abs(Q$dH)>0.1 , "red" , "black" ) )
}


#Break down of HMC2 functino
#R code 9.8

HMC2 <- function (U, grad_U, epsilon, L, current_q) {
  q = current_q
  p = rnorm(length(q),0,1) # random flick - p is momentum.
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA,nrow=L+1,ncol=length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p
  
  #R code 9.9
  # Alternate full steps for position and momentum
  for ( i in 1:L ) {
    q = q + epsilon * p # Full step for the position
    # Make a full step for the momentum, except at end of trajectory
    if ( i!=L ) {
      p = p - epsilon * grad_U(q)
      ptraj[i+1,] <- p
    }
    qtraj[i+1,] <- q
  }
  
  #R code 9.10
  # Make a half step for momentum at the end
  p = p - epsilon * grad_U(q) / 2
  ptraj[L+1,] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept <- 0
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    new_q <- q # accept
    accept <- 1
  } else new_q <- current_q # reject
  return(list( q=new_q, traj=qtraj, ptraj=ptraj, accept=accept ))
}

######## Ulam using rugged data

#R code 9.11
library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

#Quad aproximation model
#R code 9.12
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )
precis( m8.3 , depth=2 )


#When using ulam (HMC):
#Preprocc all variables before
#Trim down data to only the variable of interest

#R code 9.13 #list instead of Df because elements in a list can be of differetn lenghts
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer( dd$cid )
)
str(dat_slim)


#Using HCM Ulam
#R code 9.14

m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=1 )

#R code 9.15
precis( m9.1 , depth=2 )

#Multiple chains in paralell
#R code 9.16
m9.1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=4 , cores=4 )


#R code 9.17
show( m9.1 )

#R code 9.18
precis( m9.1 , 2 )


#R code 9.19
pairs( m9.1 )


#tace plots
#R code 9.20
traceplot( m9.1 )

#Trace rank plots, histograms of each individual parameter ranked values
#R code 9.21
trankplot( m9.1 )

#Using really flat priors canmake the chain erratically sample extreme parameter values

#R code 9.22
y <- c(-1,1)
set.seed(11)
m9.2 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- alpha ,
    alpha ~ dnorm( 0 , 1000 ) ,
    sigma ~ dexp( 0.0001 )
  ) , data=list(y=y) , chains=3 )


#R code 9.23
precis( m9.2, )

pairs(m9.2@stanfit)

traceplot(m9.2)


#The chains in this case can be fixed with more informative priors

#R code 9.24
set.seed(11)
m9.3 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- alpha ,
    alpha ~ dnorm( 1 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )
precis( m9.3 )

traceplot(m9.3)


#Non-indentifiable model
#Gaussian with mean 0 and sd 1
#R code 9.25
set.seed(41)
y <- rnorm( 100 , mean=0 , sd=1 )


#R code 9.26
set.seed(384)
m9.4 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm( 0 , 1000 ),
    a2 ~ dnorm( 0 , 1000 ),
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )
precis( m9.4 )


#With weakly informative priors
#R code 9.27
m9.5 <- ulam(
  alist(
    y ~ dnorm( mu , sigma ) ,
    mu <- a1 + a2 ,
    a1 ~ dnorm( 0 , 10 ),
    a2 ~ dnorm( 0 , 10 ),
    sigma ~ dexp( 1 )
  ) , data=list(y=y) , chains=3 )
precis( m9.5 )

traceplot(m9.5)


##### PRATICE ############
e1 = "2 The likelihood function must be Gaussian
      3 the proposl distrbution must be symmetric"

e2 = "Gibbs sampling is more effecient by using more clever proposal for next 
      sample. This make it more likely to find a good representaion of the 
      posterior in fewer steps then the metreopolis. Gibbs samplings
      stratery does not work good for high-dimensional posteriors where the Gibbs
      tend to get stuck in local area without exploring the full 
      parameter space."

e3 = "The HCM rewuires continious parameters and does not work on discrete
      parameter values. This is becuase the physcal engine that HCM is buid upon 
      rewuries a continues space to be able to glide on that space. And as any other
      algorithm, the HCM cannot correctly approximate really complex posteriors."

e4 = "The effective number of samples is the length of a hypothetical 
      Markov chain of independetn steps with no autocorrelation that is needed to 
      gain the quality of estmiate as the current chain. The effective number of 
      samples is often a smaller number than the actual samples but it can also 
      bigger. This is often a better meausrument of the quality of the chain than
      the actual number of samples."

e5 = "Rhat should apporach 1, indication convergence of the chains."

#e6
e6 = "A good trace plot is centered around a range of values where it does not 
    go to the extremes. There should also be an ovlap of the chains. A malfunctioning
    chain would have extreme values and do not shoul overlap."
e7 = "In the trace rank plots, the differetn bins should be crossing each other
      where the same chain should not have prololonged maximum ranked sample
      for longer periods, several samples."

#M1
#Re-estimate the terrain ruggedness model from the chapter, but now using a 
# uniform prior for the standard deviation, sigma. The uniform prior should be 
# dunif(0,1). Use ulam to estimate the posterior. Does the different prior have 
# any detectible influence on the posterior distribution of sigma? Why or why not?"


library(rethinking)
data(rugged)
d <- rugged
d$log_gdp <- log(d$rgdppc_2000)
dd <- d[ complete.cases(d$rgdppc_2000) , ]
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )


#When using ulam (HMC):
#Preprocc all variables before
#Trim down data to only the variable of interest

#R code 9.13 #list instead of Df because elements in a list can be of differetn lenghts
dat_slim <- list(
  log_gdp_std = dd$log_gdp_std,
  rugged_std = dd$rugged_std,
  cid = as.integer( dd$cid )
)

m1 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dunif( 0, 1 )
  ) , data=dat_slim , chains=4 , cores=4 )

precis(m9.1, depth = 2)
precis(m1, depth = 2)

#There is no difference in the posterior of sigma. The are enough observations 
#to overthrow the flatter prior of sigma.

#M2
# Modify the terrain ruggedness model again. This time, change the prior for b[cid]
# to dexp(0.3). What does this do to the posterior distribution? Can you explain it?

m2 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dexp( 0.3 ) ,
    sigma ~ dunif( 0, 1 )
  ) , data=dat_slim , chains=4 , cores=4 )

precis(m1, depth = 2)
precis(m2, depth = 2)
#This forces the slope to be positive despite the negative correlation between
#ruggedness and non-african contries gdp. Insted the slope becomes close to zero 


#M4

#R code 9.16
m4 <- ulam(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dat_slim , chains=4 , cores=4, warmup = 200 )

traceplot(m4, chains = 4)
precis(m9.1, depth = 2)
precis(m4, depth = 2)

# The chains converges pretty quickly and there is not a need for an extensive 
# warmup period in this case. With fewer iterations for the warmup the chains
# get an increased number of effective number of paratemers.


# H1
# Run the model below and then inspect the posterior distribution and explain what
# it is accomplishing.
mp <- ulam(
  alist(
    a ~ dnorm(0,1),
    b ~ dcauchy(0,1)
  ), data=list(y=1) , chains=1 )

precis(mp)
traceplot(mp)
#The model is unidentifiable. The cuachy distribution used
#to estimate the slope have undefinable expected value and variance.
# This is shown by the noncentered sampling in the trace plots.

## H2

library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce

# standardize variables
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )
d$A <- standardize( d$MedianAgeMarriage )

d_slim <- list(D= d$D, M=d$M, A= d$A)
#R code 5.3
m5.1 <- ulam(
  alist(
    D ~ dnorm( mu , sigma ),
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )) , data = d_slim, log_lik = TRUE )


#R code 5.6
# Model of marriage rate and divorce
m5.2 <- ulam(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d_slim, log_lik = TRUE )

m5.3 <- ulam(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d_slim , log_lik = T)

compare(m5.1, m5.2, m5.3, func = WAIC)

precis(m5.1)
precis(m5.3)


h2 = "The first model which include only age of marige as a predictor has the 
      best predictive power for out of sample. It seems to be reliably better 
      than the other models. The younger the couple is when getting married,
      the more likely they are to end up divorced. 
      "

##### H3
#The data 
#R code 6.2
N <- 100# number of individuals
set.seed(909)
height <- rnorm(N,10,2)# sim total height of each
leg_prop <- runif(N,0.4,0.5)# leg as proportion of height
leg_left <- leg_prop*height +# sim left leg as proportion + error
  rnorm( N , 0 , 0.02 )
leg_right <- leg_prop*height +# sim right leg as proportion + error
  rnorm( N , 0 , 0.02 )
# combine into data frame
d <- data.frame(height,leg_left,leg_right)



# Model 1
#R code 9.29
m5.8s <- ulam(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d, chains=4,
  start=list(a=10,bl=0,br=0.1,sigma=1) )

#Model 2 with contrained prior for right leg

 # R code 9.30
m5.8s2 <- ulam(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d, chains=4, cores = 4,
  constraints=list(br="lower=0"),
  start=list(a=10,bl=0,br=0.1,sigma=1) )


precis(m5.8s)
precis(m5.8s2)

H3 = "When the influence of the right leg parameter is contrained to be positive
      it also affect the posterior of the left leg parameter which becomes negative
      to compensate for this. This is because the two parameters are highly 
      correlated."

### H4

m5.8s <- ulam(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d, chains=4, cores = 4,
  start=list(a=10,bl=0,br=0.1,sigma=1), log_lik = T )

#Model 2 with contrained prior for right leg

# R code 9.30
m5.8s2 <- ulam(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d, chains=4, cores = 4,
  constraints=list(br="lower=0"),
  start=list(a=10,bl=0,br=0.1,sigma=1), log_lik = T )

compare(m5.8s, m5.8s2)
precis(m5.8s)
precis(m5.8s2)

H4 = "The first unconstrained model has the more effective number of 
      parameters for the betaparameters. This because when contrained, it shrinks
      the possible paramter space to match up with the posterior."

### H5
#The metropolis algorithm
#R code 9.1 #The king Markov example
num_weeks <- 1e5
positions <- rep(0,num_weeks)
current <- 10
for ( i in 1:num_weeks ) {
  ## record current position
  positions[i] <- current
  ## flip coin to generate proposal
  proposal <- current + current*0.76
  ## now make sure he loops around the archipelago
  if ( proposal < 1 ) proposal <- 10
  if ( proposal > 10 ) proposal <- 1
  ## move?
  prob_move <- proposal/current
  current <- ifelse( runif(1) < prob_move , proposal , current )
}

par(mfrow = c(1,1))
#R code 9.2
plot( 1:100 , positions[1:100] )

#R code 9.3
plot( table( positions ) )


#H6

#R code 2.8
n_samples <- 1000
tosses <- rep( NA , n_samples )
current <- 0.5
W <- 6
L <- 3
for ( i in 1:n_samples ) {
  tosses[i] <- current
  proposal <- rnorm( 1 , current , 0.1 )
  if ( proposal < 0 ) proposal <- abs( proposal )
  if ( proposal > 1 ) proposal <- 2 - proposal
  q0 <- dbinom( W , W+L , current )
  q1 <- dbinom( W , W+L , proposal )
  prob_move <- q1/q0
  current <- ifelse( runif(1) < prob_move , proposal , current )
}
dens( tosses , xlim=c(0,1) )
curve( dbeta( x , W+1 , L+1 ) , lty=2 , add=TRUE )

#H7
# U needs to return neg-log-probability
U <- function( q) {
  U <- dbinom(w,n,q , log = T + dunif(q, log = T))
  return( -U )
}


#The derivative of of the logarithm of any univariate Gaussian with mean a and 
#standard deviation b with respect to a is:
# ∂logN(y|a,b) = y−a/b2
#     ∂a

#R code 9.6

# gradient function
# need vector of partial derivatives of U with respect to vector q
U_gradient <- function( q) {
  # calculate the derivative of the binomial log-probability
  G <- (w - n*q) / (q*(1 - q))
  
  return( -G )# negative bc energy is neg-log-prob
}

w = 6
n = 9

Q <- list()
Q$q = 0.5

n_samples = 1e3

epsilon = 0.3

L = 10

# initalize empty vectors to store data
samples <- rep(NA, n_samples)
accepted <- rep(NA, n_samples)
HMC2 <- function (U, grad_U, epsilon, L, current_q) {
  q = current_q
  p = rnorm(length(q),0,1) # random flick - p is momentum.
  current_p = p
  # Make a half step for momentum at the beginning
  p = p - epsilon * grad_U(q) / 2
  # initialize bookkeeping - saves trajectory
  qtraj <- matrix(NA,nrow=L+1,ncol=length(q))
  ptraj <- qtraj
  qtraj[1,] <- current_q
  ptraj[1,] <- p
  
  #R code 9.9
  # Alternate full steps for position and momentum
  for ( i in 1:L ) {
    q = q + epsilon * p # Full step for the position
    # Make a full step for the momentum, except at end of trajectory
    if ( i!=L ) {
      p = p - epsilon * grad_U(q)
      ptraj[i+1,] <- p
    }
    qtraj[i+1,] <- q
  }
  
  #R code 9.10
  # Make a half step for momentum at the end
  p = p - epsilon * grad_U(q) / 2
  ptraj[L+1,] <- p
  # Negate momentum at end of trajectory to make the proposal symmetric
  p = -p
  # Evaluate potential and kinetic energies at start and end of trajectory
  current_U = U(current_q)
  current_K = sum(current_p^2) / 2
  proposed_U = U(q)
  proposed_K = sum(p^2) / 2
  # Accept or reject the state at end of trajectory, returning either
  # the position at the end of the trajectory or the initial position
  accept <- 0
  if (runif(1) < exp(current_U-proposed_U+current_K-proposed_K)) {
    new_q <- q # accept
    accept <- 1
  } else new_q <- current_q # reject
  return(list( q=new_q, traj=qtraj, ptraj=ptraj, accept=accept ))
}

for (i in 1:n_samples){
  Q = HMC2(U = U,
           grad_U = U_gradient,
           epsilon = epsilon,
           L = L, 
           current_q = Q$q)
  samples[i] = Q$q
  accepted[i] = Q$accepted
}

