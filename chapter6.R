
# rcode 6.1
set.seed(1914)
N <- 200 # num grant proposals
p <- 0.1 # proportion to select
# uncorrelated newsworthiness and trustworthiness
nw <- rnorm(N)
tw <- rnorm(N)
# select top 10% of combined scores
s <- nw + tw # total score
q <- quantile( s , 1-p ) # top 10% threshold
selected <- ifelse( s >= q , TRUE , FALSE )
cor( tw[selected] , nw[selected] )
plot(nw,tw)

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

library(rethinking)

#R code 6.3
m6.1 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left + br*leg_right ,
    a ~ dnorm( 10 , 100 ),
    bl ~ dnorm( 2 , 10 ) ,
    br ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
precis(m6.1)

#R code 6.4
plot(precis(m6.1))

#R code 6.5
post <- extract.samples(m6.1)
plot( bl ~ br , post , col=col.alpha(rangi2,0.1) , pch=16 )

#R code 6.6
sum_blbr <- post$bl + post$br
dens( sum_blbr , col=rangi2 , lwd=2 , xlab="sum of bl and br" )

#R code 6.7
m6.2 <- quap(
  alist(
    height ~ dnorm( mu , sigma ) ,
    mu <- a + bl*leg_left,
    a ~ dnorm( 10 , 100 ) ,
    bl ~ dnorm( 2 , 10 ) ,
    sigma ~ dexp( 1 )
  ) , data=d 
)  
precis(m6.2)  


#R code 6.8
library(rethinking)
data(milk)
d <- milk
d$K <- standardize( d$kcal.per.g )
d$F <- standardize( d$perc.fat )
d$L <- standardize( d$perc.lactose )


#R code 6.9
# kcal.per.g regressed on perc.fat
m6.3 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bF*F ,
    a ~ dnorm( 0 , 0.2 ) ,
    bF ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
# kcal.per.g regressed on perc.lactose
m6.4 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bL*L ,
    a ~ dnorm( 0 , 0.2 ) ,
    bL ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )
precis( m6.3 )
precis( m6.4 )

#R code 6.10
m6.5 <- quap(
  alist(
    K ~ dnorm( mu , sigma ) ,
    mu <- a + bF*F + bL*L ,
    a ~ dnorm( 0 , 0.2 ) ,
    bF ~ dnorm( 0 , 0.5 ) ,
    bL ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) ,
  data=d )
precis( m6.5 )

#R code 6.11
pairs( ~ kcal.per.g + perc.fat + perc.lactose , data=d , col=rangi2 )


#The code below makes a function that generates correlated predictors, 
#fits a model, and returns the standard deviation of the posterior
#distribution for the slope relating perc.fat to kcal.per.g. 
#Then the code repeatedly calls this function, with different degrees 
#of correlation as input, and collects the results.

#R code 6.12
library(rethinking)
data(milk)
d <- milk
sim.coll <- function( r=0.9 ) {
  d$x <- rnorm( nrow(d) , mean=r*d$perc.fat ,
                sd=sqrt( (1-r^2)*var(d$perc.fat) ) )
  m <- lm( kcal.per.g ~ perc.fat + x , data=d )
  sqrt( diag( vcov(m) ) )[2] # stddev of parameter
}
rep.sim.coll <- function( r=0.9 , n=100 ) {
  stddev <- replicate( n , sim.coll(r) )
  mean(stddev)
}
r.seq <- seq(from=0,to=0.99,by=0.01)
stddev <- sapply( r.seq , function(z) rep.sim.coll(r=z,n=100) )
plot( stddev ~ r.seq , type="l" , col=rangi2, lwd=2 , xlab="correlation" )


########### Post-treatment bias ######################

#R code 6.13
set.seed(71)
# number of plants
N <- 100

# simulate initial heights
h0 <- rnorm(N,10,2)

# assign treatments and simulate fungus and growth
treatment <- rep( 0:1 , each=N/2 )
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 )
h1 <- h0 + rnorm(N, 5 - 3*fungus)

# compose a clean data frame
d <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )
precis(d)

#R code 6.14
sim_p <- rlnorm( 1e4 , 0 , 0.25 )
precis( data.frame(sim_p) )

#R code 6.15
m6.6 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0*p,
    p ~ dlnorm( 0 , 0.25 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.6)

#R code 6.16
m6.7 <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm( 0 , 0.2 ) ,
    bt ~ dnorm( 0 , 0.5 ),
    bf ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.7)


#R code 6.17, without fungus as predictor
m6.8 <- quap(
        alist(
          h1 ~ dnorm( mu , sigma ),
          mu <- h0 * p,
          p <- a + bt*treatment,
          a ~ dlnorm( 0 , 0.2 ),
          bt ~ dnorm( 0 , 0.5 ),
          sigma ~ dexp( 1 )
        ), data=d )
precis(m6.8)


#R code 6.18
library(dagitty)
plant_dag <- dagitty( "dag {
  H_0 -> H_1
  F -> H_1
  T -> F
}")
coordinates( plant_dag ) <- list( x=c(H_0=0,T=2,F=1.5,H_1=1) ,
                                  y=c(H_0=0,T=0,F=0,H_1=0) )
drawdag( plant_dag )


#R code 6.19
impliedConditionalIndependencies(plant_dag)


#R code 6.20, new data where fungus does not have an effect but moisture 
# has an effect on both fungus and growth
set.seed(71)
N <- 1000
h0 <- rnorm(N,10,2)
treatment <- rep( 0:1 , each=N/2 )
M <- rbern(N)
fungus <- rbinom( N , size=1 , prob=0.5 - treatment*0.4 + 0.4*M )
h1 <- h0 + rnorm( N , 5 + 3*M )
d2 <- data.frame( h0=h0 , h1=h1 , treatment=treatment , fungus=fungus )


#Re-running the models with the new data

m6.7b <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment + bf*fungus,
    a ~ dlnorm( 0 , 0.2 ) ,
    bt ~ dnorm( 0 , 0.5 ),
    bf ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data=d2 )

m6.8b <- quap(
  alist(
    h1 ~ dnorm( mu , sigma ),
    mu <- h0 * p,
    p <- a + bt*treatment,
    a ~ dlnorm( 0 , 0.2 ),
    bt ~ dnorm( 0 , 0.5 ),
    sigma ~ dexp( 1 )
  ), data=d2 )

precis(m6.7b)
precis(m6.8b)


############3 Colider bias #######################3

#Simulating data for age, happiness and marriage where
# happiness is static and influences marriage rate which is also influenced 
# by age
#R code 6.21
library(rethinking)
d <- sim_happiness( seed=1977 , N_years=1000 )
precis(d)

#rescale age
#R code 6.22
d2 <- d[ d$age>17 , ] # only adults
d2$A <- ( d2$age - 18 ) / ( 65 - 18 )

#R code 6.23
d2$mid <- d2$married + 1 # New index for marriage 1 - single, 2 -married

m6.9 <- quap(
  alist(
    happiness ~ dnorm( mu , sigma ),
    mu <- a[mid] + bA*A,
    a[mid] ~ dnorm( 0 , 1 ) ,
    bA ~ dnorm( 0 , 2 ) ,
    sigma ~ dexp(1)
  ) , data=d2 )
precis(m6.9,depth=2)

#The model below omits marriage status to avoid collider bias
#R code 6.24
m6.10 <- quap(
  alist(
    happiness ~ dnorm( mu , sigma ),
    mu <- a + bA*A,
    a ~ dnorm( 0 , 1 ),
    bA ~ dnorm( 0 , 2 ),
    sigma ~ dexp(1)
  ) , data=d2 )
precis(m6.10)


#R code 6.25
#These parameters are like slopes in a regression model
N <- 200 # number of grandparent-parent-child triads
b_GP <- 1 # direct effect of G on P
b_GC <- 0 # direct effect of G on C
b_PC <- 1 # direct effect of P on C
b_U <- 2 # direct effect of U on P and C

#Now we use these slopes to draw random observations:
#R code 6.26
set.seed(1)
U <- 2*rbern( N , 0.5 ) - 1
G <- rnorm( N )
P <- rnorm( N , b_GP*G + b_U*U )
C <- rnorm( N , b_PC*P + b_GC*G + b_U*U )
d <- data.frame( C=C , P=P , G=G , U=U )


#R code 6.27, checking Geffect on c while controlling P 
m6.11 <- quap(
  alist(
    C ~ dnorm( mu , sigma ),
    mu <- a + b_PC*P + b_GC*G,
    a ~ dnorm( 0 , 1 ),
    c(b_PC,b_GC) ~ dnorm( 0 , 1 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.11)

#R code 6.28, including the unobserved variable
m6.12 <- quap(
  alist(
    C ~ dnorm( mu , sigma ),
    mu <- a + b_PC*P + b_GC*G + b_U*U,
    a ~ dnorm( 0 , 1 ) ,
    c(b_PC,b_GC,b_U) ~ dnorm( 0 , 1 ),
    sigma ~ dexp( 1 )
  ), data=d )
precis(m6.12)

#R code 6.29
library(dagitty)
dag_6.1 <- dagitty( "dag {
  U [unobserved]
  X -> Y
  X <- U <- A -> C -> Y
  U -> B <- C
}")
adjustmentSets( dag_6.1 , exposure="X" , outcome="Y" )

#R code 6.30
library(dagitty)
dag_6.2 <- dagitty( "dag {
  A -> D
  A -> M -> D
  A <- S -> M
  S -> W -> D
}")

adjustmentSets( dag_6.2 , exposure="W", outcome="D" )

#R code 6.31
impliedConditionalIndependencies( dag_6.2 )


################## Practice problems ##################

e1 = "multicolinearity, post-treatment bias and collider bias"

e2= "Predicting langugae learning performance based on a battery of cognitve
    tests, such as n-back, stroop, and ravens. Performance on these tests are
    likely to be correlated. Using all of them will affect the predictive power
    of an individual test on language learning. If there is a negative correlation
    it can even remove the one predictors influence of the posterior distribution."

e3 = "Fork, Pipe, Collier, Descendant
      Coditioning of a fork or a pipe will close the relationship betweem two varibles 
      connected by the confounder. Conditioning on a collider will remove the 
      independence, it will open up a closed path. The effect on conditioning on
      on a descendant is determined by its parent. Conditioning on the descendant
      will partly give information about the parent"

#m1
library(dagitty)
m1 <- dagitty("dag{ 
  U [unobserved]
  V [unobserved]
  X -> Y
  X <- U <- A -> C -> Y
  U -> B <- C
  C <- V -> Y}")

adjustmentSets(m2, exposure = "X", outcome = "Y")
#This does not open any new paths because it is a collider, no further action 
#necessary.

#M2
X <- rnorm(1000,0,1)
Z <- X + rnorm(1000, 0, 0.3)
Y <- rnorm(1000,Z, 0.5)
cor(X,Z)
cor(X,Y)
cor(Z,Y)
d.m2 = data.frame(X=X, Z=Z, Y=Y)

library(rethinking)
m2 <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bx*X + bz*Z,
    a ~ dnorm(0,0.2),
    bx ~dnorm(0,0.4),
    bz ~dnorm(0,0.4),
    sigma ~ exp(1)
  ), data = d.m2
)
#Model with just X as predictor
m2.b <- quap(
  alist(
    Y ~ dnorm(mu, sigma),
    mu <- a + bx*X ,
    a ~ dnorm(0,0.2),
    bx ~dnorm(0,0.4),
    sigma ~ exp(1)
  ), data = d.m2
)
precis(m2)
precis(m2.b)
#There is no multicollinearity, Z act as a mediator. Not including Z in the model
#would lead to a spurious association between X and Y. In contrast to the left
# right leg example, the outcome variable is only caused by one of the predictors.

#M3
dag1 = dagitty("dag{
               x -> y
               x <- z -> y 
               z <- a -> y}")

dag2= dagitty("dag{
               x -> y
               x -> z -> y 
               z <- a -> y}")

dag3 = dagitty("dag{
               x -> y
               x -> z <- y 
               z <- a -> z}")

dag4 = dagitty("dag{
               x -> y
               x -> z -> y 
               z <- a -> z}")

#My answers
#dag 1 you could condition on a
#dag 2 you should condition on z
#dag 3 nothing
#dag 4 condition on z

####H1
library(rethinking)
data(WaffleDivorce)
d = WaffleDivorce

d$A = standardize(d$MedianAgeMarriage)
d$M = standardize(d$Marriage)
d$D = standardize(d$Divorce)
d$W = standardize(d$WaffleHouses)
d$S = standardize(d$South)


library(dagitty)
dag_H1 <- dagitty("dag{
                  S -> W -> D
                  A <- S -> M
                  A -> D
                  A -> M -> D
                  }") 

adjustmentSets(dag_e1, exposure = "W", outcome = "D")

H1 = quap(
    alist(
      D ~ dnorm(mu, sigma),
      mu <- a + bS*S + bW*W ,
      a ~ dnorm(0, 0.1),
      bS ~ dnorm(0, 0.5),
      bW ~ dnorm(0, 0.5),
      sigma ~ dexp(1)
    ), data = d
)
plot(precis(e1))

### H2

impliedConditionalIndependencies(dag_H1)


#A _||_ W | S
H2a = quap(
  alist(
    A ~ dnorm(mu, sigma),
    mu <- a + bS*S + bW*W ,
    a ~ dnorm(0, 0.1),
    bS ~ dnorm(0, 0.5),
    bW ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
plot(precis(H2a))

#D _||_ S | A, M, W
H2b = quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + bS*S + bW*W + bA*A + bM*M ,
    a ~ dnorm(0, 0.1),
    c(bS, bW, bA, bM) ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
plot(precis(H2b))

#M _||_ W | S
H2c = quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bS*S + bW*W ,
    a ~ dnorm(0, 0.1),
    bS ~ dnorm(0, 0.5),
    bW ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d
)
plot(precis(H2c))




##### Home ẃork week 3 ###
#https://github.com/rmcelreath/statrethinking_winter2019/blob/master/homework/week03.pdf

#data
library(rethinking)
data(foxes)


#Dag for all problems
library(dagitty)
foxes_dag <- dagitty("dag{
                     F -> G -> W
                     A -> F -> W
}")

#standardise variables
foxes$G = standardize(foxes$groupsize)
foxes$W = standardize(foxes$weight)
foxes$F = standardize(foxes$avgfood)
foxes$A = standardize(foxes$area)

with(foxes, hist(G))
with(foxes, hist(W))
with(foxes, hist(F))
with(foxes, hist(A))

#1. Model to infer the total causal influence of area on weights. Would increasing
# the area for each fox make them heavier (positive correlation)

adjustmentSets(foxes_dag, exposure = "A", outcome = "W")

#Only front door paths, no need for control



#Prior predictive simulations
sim_a <- rnorm(1e3, 0, 0.6)
sim_w <- rnorm(1e3, 0, 0.8)
dens(sim_a)
dens(foxes$A)
dens(sim_w)
dens(foxes$W)

# Model 
area_model <- quap(
  alist(
    W ~ dnorm(mu, sigma),
    mu <- a + bA*A,
    a ~ dnorm(0, 0.2),
    bA  ~ dnorm(0, 0.6),
    sigma ~ dexp(1)
  ), data = foxes
)

precis(area_model)


## 2. Infer the causal impact of food, which covariates needs to be adjust?

adjustmentSets(foxes_dag, exposure = "F", outcome = "W")

sim_food = runif(1000, -2,2)

dens(foxes$F)
dens(standardize(sim_food))
food_model <- quap(
  alist(
    W ~ dnorm(mu, sigma),
    mu <- a + bF*F ,
    a ~ dnorm(0, 0.2),
    bF  ~ dnorm(0, 0.8),
    sigma ~ dexp(1)
  ), data = foxes
)

precis(food_model)


# 3 Infer causal impact of group size

adjustmentSets(foxes_dag, exposure = "G", outcome = "W")

dens(foxes$G)
sim_G = rnorm(1000, 0, 0.9)
dens(sim_G)
group_model <- quap(
  alist(
    W ~ dnorm(mu, sigma),
    mu <- a + bG *G + bF*F ,
    a ~ dnorm(0, 0.2),
    bF  ~ dnorm(0, 0.8),
    bG  ~ dnorm(0, 0.9),
    sigma ~ dexp(1)
  ), data = foxes
)

precis(group_model)
