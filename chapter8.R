#R code 8.1
library(rethinking)
data(rugged)
d <- rugged

# make log version of outcome
d$log_gdp <- log( d$rgdppc_2000 )

# extract countries with GDP data
dd <- d[ complete.cases(d$rgdppc_2000) , ]

# rescale variables
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)

#R code 8.2
m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a + b*( rugged_std - 0.215 ) ,
    a ~ dnorm( 1 , 1 ) ,
    b ~ dnorm( 0 , 1 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )


# Looking at prior predictions to evaluate priors
#R code 8.3
set.seed(7)
prior <- extract.prior( m8.1 )

# set up the plot dimensions
plot( NULL , xlim=c(0,1) , ylim=c(0.5,1.5) ,
      xlab="ruggedness" , ylab="“log GDP” ")
abline( h=min(dd$log_gdp_std) , lty=2 )
abline( h=max(dd$log_gdp_std) , lty=2 )

# draw 50 lines from the prior
rugged_seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
mu <- link( m8.1 , post=prior , data=data.frame(rugged_std=rugged_seq) )
for ( i in 1:50 ) lines( rugged_seq , mu[i,] , col=col.alpha("black",0.3) )


#The prior of a is too vague and many of the lines' slope are to steap

#More than half of the slopes are larger than .6
#R code 8.4 
sum( abs(prior$b) > 0.6 ) / length(prior$b)

#New model with updated priors
#R code 8.5
m8.1 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a + b*( rugged_std - 0.215 ) ,
    a ~ dnorm( 1 , 0.1 ) ,
    b ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp(1)
  ) , data=dd )

#R code 8.6
precis( m8.1 )


#Including information about African and non-african continent in the model
#R code 8.7
# make variable to index Africa (1) or not (2)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

#R code 8.8
m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )

#R code 8.9
compare( m8.1 , m8.2 )

#R code 8.10
precis( m8.2 , depth=2 )

# The difference between the two intercepts (when ruggedness is the same as the sample mean)
#R code 8.7
post <- extract.samples(m8.2)
diff_a1_a2 <- post$a[,1] - post$a[,2]
PI( diff_a1_a2 )


#R code 8.12
rugged.seq <- seq( from=-0.1 , to=1.1 , length.out=30 )
# compute mu over samples, fixing cid=2 and then cid=1
mu.NotAfrica <- link( m8.2 ,
                      data=data.frame( cid=2 , rugged_std=rugged.seq ) )
mu.Africa <- link( m8.2 ,
                   data=data.frame( cid=1 , rugged_std=rugged.seq ) )
# summarize to means and intervals
mu.NotAfrica_mu <- apply( mu.NotAfrica , 2 , mean )
mu.NotAfrica_ci <- apply( mu.NotAfrica , 2 , PI , prob=0.97 )
mu.Africa_mu <- apply( mu.Africa , 2 , mean )
mu.Africa_ci <- apply( mu.Africa , 2 , PI , prob=0.97 )

#Add prediction lines to the plot
lines(rugged.seq, mu.NotAfrica_mu, col="blue")
shade(mu.NotAfrica_ci, rugged.seq, col = col.alpha( "blue", .15))
lines(rugged.seq, mu.Africa_mu, col="red")
shade(mu.Africa_ci, rugged.seq, col = col.alpha("red", .15))

#Indexing the slope as well to see if the interaction effect is captured
#R code 8.13
m8.3 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b[cid]*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b[cid] ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )

# R code  8.14
precis(m8.3, depth = 2)

#R code 8.15
compare( m8.1 , m8.2 , m8.3 , func=PSIS )

# The third model might be a little overfit, the standard deviation of the difference
# between the senond and third model is as high as the difference itself

#Looking at PSIS values for individual points
#R code 8.16
plot( PSIS( m8.3 , pointwise=TRUE )$k )

#There are some influential outliers, robust regression using the student t
# distribution could be a good idea

#Plotting the slopes of African and non-african GDP
#R code 8.17
# plot Africa - cid=1
d.A1 <- dd[ dd$cid==1 , ]
plot( d.A1$rugged_std , d.A1$log_gdp_std , pch=16 , col=rangi2 ,
      xlab="ruggedness (standardized)" , ylab="log GDP (as proportion of mean)" ,
      xlim=c(0,1) )
mu <- link( m8.3 , data=data.frame( cid=1 , rugged_std=rugged_seq ) )
mu_mean <- apply( mu , 2 , mean )
mu_ci <- apply( mu , 2 , PI , prob=0.97 )
lines( rugged_seq , mu_mean , lwd=2 )
shade( mu_ci , rugged_seq , col=col.alpha(rangi2,0.3) )
mtext("African nations")
# plot non-Africa - cid=2
d.A0 <- dd[ dd$cid==2 , ]
plot( d.A0$rugged_std , d.A0$log_gdp_std , pch=1 , col="black" ,
      xlab="ruggedness (standardized)", ylab="log GDP (as proportion of mean)" ,
      xlim=c(0,1) )
mu <- link( m8.3 , data=data.frame( cid=2 , rugged_std=rugged_seq ) )
mu_mean <- apply( mu , 2 , mean )
mu_ci <- apply( mu , 2 , PI , prob=0.97 )
lines( rugged_seq , mu_mean , lwd=2 )
shade( mu_ci , rugged_seq )
mtext("Non-African nations")


#### SYMETRY OF INTERACTION ##########

#The interaction can be stated as:
# 
# (1)  How much does the association between ruggedness and log GDP depend upon whether the nation is in Africa?
# (2)  How much does the association of Africa with log GDP depend upon ruggedness?

#To a statistical model, these rephrasings are the same
#To plot the second phrasing, the inverted reverse interpretation
#R code 8.18
rugged_seq <- seq(from=-0.2,to=1.2,length.out=30)
muA <- link( m8.3 , data=data.frame(cid=1,rugged_std=rugged_seq) )
muN <- link( m8.3 , data=data.frame(cid=2,rugged_std=rugged_seq) )
delta <- muA - muN

mu_delta = apply(delta, 2, mean)
pi_delta = apply(delta, 2, PI)


#Plto the differnce
plot(rugged_seq, mu_delta, type ="n",
     ylim = c(-0.3, 0.2))

lines(rugged_seq, mu_delta)
shade(pi_delta, rugged_seq,col=col.alpha("black", .15))
abline(h = 0, lty = "dashed")

##### CONTINUOUS INTERACTIONS ########

#R code 8.19
library(rethinking)
data(tulips)
d <- tulips
str(d)

#Scaling variables by their max
#R code 8.20
d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

#R code 8.21
#Too wide standard deviation
a <- rnorm( 1e4 , 0.5 , 1 ); sum( a < 0 | a > 1 ) / length( a )

#R code 8.22
a <- rnorm( 1e4 , 0.5 , 0.25 ); sum( a < 0 | a > 1 ) / length( a )

#R code 8.2
m8.4 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bw*water_cent + bs*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dnorm( 0 , 0.25 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )

#R code 8.24
#Model with interaction
m8.5 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dnorm( 0 , 0.25 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )

precis(m8.4)
precis(m8.5)

#Plot posterior predictions using triptych (Three plots for three values)
#R code 8.25
par(mfrow=c(2,3)) # 3 plots in 1 row
models = c(m8.4, m8.5)
for (model in models){
  for ( s in -1:1 ) {
    idx <- which( d$shade_cent==s )
    plot( d$water_cent[idx] , d$blooms_std[idx] , xlim=c(-1,1) , ylim=c(0,1) ,
          xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
    mu <- link( model , data=data.frame( shade_cent=s , water_cent=-1:1 ) )
    for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
  }
}

#Plot piror predictions using triptych
#R code 8.26
set.seed(7)
prior <- extract.prior(m8.5)
par(mfrow=c(2,3)) # 3 plots in 1 row
for (model in models){
  for ( s in -1:1 ) {
    idx <- which( d$shade_cent==s )
    plot( d$water_cent[idx] , d$blooms_std[idx] , type = "n", 
          xlim=c(-1,1) , ylim=c(-0.1,1.1) ,
          xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
    mu <- link( model , post = prior, data=data.frame( shade_cent=s , water_cent=-1:1 ) )
    for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",i/20) )
    mean_mu = apply(mu, 2, max)
    lines(-1:1, mean_mu, lwd = 2)
    abline(h = c(-0, 1), lty = 'dashed' )
  }
}

######### PRATICE PROBLEMS #############
e1 = " 1: temperature, 2: Sex, 3: Functional tank"

e2 = "(1)  Caramelizing onions requires cooking over low heat and making sure the onions do not dry out."

#M1
M1 = "There is an interaction where the temperature 
      restric the photosyntetic process from happening,
      either by imparing the effect of water, light
      or both."
#M2
M2 = "mu = a + Bt(Bw + BS + Bws)"

par(mfrow=c(1,1))


#M3
wolves = rnorm(1e4, 5, 0.2)
pray = rnorm(1e4,10,5) * wolves ^ 2
seeds = rnorm(1e4,10,1) - pray/2
mu = pray + seeds
ravens = rnorm(1e4, mu, 2)

plot(wolves, ravens)
cor(wolves, ravens)

#M4
# Repeat the tulips analysis, but this time use priors that constrain
# the effect of water to be positive and the effect of shade to be negative. 
# Use prior predictive simulation. What do these prior assumptions mean for the 
# interaction prior, if anything


library(rethinking)
data(tulips)
d <- tulips
str(d)

#Scaling variables by their max
#R code 8.20
d$blooms_std <- d$blooms / max(d$blooms)
d$water_cent <- d$water - mean(d$water)
d$shade_cent <- d$shade - mean(d$shade)

#Model with interaction
M4 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dlnorm(0.5 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )


#Plot piror predictions using triptych
#R code 8.26
set.seed(7)
par(mfrow=c(2,3)) # 3 plots in 1 row
models = c(m8.5, M4)
for (model in models){
  prior <- extract.prior(model)
  for ( s in -1:1 ) {
    idx <- which( d$shade_cent==s )
    plot( d$water_cent[idx] , d$blooms_std[idx] , type = "n", 
          xlim=c(-1,1) , ylim=c(-0.1,1.1) ,
          xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
    mu <- link(model , post = prior, data=data.frame( shade_cent=s , water_cent=-1:1 ) )
    for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
    abline(h = c(-0, 1), lty = 'dashed' )
  }
  
}

par(mfrow=c(2,3)) # 3 plots in 1 row
for (model in models){
  for ( s in -1:1 ) {
    idx <- which( d$shade_cent==s )
    plot( d$water_cent[idx] , d$blooms_std[idx] , xlim=c(-1,1) , ylim=c(0,1) ,
          xlab="water" , ylab="blooms" , pch=16 , col=rangi2 )
    mu <- link( model , data=data.frame( shade_cent=s , water_cent=-1:1 ) )
    for ( i in 1:20 ) lines( -1:1 , mu[i,] , col=col.alpha("black",0.3) )
  }
}

#When the prior for waters effect is forced to be positive all the piror slopes 
# are positive. However it does not seems to affect the posterior that much.
# The large amount of data makes the choice of prior less crucial.


##H1
#Return to the data(tulips) example in the chapter. Now include the bed variable 
#as a predictor in the interaction model. 


d$ind = as.integer(d$bed)


H1 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a + bb[ind]* ind + bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bb[ind] ~dnorm(0, 0.25),
    bw ~ dlnorm(0.5 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )


precis(H1, depth = 3)


# H2
H2 <- quap(
  alist(
    blooms_std ~ dnorm( mu , sigma ) ,
    mu <- a+ bw*water_cent + bs*shade_cent + bws*water_cent*shade_cent ,
    a ~ dnorm( 0.5 , 0.25 ) ,
    bw ~ dlnorm(0.5 ) ,
    bs ~ dnorm( 0 , 0.25 ) ,
    bws ~ dnorm( 0 , 0.25 ) ,
    sigma ~ dexp( 1 )
  ) , data=d )

compare(H1, H2, WAIC = T)
compare(H1, H2, func = PSIS)

par(mfrow=c(1,1))
plot(compare(H1, H2, WAIC = T))


#According to the posteriors, bed does not seem to be a good predictor of water
#when the predictors are included in the model. Furthermore, the WAIC score
# does not give reliable more weight to the model that includes bed. 

#H3
#R code 8.1
library(rethinking)
data(rugged)
d <- rugged

# make log version of outcome
d$log_gdp <- log( d$rgdppc_2000 )

# extract countries with GDP data
dd <- d[ complete.cases(d$rgdppc_2000) , ]

# rescale variables
dd$log_gdp_std <- dd$log_gdp / mean(dd$log_gdp)
dd$rugged_std <- dd$rugged / max(dd$rugged)



#Including information about African and non-african continent in the model
#R code 8.7
# make variable to index Africa (1) or not (2)
dd$cid <- ifelse( dd$cont_africa==1 , 1 , 2 )

#R code 8.8
m8.2 <- quap(
  alist(
    log_gdp_std ~ dnorm( mu , sigma ) ,
    mu <- a[cid] + b*( rugged_std - 0.215 ) ,
    a[cid] ~ dnorm( 1 , 0.1 ) ,
    b ~ dnorm( 0 , 0.3 ) ,
    sigma ~ dexp( 1 )
  ) , data=dd )

max(WAIC(m8.2, pointwise = T)$WAIC)
