#Integers counting


#Binomninal regression

#Logistic regression

#Chimpanzee prosociality experiment
#R code 11.1
library(rethinking)
data(chimpanzees)
d <- chimpanzees

#Index variable
#R code 11.2
d$treatment <- 1 + d$prosoc_left + 2*d$condition


#Validate the indexing
#R code 11.3
xtabs( ~ treatment + prosoc_left + condition , d )

#Quap with flat prior
#R code 11.4
m11.1 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a ,
    a ~ dnorm( 0 , 10 )
  ) , data=d )

#sample from prior
#R code 11.5
set.seed(1999)
prior <- extract.prior( m11.1 , n=1e4 )


#The parameters need to be converted to the output scale
#by using the inverse logit function
#R code 11.6
p <- inv_logit( prior$a )
dens( p , adj=0.1 )


#Flat prior for treatment effect
#R code 11.7
m11.2 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 10 )
  ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.2 , n=1e4 )
#Probility of pulling left for each treatment
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )

#R code 11.8
dens( abs( p[,1] - p[,2] ) , adj=0.1 )

#The average prior difference with informative priors
#R code 11.9
m11.3 <- quap(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a + b[treatment] ,
    a ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=d )
set.seed(1999)
prior <- extract.prior( m11.3 , n=1e4 )
p <- sapply( 1:4 , function(k) inv_logit( prior$a + prior$b[,k] ) )
mean( abs( p[,1] - p[,2] ) )

#Trimmed data list for the HMC
#R code 11.10
# trimmed data list
dat_list <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  treatment = as.integer(d$treatment) )

#R code 11.11
m11.4 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + b[treatment] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    b[treatment] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list , chains=4 , cores = 4, log_lik=TRUE )
precis( m11.4 , depth=2 )


#Each chimpazee's tendency to pull the left lever in output scale
#R code 11.12
post <- extract.samples(m11.4)
p_left <- inv_logit( post$a )
plot( precis( as.data.frame(p_left) ) , xlim=c(0,1) )

#Treatment effects on logit scale
#R code 11.13
labs <- c("R/N","L/N","R/P","L/P")
plot( precis( m11.4 , depth=2 , pars="b" ) , labels=labs )


#Calculate the difference between partner/no partner condition
#R code 11.14
diffs <- list(
  db13 = post$b[,1] - post$b[,3],
  db24 = post$b[,2] - post$b[,4] )
plot( precis(diffs) )


#Posterior prediction check

# calculate the proportion in each combination of actor and treatment
#R code 11.15
pl <- by(d$pulled_left, list(d$actor, d$treatment), mean)
pl[1,]

#Display the proportions against the predictive distributions
#R code 11.16
plot( NULL , xlim=c(1,28) , ylim=c(0,1) , xlab="" ,
      ylab="proportion left lever" , xaxt="n" , yaxt="n" )
axis( 2 , at=c(0,0.5,1) , labels=c(0,0.5,1) )
abline( h=0.5 , lty=2 )
for ( j in 1:7 ) abline( v=(j-1)*4+4.5 , lwd=0.5 )
for ( j in 1:7 ) text( (j-1)*4+2.5 , 1.1 , concat("actor",j) , xpd=TRUE )
for ( j in (1:7)[-2] ) {
  lines( (j-1)*4+c(1,3) , pl[j,c(1,3)] , lwd=2 , col=rangi2 )
  lines( (j-1)*4+c(2,4) , pl[j,c(2,4)] , lwd=2 , col=rangi2 )
}
points( 1:28 , t(pl) , pch=16 , col="white", cex=1.7 )
points( 1:28 , t(pl) , pch=c(1,1,16,16) , col=rangi2 , lwd=2 )
yoff <- 0.01
text( 1 , pl[1,1]-yoff , "R/N" , pos=1 , cex=0.8 )
text( 2 , pl[1,2]+yoff , "L/N" , pos=3 , cex=0.8 )
text( 3 , pl[1,3]-yoff , "R/P" , pos=1 , cex=0.8 )
text( 4 , pl[1,4]+yoff , "L/P" , pos=3 , cex=0.8 )
mtext( "observed proportions\n" )

#Posterior predictions using link
#R code 11.17
dat <- list( actor=rep(1:7,each=4) , treatment=rep(1:4,times=7) )
p_post <- link( m11.4 , data=dat )
p_mu <- apply( p_post , 2 , mean )
p_ci <- apply( p_post , 2 , PI )

#New model that includes interaction

#New index variables for interaction
#R code 11.18
d$side <- d$prosoc_left + 1 # right 1, left 2
d$cond <- d$condition + 1 # no partner 1, partner 2

#R code 11.19
dat_list2 <- list(
  pulled_left = d$pulled_left,
  actor = d$actor,
  side = d$side,
  cond = d$cond )

m11.5 <- ulam(
  alist(
    pulled_left ~ dbinom( 1 , p ) ,
    logit(p) <- a[actor] + bs[side] + bc[cond] ,
    a[actor] ~ dnorm( 0 , 1.5 ),
    bs[side] ~ dnorm( 0 , 0.5 ),
    bc[cond] ~ dnorm( 0 , 0.5 )
  ) , data=dat_list2 , chains=4 , cores= 4, log_lik=TRUE )

#Compare the two models
#R code 11.20
compare( m11.5 , m11.4 , func=PSIS )

plot(precis(m11.5, depth = 2))



#Why take pointwise log probability
#R code 11.21
post <- extract.samples( m11.4 , clean=FALSE )
str(post)

stancode(m11.4)

#Run the model without the ulam funciton

#R code 11.22
m11.4_stan_code <- stancode(m11.4)
m11.4_stan <- stan( model_code=m11.4_stan_code , data=dat_list , chains=4 )
compare( m11.4_stan , m11.4 )
