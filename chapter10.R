#Maximum entropy and generalised linear models

#Throwing pebbles into buckets

#R code 10.1
p <- list()
p$A <- c(0,0,10,0,0)
p$B <- c(0,1,8,1,0)
p$C <- c(0,2,6,2,0)
p$D <- c(1,2,4,2,1)
p$E <- c(2,2,2,2,2)

#Normalise so they become probaility distributions

#R code 10.2
p_norm <- lapply( p , function(q) q/sum(q))

#Compute informational entropy of each
#R code 10.3
( H <- sapply( p_norm , function(q) -sum(ifelse(q==0,0,q*log(q))) ) )

# E has the highest number of ways to be realised and thereby the highest entropy

#R code 10.4
ways <- c(1,90,1260,37800,113400)
logwayspp <- log(ways)/10

plot(logwayspp, H, 'l', ylab = "Entropy")


#The distribution that can happen in greatest number of ways is
# also the most plausible distribution. 
#This can be called the Maximum entropy distribution


#Example with the blue and white marbles

#With an exceptec value of 1 blue (on two draws with replacement)
#R code 10.5
# build list of the candidate distributions
p <- list()
p[[1]] <- c(1/4,1/4,1/4,1/4)
p[[2]] <- c(2/6,1/6,1/6,2/6)
p[[3]] <- c(1/6,2/6,2/6,1/6)
p[[4]] <- c(1/8,4/8,2/8,1/8)

# compute expected value of each
sapply( p , function(p) sum(p*c(0,1,1,2)) )

#R code 10.6
# compute entropy of each distribution
sapply( p , function(p) -sum( p*log(p) ) )


#Excpected value of 1.4 blue (70% blue, 30% white)
#R code 10.7
p <- 0.7
( A <- c( (1-p)^2 , p*(1-p) , (1-p)*p , p^2 ) )

#Entropy of the distribution above
#R code 10.8
-sum( A*log(A) )

#Function to simulate thousands distribution with expected value of 1.4
# that return the entropy of that distrubution
#R code 10.9
sim.p <- function(G=1.4) {
  x123 <- runif(3)
  x4 <- ( (G)*sum(x123)-x123[2]-x123[3] )/(2-G)
  z <- sum( c(x123,x4) )
  p <- c( x123 , x4 )/z
  list( H=-sum( p*log(p) ) , p=p )
}

#Calling that fucntion 100000 times and plot the entropy
#R code 10.10
H <- replicate( 1e5 , sim.p(1.4) )
rethinking::dens( as.numeric(H[1,]) , adj=0.1 )


#R code 10.11
entropies <- as.numeric(H[1,])
distributions <- H[2,]

#R code 10.12
max(entropies)

#R code 10.13
distributions[ which.max(entropies) ]
