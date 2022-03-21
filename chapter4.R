#clear workspace
rm(list = ls())
#clear plots
graphics.off()
#clear the console
cat("\014") 

library(rethinking)

#R code 4.1
pos <- replicate( 1000 , sum( runif(16,-1,1) ) )
plot(density(pos))

sum(pos)
mean(pos)


#R code 4.2
prod( 1 + runif(12,0,0.1) )

#R code 4.3
growth <- replicate( 10000 , prod( 1 + runif(12,0,0.1) ) )
dens( growth , norm.comp=TRUE , col = "red")

#R code 4.4
big <- replicate( 10000 , prod( 1 + runif(12,0,0.5) ) )
small <- replicate( 10000 , prod( 1 + runif(12,0,0.01) ) )
dens(big, norm.comp = TRUE)
dens(small, norm.comp = TRUE, col="red")


#R code 4.5
log.big <- replicate( 10000 , log(prod(1 + runif(12,0,0.5))) )
dens(log.big, norm.comp = TRUE, col="red")

#Probability density can exceed 1 but the area under the curve 
#(the probability mass never exceed 1)
dnorm(0,0,0.1)

#R code 4.6
w <- 6; n <- 9;
p_grid <- seq(from=0,to=1,length.out=100)
posterior <- dbinom(w,n,p_grid)*dunif(p_grid,0,1)
posterior <- posterior/sum(posterior)

plot(p_grid,posterior)
