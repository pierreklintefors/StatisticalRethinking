rm(list = ls())
graphics.off()

#2.1
ways = c(0, 3, 8, 9, 0)

ways/sum(ways)


#2.2
dbinom(6, size = 9, prob = 0.5)
dbinom(6, size = 9, prob = 0.6)
dbinom(6, size = 9, prob = 0.7)

?dbinom()

#2.3

#define grid
p_grid = seq(from=0, to=1, length.out = 100)

#define prior
prior = rep(1,20)

#compute likelihood at each value in grid 
likelihood = dbinom(6, size = 9, prob = p_grid)

#compute the product of likelihood and prior
unstd_posterior = prior * likelihood

posterior = unstd_posterior/sum(unstd_posterior)

#2.4#Ploting the posterior

plot(p_grid, posterior, type = "b", 
     xlab = "Probability of water", ylab = "Posterior probability")


#2.5## Differnt priors
prior2 <- ifelse( p_grid < 0.5 , 0 , 1)

prior3 <- exp( -5*abs( p_grid - 0.5 ) )

#compute the product of likelihood and prior
unstd_posterior2 = prior2 * likelihood

posterior2 = unstd_posterior2/sum(unstd_posterior2)

#compute the product of likelihood and prior
unstd_posterior3 = prior3 * likelihood

posterior3 = unstd_posterior3/sum(unstd_posterior3)

#2.4#Ploting the posterior

plot(p_grid, posterior2, type = "b", 
     xlab = "Probability of water", ylab = "Posterior probability")

#2.4#Ploting the posterior

plot(p_grid, posterior3, type = "b", 
     xlab = "Probability of water", ylab = "Posterior probability")

#R code 2.6

library(rethinking)

globe.qa <- quap(
  
  alist(
    
    W ~ dbinom( W+L ,p) , # binomial likelihood
    
    p ~ dunif(0,1) # uniform prior
    
  ),
  
  data=list(W=6,L=3) )

# display summary of quadratic approximation

precis( globe.qa )

#R code 2.7

# analytical calculation

W <- 6

L <- 3

curve( dbeta( x , W+1 , L+1 ) , from=0 , to=1 )

# quadratic approximation
curve( dnorm( x , 0.67 , 0.16 ) , lty=2 , add=TRUE )


#R code 2.8

n_samples <- 1000

p <- rep( NA , n_samples )

p[1] <- 0.5

W <- 6

L <- 3

for ( i in 2:n_samples ) {
  
  p_new <- rnorm( 1 , p[i-1] , 0.1 )
  
  if ( p_new < 0 ) p_new <- abs( p_new )
  
  if ( p_new > 1 ) p_new <- 2 - p_new
  
  q0 <- dbinom( W , W+L , p[i-1] )
  
  q1 <- dbinom( W , W+L , p_new )
  
  p[i] <- ifelse( runif(1) < q1/q0 , p_new , p[i-1] )
  
}

#R code 2.9

dens( p , xlim=c(0,1) )

curve( dbeta( x , W+1 , L+1 ) , lty=2 , add=TRUE )
################################################################################
#################### 2.6 Pratice ###############################################
################################################################################
e1 = 4
e2 = 3
e3 = 1
e4 = "That if we throw up the globe and catch it without watchin where the index finger land on the globe, 
      we can be 70% certain that it is in an area with water"

########M1 - grid approximation of the globe tossing

#Making the grid
globe_toss_grid = seq(from=0, to=1, length.out = 1000)

#Uniform prior
globe_toss_prior = rep(1,1000)

#First sequence of tosses
seq1_likelihood = dbinom(3, size = 3, prob = globe_toss_grid)
seq2_likelihood = dbinom(3, size = 4, prob = globe_toss_grid)
seq3_likelihood = dbinom(5, size = 7, prob = globe_toss_grid)

#compute posteriors

#compute the product of likelihood and prior
unstd_posterior_seq1 = globe_toss_prior * seq1_likelihood
posterior_seq1 = unstd_posterior_seq1/sum(unstd_posterior_seq1)

#Sequence 2
unstd_posterior_seq2 = globe_toss_prior * seq2_likelihood
posterior_seq2 = unstd_posterior_seq2/sum(unstd_posterior_seq2)

#Sequence 3
unstd_posterior_seq3 = globe_toss_prior * seq3_likelihood
posterior_seq3 = unstd_posterior_seq3/sum(unstd_posterior_seq3)

par(mfrow=c(1,3))
#Plot
plot(globe_toss_grid, posterior_seq1, type = "b", col = "green",
     xlab = "Probability of water", ylab = "Posterior probability")


plot(globe_toss_grid, posterior_seq2, type = "b", col = "blue",
     xlab = "Probability of water", ylab = "Posterior probability")

plot(globe_toss_grid, posterior_seq3, type = "b", col = "red",
     xlab = "Probability of water", ylab = "Posterior probability")

###### M2 #################
globe_toss_prior = ifelse(globe_toss_grid < .5, 0, 1)

#compute posteriors

#compute the product of likelihood and prior
unstd_posterior_seq1 = globe_toss_prior * seq1_likelihood
posterior_seq1 = unstd_posterior_seq1/sum(unstd_posterior_seq1)

#Sequence 2
unstd_posterior_seq2 = globe_toss_prior * seq2_likelihood
posterior_seq2 = unstd_posterior_seq2/sum(unstd_posterior_seq2)

#Sequence 3
unstd_posterior_seq3 = globe_toss_prior * seq3_likelihood
posterior_seq3 = unstd_posterior_seq3/sum(unstd_posterior_seq3)

par(mfrow=c(1,3))
#Plot
plot(globe_toss_grid, posterior_seq1, type = "b", col = "green",
     xlab = "Probability of water", ylab = "Posterior probability")


plot(globe_toss_grid, posterior_seq2, type = "b", col = "blue",
     xlab = "Probability of water", ylab = "Posterior probability")

plot(globe_toss_grid, posterior_seq3, type = "b", col = "red",
     xlab = "Probability of water", ylab = "Posterior probability")

###### M3 ##############
p.land.mars = 1
p.land.earth = .3
p.planet.tossed = .5
p.land = (p.planet.tossed*p.land.earth)+(p.planet.tossed*p.land.mars)

p.earth.land = (p.land.earth* p.planet.tossed) / p.land
p.earth.land

###### M4 ############
#Number of ways each card can return a black side
bb =2
bw = 1
ww = 0

#Sum up the possible ways to get a back side
draw_black_side = bb + bw + ww

#Divide the ways two black sides with the possible ways to get one black side
second_black_given_first_black = bb /draw_black_side

######## M5 #######################

#Sum up the possible ways to get a back side
draw_black_side_2 = bb*2 + bw + ww

#Divide the ways two black sides with the possible ways to get one black side
second_black_given_first_black_2 = bb*2 /draw_black_side_2

######### M6 #######################

bw = bw * 2
ww = ww * 3

#Sum up the possible ways to get a back side
draw_black_side_3 = bb/ (bb + bw + ww)

draw_black_side_3

######## M7#####################

#Number of ways each card can return a black side
black_card = 3

white_card = 3



#There are 6 ways of getting black then white if the first card is bb
bw_with_bb_first = bb + bw+ white_card

#There are 8 possible ways total of getting the black white sequence 
bw_total = 8

bb_when_second_is_white = bw_with_bb_first / bw_total
bb_when_second_is_white


########## H1 ########################

#Twins given A
A_twins = 0.10

#Twins given B
B_twins = 0.20

#50% probability of the bear being A
p_A = .5
p_B = .5

#Probability of twins regardless species 
p_twins = A_twins* p_A + B_twins*p_A

#Probability of A given born twins
A_given_twins = (A_twins * p_A )/p_twins

#Probability of B given twins
B_given_twins = (B_twins * p_B)/ p_twins

#Probability of twins given that twins has already been born.´
p_twinsX2 = A_twins*A_given_twins + B_twins*B_given_twins
p_twinsX2

########### H2 ###############

#The probabilty of A given twins 
A_given_twins


########### H3 ###############

p_singleton = (1-A_twins)* A_given_twins + (1- B_twins)*B_given_twins

A_given_singleton = ((1- A_twins) *A_given_twins)/p_singleton
A_given_singleton

########### H4 ###################
test_A = .8
test_B = .65

test_pos= test_A*p_A + test_B*p_B

p_a_test_pos = (test_A * p_A) / test_pos



##With birth info
p_a = A_given_singleton
p_b = 1 -A_given_singleton

#Recalculate probability of positive test based on species distribution
test_pos_birth = test_A*p_a + test_B*p_b

p_a_test_pos_birth_infor = test_A*p_a / test_pos_birth

p_a_test_pos_birth_infor
