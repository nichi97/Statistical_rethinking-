###################################
#       Problem 1                 #
###################################


# So, let's try grid with 100 points
grid <- seq(0,1,length.out=100)

# Then compute value of prior
prior <- rep(1, 100)

# calculate likelihood using binomial model, calculate prob at each grid point to find 
# an approximation of the likelihood function
likelihood <- dbinom(8,15, grid)

# calculate posterior 
posterior_unstd <- likelihood * prior

posterior <- posterior_unstd / sum(posterior_unstd)

plot(grid, posterior, type="b")



###########################################
#               Problem 2                 #
###########################################



# So, let's try grid with 100 points
grid <- seq(0,1,length.out=100)

# Then compute value of prior
prior <- ifelse(grid > 0.5, 2, 0)

# calculate likelihood using binomial model, calculate prob at each grid point to find 
# an approximation of the likelihood function
likelihood <- dbinom(8,15, grid)

# calculate posterior 
posterior_unstd <- likelihood * prior

posterior <- posterior_unstd / sum(posterior_unstd)

plot(grid, posterior, type="b")


##################################################
#                Problem 3                       #
##################################################




# So, let's try grid with 100 points
grid <- seq(0,1,length.out=1000)

# Then compute value of prior
prior <- ifelse(grid > 0.5, 2, 0)

# calculate likelihood using binomial model, calculate prob at each grid point to find 
# an approximation of the likelihood function
likelihood <- dbinom(1600,2000, grid)

# calculate posterior 
posterior_unstd <- likelihood * prior

posterior <- posterior_unstd / sum(posterior_unstd)

plot(grid, posterior, type="b")


