# define grid
data_len = 50

p_grid <- seq(0,1,length.out=data_len)

# define prior; We are using uniform distribution here
# prior <- rep(1,1000)
# prior <- ifelse(p_grid < 0.5, 0, 1)  
prior <- exp(-5*abs(p_grid - 0.5))

# compute likelihood at each value in grid
likelihood <- dbinom(6, 9, prob=p_grid)

# compute product of likelihood and prior
unstd_posterior <- likelihood * prior

# standardize the posterior, so it sums to 1
posterior <- unstd_posterior / sum(unstd_posterior)

# plot it 
plot(p_grid, posterior, type="b")


