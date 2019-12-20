#2M1
# define the grid
p_grid <- seq(0,1,length.out=1000)

# compute the prior at each grid point
# p ~ unif(0,1)
# W ~ bino(N,p)
p <- rep(1,1000)
W <- dbinom(5, 7, prob=p_grid)

# now compute the likelihood
posterior_raw <- p * W
posterior <- posterior_raw / sum(posterior_raw)

# then sample the posterior
samples <- sample(p_grid, size=1e4, replace = TRUE, prob = posterior)

dens(samples)



#2M2
# define the grid
p_grid <- seq(0,1,length.out=1000)

# compute the prior at each grid point
# p ~ unif(0,1)
# W ~ bino(N,p)
p <- ifelse(p_grid < 0.5, 0, 2)
W <- dbinom(3, 3, prob=p_grid)

# now compute the likelihood
posterior_raw <- p * W
posterior <- posterior_raw / sum(posterior_raw)

# then sample the posterior
samples <- sample(p_grid, size=1e4, replace = TRUE, prob = posterior)

dens(samples)

