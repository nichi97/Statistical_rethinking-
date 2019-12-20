curve(exp(-x^2), from=-3, to=3)

library(rethinking)

data(Howell1)
d <- Howell1

str(d)
precis(d)

d2 <- d[d$age >= 18, ]

dens(d2$height)

curve(dnorm(x,178,20), from=100, to=250)

sample_mu <- rnorm(1e4, 178, 100)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)

# using grid_approximation 
mu_list <- seq(from=150, to=160, length.out=100)
sigma_list <- seq(from=7, to=9, length.out=100)
post <- expand.grid(mu=mu_list, sigma=sigma_list)
post$LL <- sapply( 1:nrow(post), function(i) sum(
  dnorm(d2$height, post$mu[i], post$sigma[i], log=TRUE)))
post$prod <- post$LL + dnorm(post$mu, 178, 20, TRUE) + 
  dunif(post$sigma, 0, 50, TRUE)
post$prob <- exp(post$prod - max(post$prod))

sample_rows <- sample(1:nrow(post), size=1e4, replace=TRUE, prob=post$prob)
sample_mu <- post$mu[sample_rows]
sample_sigma <- post$sigma[sample_rows]

plot(sample_mu, sample_sigma, cex=0.5, pch=16, col=col.alpha(rangi2, 0.1))


# Quadratic Approximation
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- d[d$age >= 18, ]

flist <- alist(
  height ~ dnorm(mu, sigma),
  mu ~ dnorm(178, 20),
  sigma ~ dunif(0,50)
)

m4.1 <- quap(flist, data=d2)

vcov(m4.1)

post <- extract.samples(m4.1, n=1e4)
head(post)


# Linear Regression
plot(d2$height ~ d2$weight)


set.seed(2971)
N <- 100
a <- rnorm(N, 178,20)
b <- rnorm(N, 0, 10)

plot(NULL, xlim=range(d2$weight), ylim=c(-100,400), xlab="Weight", ylab="Height")
abline(h=0, lty=2) # dashed line
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)

for (i in 1:N) {
  curve( a[i] + b[i]*(x- xbar), from=min(d2$weight), to=max(d2$weight), 
        add=TRUE, col=col.alpha("black", 0.2))
}



# try out log normal 
b <- rlnorm(1e4, 0,1)
dens(b, xlim=c(0,5), adj=0.1)




set.seed(2971)
N <- 100
a <- rnorm(N, 178,20)
b <- rlnorm(N, 0, 1)

plot(NULL, xlim=range(d2$weight), ylim=c(-100,400), xlab="Weight", ylab="Height")
abline(h=0, lty=2) # dashed line
abline(h=272, lty=1, lwd=0.5)
mtext("b ~ dnorm(0,10)")
xbar <- mean(d2$weight)

for (i in 1:N) {
  curve( a[i] + b[i]*(x- xbar), from=min(d2$weight), to=max(d2$weight), 
        add=TRUE, col=col.alpha("black", 0.2))
}


# compute the posterior distribution
library(rethinking)
data(Howell1)
d <- Howell1
d2 <- Howell1[Howell1$age >= 18, ]

# Define the average weight, x-bar
xbar <- mean(d2$weight)

#fit model
m4.3 <- quap(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b * (weight - xbar),
    a ~ dnorm(178, 20),
    b ~ dlnorm(0,1),
    sigma ~ dunif(0,50)
  ), data = d2
)

extract.samples(m4.3)


plot(height ~ weight, data=d2, col=rangi2)
post <- extract.samples(m4.3)
a_map <- mean(post$a)
b_map <- mean(post$b)
curve(a_map + b_map * (x-xbar), add=TRUE)

# start with the first 10 cases in d2
N <- 352
dN <- d2[1:N, ]
mN <- quap(
  alist(
    height ~ dnorm( mu, sigma ),
    mu <- a + b * (weight - mean(weight)),
    a ~ dnorm(178, 20), 
    b ~ dlnorm(0,1),
    sigma ~ dunif( 0 , 50)
  ), data=dN )

post <- extract.samples(mN, n=20)

# display raw data and sample size
plot( dN$weight, dN$height, 
      xlim=range(d2$weight), ylim=range(d2$height),
      col=rangi2, xlab="weight", ylab="height")
mtext(concat("N = ", N))

for( i in 1:20)
  curve( post$a[i] + post$b[i] * (x-mean(dN$weight)), 
         col=col.alpha("black", 0.3), add=TRUE)

post <- extract.samples(m4.3)
mu_at_50 <- post$a + post$b * (50 - xbar)

dens(mu_at_50, col=rangi2, lwd=2, xlab="mu|weight=50")
mu_pi <- PI(mu_at_50)

# use link function to repeat the above steps 
mu <- link( m4.3 )

# define sequence of weights to compute predictions for 
# these values will be on the horizontal axis
weight_seq <- seq(from=25, to=70, by=1)

# use link to compute mu
# for each sample from posterior
# and for each weight in weight_seq
mu <- link( m4.3, data=data.frame(weight=weight_seq))
str(mu)

plot(height ~ weight, d2, type="n")

for(i in 1:100)
  points(weight_seq, mu[i, ], pch=16, col=col.alpha(rangi2, 0.1))

# summarize the distribution of mu
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI, prob=0.89)

# plot raw data
# fade out points to maek line and interval more visible
plot( height ~ weight, data=d2, col=col.alpha(rangi2, 0.5))

# plot the MAP line, aka the mean mu for each weight
lines(weight_seq, mu_mean)
shade(mu_PI, weight_seq)

sim_height <- sim(m4.3, data=list(weight=weight_seq, n=1e4))
str(sim_height)

height_PI <- apply(sim_height, 2, PI, prob=0.89)
  
plot(height ~ weight, d2, col=col.alpha(rangi2, 0.5))
  
lines(weight_seq, mu_mean)
shade(height_PI, weight_seq)  


# 4.5.1
