library(rethinking)

globe.qa <- quap(
  alist(
    W ~ dbinom(W + L, p),
    p ~ dunif(0,1)
  ),
  data=list(W=6, L=3)
)

precis(globe.qa)

W <- 6
L <- 3
curve(dbeta(x, W+1, L+1), from=0, to=1)
curve(dnorm(x, 0.67, 0.16), lty=2, add=TRUE)

# posterier predictive distribution

w <- rbinom(1e4, size=9, prob=samples)

simplehist(w)

