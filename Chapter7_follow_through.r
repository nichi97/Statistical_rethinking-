library(rethinking)
sppnames <- c("afarensis", "africanus", "habilis", "boisei", "rudolfensis", "ergaster", "sapiens")
brainvolcc <- c(438, 452, 612, 521, 752, 871, 1350)
masskg <- c(37.0, 35.5, 34.5, 41.5, 55.5, 61.0, 53.5)
d <- data.frame(species = sppnames, brain = brainvolcc, mass = masskg)


# standardize parameters
d$mass_std <- (d$mass - mean(d$mass)) / sd(d$mass)
d$brain_std <- d$brain / max(d$brain)

m7.1 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b * mass_std,
    a ~ dnorm(0.5, 1), 
    b ~ dnorm(0, 10),
    log_sigma ~ dnorm(0,1)
  ), data = d)

set.seed(12)

s <- sim(m7.1) # generate prediction
r <- apply(s, 2, mean) - d$brain_std # get residual for each prediction
resid_var <- var2(r) # find out the variance of the residual
outcome_var <- var2(d$brain_std) # find out the variance
1 - resid_var/outcome_var # r^2 is 0.4774589

R2_is_bad <- function(quap_fit){
  s <- sim(quap_fit, refresh=0)
  r <- apply(s,2,mean) - d$brain_std
  l - var2(r)/var2(d$brain_std)
}

m7.2 <- quap(
  alist(
    brain_std ~ dnorm(mu, exp(log_sigma)), 
    mu <- a + b[1]*mass_std + b[2]*mass_std^2,
    a ~ dnorm(0.5,1),
    b ~ dnorm(0,10), 
    log_sigma ~ dnorm(0,1)
  ), data = d, start=list(b=rep(0,2))
)


m7.3 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,3)) )
m7.4 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,4)) )
m7.5 <- quap(
  alist(
    brain_std ~ dnorm( mu , exp(log_sigma) ),
    mu <- a + b[1]*mass_std + b[2]*mass_std^2 +
      b[3]*mass_std^3 + b[4]*mass_std^4 +
      b[5]*mass_std^5,
    a ~ dnorm( 0.5 , 1 ),
    b ~ dnorm( 0 , 10 ),
    log_sigma ~ dnorm( 0 , 1 )
  ), data=d , start=list(b=rep(0,5)))


m7.6 <- quap(
  alist(
    brain_std ~ dnorm( mu, 0.001), 
    mu <- a + b[1] * mass_std + b[2] * mass_std^2 + b[3] * mass_std^3 + b[4] * mass_std^4 +
      b[5] * mass_std^5 + b[6] * mass_std^6, 
    a ~ dnorm(0.5, 1), 
    b ~ dnorm(0, 10) ), 
  data = d, start = list(b=rep(0,6))
)

lppd(m7.1)

sapply(list(m7.1, m7.2, m7.3, m7.4, m7.5, m7.6), function(m) sum(lppd(m)))


library(rethinking)
data(WaffleDivorce)
d <- WaffleDivorce
d$A <- standardize( d$MedianAgeMarriage )
d$D <- standardize( d$Divorce )
d$M <- standardize( d$Marriage )
m5.1 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bA * A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )
m5.2 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM * M ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )
m5.3 <- quap(
  alist(
    D ~ dnorm( mu , sigma ) ,
    mu <- a + bM*M + bA*A ,
    a ~ dnorm( 0 , 0.2 ) ,
    bM ~ dnorm( 0 , 0.5 ) ,
    bA ~ dnorm( 0 , 0.5 ) ,
    sigma ~ dexp( 1 )
  ) , data = d )

set.seed(24071847)
compare(m5.1, m5.2, m5.3, func=PSIS)

PSIS_m5.3 <- PSIS(m5.3, pointwise=TRUE)
WAIC_m5.3 <- WAIC(m5.3, pointwise=TRUE)
plot(PSIS_m5.3$k, WAIC_m5.3$penalty, xlab="PSIS Pareto k", ylab="WAIC penalty", col=rangi2, lwd=2)


m5.3t <- quap(
  alist(
    D ~ dstudent(2, mu, sigma), 
    mu <- a + bM * M + bA * A, 
    a ~ dnorm(0, 0.2), 
    bM ~ dnorm(0, 0.5), 
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d
)

PSIS(m5.3t)

PSIS(m5.3)







