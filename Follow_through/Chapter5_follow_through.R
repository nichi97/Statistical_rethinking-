library(rethinking)
library(tidyverse)
data(WaffleDivorce)
d <- WaffleDivorce

d$A <- scale(d$MedianAgeMarriage)
d$D <- scale(d$Divorce)

m5.1 <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + b_a * A,
    a ~ dnorm(0,0.2),
    b_a ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data=d
)

plot(d$A, d$D, xlab="Median Age Marriage", ylab="Divorce Rate")
title("Divorce Rate vs Median Age Marriage")

samples <- extract.samples(m5.1)
curve(mean(samples$a) + mean(samples$b_a) * (x - mean(d$A)), add=TRUE)

library(dagitty)
dag5.1 <- dagitty("dag{
  A -> D
  A -> M
  M -> D
}")
coordinates(dag5.1) <- list(x=c(A=0, D=1, M=2), y=c(A=0, D=1, M=0))
drawdag(dag5.1)

d$M <- scale(d$Marriage)
m5.3 <- quap(
  alist(
    D ~ dnorm(mu, sigma), 
    mu <- a + bM * M + bA * A,
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5),
    bA ~ dnorm(0, 0.5),
    sigma ~ dexp(1)
  ), data = d)
)
precis(m5.3)

plot(coeftab(m5.1, m5.3), par=c("bA", "bM"))


# Predictor residual plots

# To compute predictor residuals for either, we just use the other predictor to model it

# so for marriage rate, we use median marriage age to model it

m5.4 <- quap(
  alist(
    M ~ dnorm(mu, sigma),
    mu <- a + bAM * A, 
    a ~ dnorm(0,0.2),
    bAM ~ dnorm(0,0.5), 
    sigma ~ dexp(1)
  ), data = d
)

mu <- link(m5.4)
mu_mean <- apply(mu, 2, mean)
mu_resid <- d$M - mu_mean
d$M_res <- mu_resid

samples <- extract.samples(m5.4)

# Now do another linear regression with x as the marriage rate residuals and 
# y as the standardized divorce rate

m_res <- quap(
  alist(
    D ~ dnorm(mu, sigma),
    mu <- a + b * ( M_res - mean(M_res)), 
    a ~ dnorm(0,0.2),
    b ~ dnorm(0,0.5),
    sigma ~ dexp(1)
  ), data = d
)

samples <- extract.samples(m_res)

a <- mean(samples$a)
b <- mean(samples$b)

x_bar <- mean(d$M_res)

ggplot(data = d, aes(M_res, D)) + 
  geom_point() + 
  geom_abline(aes(intercept=a - b*x_bar, slope=b)) +
  xlab("Marriage rate residuals") + 
  ylab("Divorce rate (std)") 

# Now, turn to Posterior prediction plots

mu <- link(m5.3)

mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)

D_sim <- sim(m5.3, n=1e4)
D_PI <- apply(D_sim, 2, PI)

plot(mu_mean ~ d$D, col=rangi2, ylim=range(mu_PI), 
     xlab="Observed divorce", ylab="Predicted divorce")
abline(a=0, b=1, lty=2)
for(i in 1:nrow(d)) lines(rep(d$D[i],2), mu_PI[,i], col=rangi2)

identify(x=d$D, y=mu_mean, labels=d$Loc)


# counterfactual model

data("WaffleDivorce")
d <- list()
d$A <- standardize(WaffleDivorce$MedianAgeMarriage)
d$D <- standardize(WaffleDivorce$Divorce)
d$M <- standardize(WaffleDivorce$Marriage)

m5.3_A <- quap(
  alist(
    ## A -> D <- M
    D ~ dnorm(mu, sigma), 
    mu <- a + bM*M + bA*A, # because there is no interaction between A and M
    a ~ dnorm(0, 0.2),
    bM ~ dnorm(0, 0.5), 
    bA ~ dnorm(0, 0.5), 
    sigma ~ dexp(1), 
    
    ## A -> M
    M ~ dnorm(mu_M, sigma_M), 
    mu_M <- aM + bAM*A, 
    aM ~ dnorm(0,0.2), 
    bAM ~ dnorm(0, 0.5), 
    sigma_M ~ dexp(1)
    
  ), data =d )

A_seq <- seq(from=-2, to=2, length.out=30)

sim_dat <- data.frame( A=A_seq)
s <- sim(m5.3_A, data=sim_dat, var=c("M", "D"))

str(s)

# display counterfactual predictions
plot(sim_dat$A, colMeans(s$D), ylim=c(-2,2), type="l", 
     xlab="mnipulated A", ylab="counterfactual D")
shade(apply(s$D, 2, PI), sim_dat$A)
mtext("Total counterfactual effect of A on D")


plot(sim_dat$A, colMeans(s$M), ylim=c(-2,2), type="l", 
     xlab="manipulated A", ylab="counterfactual M")
shade(apply(s$M, 2, PI), sim_dat$A)
mtext("Counterfactual effect A -> M")


# simulate the effect of A -> D <- M, but assume no relationship between 
# A and M

sim_dat <- data.frame(M=seq(from=-2, to=2, length.out=30), A=0)
s <- sim(m5.3_A, data=sim_dat, vars="D")

plot(sim_dat$M, colMeans(s), ylim=c(-2,2), type="l", 
     xlab="manipulated M", ylab="counterfactual D")
shade(apply(s,2,PI), sim_dat$M)
mtext("Total counterfactual effect of M on D")





######################### 5.2 Masked relationship

library(rethinking)
data(milk)
d <- milk
str(d)

d$K <- scale(d$kcal.per.g)
d$N <- scale(d$neocortex.perc)
d$M <- scale(log(d$mass))

# first model to consider -- bivariate regression 

# deal with NA value, use complete.cases
dcc <- d[complete.cases(d$K, d$N, d$M), ]

m5.5_draft <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <-  a + bN * N,
    a ~ dnorm(0,1), 
    bN ~ dnorm(0,1), 
    sigma ~ dexp(1)
    
  ), data=dcc )

# do prior prediction to see if things makes sense

prior <- extract.prior(m5.5_draft)
xseq <- c(-2,2)
mu <- link(m5.5_draft, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for(i in 1:50) lines(xseq, mu[i, ], col=col.alpha("black", 0.3))

# but the prior graph above is absolute awful. Change the prior so that it 
# looks better


m5.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <-  a + bN * N,
    a ~ dnorm(0,0.2), 
    bN ~ dnorm(0,0.5), 
    sigma ~ dexp(1)
    
  ), data=dcc )

# do prior prediction to see if things makes sense

prior <- extract.prior(m5.5)
xseq <- c(-2,2)
mu <- link(m5.5_draft, post=prior, data=list(N=xseq))
plot(NULL, xlim=xseq, ylim=xseq)
for(i in 1:50) lines(xseq, mu[i, ], col=col.alpha("black", 0.3))

precis(m5.5)

xseq <- seq(from=min(dcc$N)-0.15, to=max(dcc$N)+0.15, length.out=30)
mu <- link(m5.5, data=list(N=xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(K ~ N, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

m5.6 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bM * M ,
    a ~ dnorm(0, 0.2) , 
    bM ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ) , data=dcc
)
precis(m5.6)

# this graph has a negative slope
xseq <- seq(from=min(dcc$M)-0.15, to=max(dcc$M)+0.15, length.out=30)
mu <- link(m5.6, data=list(M=xseq))
mu_mean <- apply(mu, 2, mean)
mu_PI <- apply(mu, 2, PI)
plot(K ~ M, data=dcc)
lines(xseq, mu_mean, lwd=2)
shade(mu_PI, xseq)

# Now let's see what happens when we add both predictors
m5.7 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bN * N + bM * M, 
    a ~ dnorm(0, 0.2), 
    bN ~ dnorm(0, 0.5), 
    bM ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data=dcc
)
precis(m5.7)

plot(coeftab(m5.5, m5.6, m5.7), pars=c("bM", "bN"))

pairs( ~K + M + N, dcc)


# Categorical data analysis

data("Howell1")
d <- Howell1

d$sex <- ifelse(d$male==1, 2, 1)
str(d$sex)

m5.8 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a[sex], 
    a[sex] ~ dnorm(178, 20), 
    sigma ~ dunif(0, 50)
  ), data=d )

precis(m5.8, depth=2)

post <- extract.samples(m5.8)  
post$diff_fm <- post$a[,1] - post$a[,2]
precis(post, depth=2, hist=FALSE)
  
# multipel categories using index method

data(milk)
d <- milk
unique(d$clade)

d$clade_id <- as.integer(d$clade )
  

d$K <- scale(d$kcal.per.g)
m5.9 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a[clade_id], 
    a[clade_id] ~ dnorm(0,0.5), 
    sigma ~ dexp(1)
  ), data=d )

labels <- paste( "a[" , 1:4 , "]:" , levels(d$clade) , sep="" )
plot( precis( m5.9 , depth=2 , pars="a" ) , labels=labels ,
      xlab="expected kcal (std)" ) 

  
  
  
  





