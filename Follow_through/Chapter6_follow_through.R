library(rethinking)
# simulate newsworthiness and trustworthiness

set.seed(1914)
N <- 200
p <- 0.1

nw <- rnorm(N)
tw <- rnorm(N)

s <- nw + tw
q <- quantile(s, 1-p)
selected <- ifelse(s >= q, TRUE, FALSE)
cor(tw[selected], nw[selected])

# multicollinear legs

N <- 100
#set.seed(909)
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5)
leg_left <- leg_prop * height + rnorm(N, 0, 0.02)
leg_right <- leg_prop * height + rnorm(N, 0, 0.02)

d <- data.frame(height, leg_left, leg_right)

m6.1 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + bl*leg_left + br*leg_right, 
    a ~ dnorm(10, 100), 
    bl ~ dnorm(2, 10), 
    br ~ dnorm(2, 10), 
    sigma ~ dexp(1)
  ), data = d
)

precis(m6.1)

plot(precis(m6.1))


post <- extract.samples(m6.1)
plot(bl ~ br, post, col=col.alpha(rangi2, 0.1), pch=16)

sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd=2, xlab="sum of bl and br")


m6.2 <- quap(
  alist(
    height ~ dnorm(mu, sigma), 
    mu <- a + bl*leg_left, 
    a ~ dnorm(10, 100), 
    bl ~ dnorm(2, 10), 
    sigma ~ dexp(1)
  ), data = d
)
precis(m6.2)



####### Milk example #########

library(rethinking)
data(milk)
d <- milk
d$K <- scale( d$kcal.per.g)
d$F <- scale(d$perc.fat)
d$L <- scale(d$perc.lactose)

m6.3 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bF * F,
    a ~ dnorm(0, 0.2), 
    bF ~ dnorm(0, 0.5), 
    sigma ~ dexp(1) 
  ), data = d
)

m6.4 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bL * L, 
    a ~ dnorm(0, 0.2), 
    bL ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d
)

precis(m6.3)
precis(m6.4)

pairs(kcal.per.g ~ perc.fat + perc.lactose, data=d)


m6.5 <- quap(
  alist(
    K ~ dnorm(mu, sigma), 
    mu <- a + bF*F + bL*L, 
    a ~ dnorm(0, 0.2), 
    bF ~ dnorm(0, 0.5), 
    bL ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), 
  data = d
)

precis(m6.5)

pairs( ~ kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)

################## Post-Treatment Bias #################
set.seed(71)

N <- 100

h0 <- rnorm(N, 10, 2)

treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4)
h1 <- h0 + rnorm(N, 5-3*fungus)

# compose a clean data frame
d <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)
precis(d)

sim_p <- rlnorm(1e4, 0, 0.25)
precis(data.frame(sim_p))

m6.6 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0*p, 
    p ~ dlnorm(0, 0.25), 
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.6)

m6.7 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt * treatment + bf * fungus, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    bf ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data=d
)

precis(m6.7)

library(dagitty) 
plant_dag <- dagitty("dag{
                     H_0 -> H_1 
                     F -> H_1
                     T -> F
}")
coordinates(plant_dag) <- list(x=c(H_0=0, T=2, F=1.5, H_1=1), 
                               y=c(H_0=0, T=0, F=0, H_1=0))
drawdag(plant_dag)

impliedConditionalIndependencies(plant_dag)


# with moisture as the new variable
set.seed(71)
M <- 1000
h0 <- rnorm(N, 10, 2)
treatment <- rep(0:1, each=N/2)
M <- rbern(N)
fungus <- rbinom(N, size=1, prob=0.5-treatment*0.4 + 0.4*M)
h1 <- h0 + rnorm(N, 5 + 3*M)
d2 <- data.frame(h0=h0, h1=h1, treatment=treatment, fungus=fungus)

m6.8 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt*treatment, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d
)


m6.9 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt * treatment + bf * fungus, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    bf ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data=d2
)

precis(m6.8)


m6.10 <- quap(
  alist(
    h1 ~ dnorm(mu, sigma), 
    mu <- h0 * p, 
    p <- a + bt*treatment, 
    a ~ dlnorm(0, 0.2), 
    bt ~ dnorm(0, 0.5), 
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.10)

#################### collider #############################

library(rethinking)
d <- sim_happiness(seed=1997, N_years=1000)
precis(d, hist=F)

d2 <- d[d$age > 17, ]
d2$A <- (d2$age - 18) / (65-18)

d2$mid <- d2$married + 1
m6.11 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma), 
    mu <- a[mid] + bA * A, 
    a[mid] ~ dnorm(0,1), 
    bA ~ dnorm(0,2), 
    sigma ~ dexp(1)
  ), data = d2
)

precis(m6.11, depth=2)

m6.12 <- quap(
  alist(
    happiness ~ dnorm(mu, sigma), 
    mu <- a + bA * A, 
    a ~ dnorm(0, 1), 
    bA ~ dnorm(0,2), 
    sigma ~ dexp(1)
  ), data = d2
)
precis(m6.12)

# Bad grandparents

N <- 200
b_GP <- 1
b_GC <- 0
b_PC <- 1
b_U <- 2

set.seed(1)
U <- 2 * rbern(N, 0.5) - 1
G <- rnorm(N)
P <- rnorm(N, b_GP*G + b_U*U)
C <- rnorm(N, b_PC * P + b_GC*G + b_U*U)
d <- data.frame(C=C, P=P, G=G, U=U)

m6.13 <- quap(
  alist(
    C ~ dnorm(mu, sigma), 
    mu <- a + b_PC*P + b_GC*G, 
    a ~ dnorm(0, 1), 
    c(b_PC, b_GC) ~ dnorm(0,1), 
    sigma ~ dexp(1)
  ), data=d
)
precis(m6.13)









