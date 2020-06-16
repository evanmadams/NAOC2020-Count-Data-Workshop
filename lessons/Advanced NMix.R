# More advanced tools for N-mixture models

#### DISTRIBUTIONS FOR ABUNDANCE ####
# Can use zero-inflated Poisson (ZIP) or negative binomial as well

fm.nmix2 <- pcount(~wind ~vegHt, data=umf,mixture="NB")
summary(fm.nmix2)
fm.nmix3 <- pcount(~wind ~vegHt, data=umf,mixture="ZIP")
summary(fm.nmix3)

# Example with real data
nobo.nb <- pcount(~sky ~BA, data=nobo.umf, mixture="NB")
summary(nobo.nb)
summary(nobo.1)

nobo.zip <- pcount(~sky ~BA, data=nobo.umf, mixture="ZIP")
summary(nobo.zip)

nobo.pois <- pcount(~sky ~BA, data=nobo.umf, mixture="P")

#### MODEL COMPARISON ####
## We can compare models with different covariates using AIC
cbind(AIC.1=nobo.1@AIC, AIC.pois=nobo.pois@AIC)

## We can also compare between different data models
cbind(AIC.zip=nobo.zip@AIC, AIC.nb=nobo.nb@AIC, AIC.pois=nobo.pois@AIC)

# Goodness-of-fit statistics
# Is the best model from AIC comparison also a good model?

library(AICcmodavg)
Nmix.gof.test(nobo.1,nsim=10) #looking for p > 0.05

#### OPEN N-MIXTURE MODEL ####

# Simulate data under the "constant" model.
#    Basic birth-death process but birth rate isn't affected by
#    abundance in previous year. Not realistic, but it is easy to extend

lam <- 3
omega <- 0.5
gamma <- 2
p <- 0.8

nSites <- 100
nYears <- 20
y <- N <- matrix(NA, nSites, nYears)
S <- G <- matrix(NA, nSites, nYears-1)

N[,1] <- rpois(nSites, lam)
for(t in 2:nYears) {
  S[,t-1] <- rbinom(nSites, N[,t-1], omega)
  G[,t-1] <- rpois(nSites, gamma)
  N[,t] <- S[,t-1] + G[,t-1]
}
y[] <- rbinom(nSites*nYears, N, p)

plot(1:nYears, colSums(N), xlab="Year", ylab="Population size",
     type="o")

# Fit the model in unmarked

umf.const <- unmarkedFramePCO(y=y, numPrimary=nYears)
fm.const <- pcountOpen(~1, ~1, ~1, ~1, umf.const, dynamics="constant",K=50)
summary(fm.const)

