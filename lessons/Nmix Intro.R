
# We'll start by simulating the data for a basic N-mixture model. In our simulation, abundance
# will be affected by a variable called "vegHt" (vegetation height) and detection will be affected
# by a variabled called "wind" (wind speed).

#### SIMULATING ABUNDANCE ####  

# We start by simulating the number of sites we need
nSites <- 100
set.seed(443)  # we use this function so that we all get the same values for vegHt when simulating random variables

# Create a covariate called vegHt
vegHt <- rnorm(nSites, 10, 3) # Normal distribution with mean 10 and variation 3

# We now standardize vegHt for use in our model (set to mean 0 with sd 1)
vegHt <- scale(vegHt)

# We now need to create our linear model for abundance with vegHt
# The relationship is described by an intercept of -1 and
#    a slope parameter of 2 on the log scale
lambda <- exp(-1 + 2*vegHt)

# Now we simulate abundace at each site. The 
N <- rpois(nSites, lambda)

plot(N)
plot(vegHt,N)

# We can fit a model without detection probability that relates abundance to vegHt 
# using the glm() function with "family=Poisson":

N.glm <- summary(fm.glm1 <- glm(N ~ vegHt, family=poisson))
N.glm

# Do some analysis of the results
plot(vegHt, N, xlab="Vegetation height", ylab="Abundance (N)")

glm1.est <- coef(fm.glm1) #returns intercept and vegHt beta coefficients. should be around -1 and 2.


plot(function(x) exp(-1 + 2*x), 1, 3, add=TRUE, lwd=3)
plot(function(x) exp(glm1.est[1] + glm1.est[2]*x), 1, 3, add=TRUE,
     lwd=3, col="blue")
legend(1, 20, c("Truth", "Estimate"), col=c("black", "blue"), lty=1,
       lwd=3)


#### ADDING DETECTION PROBABILITY ####

#Now we incorporate the detection/observation portion
nVisits <- 3    # number of repeat samples at each site

wind <- rnorm(300,15,5)
wind <- scale(wind)

#now use inverse logit to calculate linear model for detection probability
p <- 1/(1+exp(-(0.5 + 1*wind)))
plot(p,wind)

#arrange p and wind into dataframes
p <- data.frame(J1=p[1:100],J2=p[101:200],J3=p[201:300])
wind.df <- data.frame(wind1=wind[1:100],wind2=wind[101:200],wind3=wind[201:300])

#create matrix to hold each observation at each site for each visit
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
  for(j in 1:nVisits){
  y[i,j] <- rbinom(1, N[i], p[i,j])
}}

# Look at the data
cbind(N=N, y1=y[,1], y2=y[,2], y3=y[,3] ,p1=p[,1], p2=p[,2],p3=p[,3])

# Load library, format data and summarize
library(unmarked)

#Create data structure unique to unmarked package. Inputs data (y) and site covariates (siteCovs)
umf <- unmarkedFramePCount(y=y, siteCovs=as.data.frame(vegHt), 
                           obsCovs=list(wind=wind.df))
summary(umf)

# Fit a model and extract estimates
# When writing formula, detection covariates follow first tilde, then come abundance covariates

#Note that we're using the pcount function for an N-mixture model
fm.nmix1 <- pcount(~wind ~vegHt, data=umf) 
summary(fm.nmix1)

# Note, estimates of detection coefficients are on the logit-scale
# When covariates are in the model we can use the following to
# backtransform them
beta1 <- coef(fm.nmix1) #coef() extracts coefficients from the nmix model
beta1

#We can now plot the predicted response between the covariates and 
#abundance
veg.plot=seq(min(vegHt),max(vegHt),length=100)

plot(veg.plot,exp(beta1[1] + beta1[2]*veg.plot),type="l",
     xlab="vegetation height", ylab="Expected Abundance")

#or covariates and detection probability (logit link)
wind.plot=seq(min(wind),max(wind),length=300)
plot(wind.plot,(1/(1+exp(-(beta1[3]+beta1[4]*wind.plot)))),type="l",
     xlab="Wind Speed",ylab="Detection Probability")

# Or suppose you want predictions for new values of vegHt, say 1.2 and 3.1
newdat <- data.frame(vegHt=c(3.1,5))
predict(fm.nmix1, type="state", newdata=newdat)

# We can calculate posterior distributions for latent abundance at
# each site as well.
site.N <- ranef(fm.nmix1)
# and get confidence intervals for each of those estimates
confint(site.N,level=0.95)

# Or compare the estimate the overall total abundance to truth
sum(N) #original simulated abundance
sum(bup(site.N)) #total estimated abundance

# Analyze with real data
nobo <- read.csv("nobo_abund.csv")
head(nobo)

# Create unmarked dataframe
#Create data structure unique to unmarked package. Inputs data (y) and site covariates (siteCovs)
nobo.umf <- unmarkedFramePCount(y=nobo[,2:4], 
            siteCovs=as.data.frame(scale(nobo[,5:10])), 
            obsCovs=list(sky=scale(nobo[,11:13]),jdate=scale(nobo[,14:16]),
                         time=scale(nobo[,17:19])))
summary(nobo.umf)

# Fit a model to evaluate covariate affects on NOBO abundance
nobo.1 <- pcount(~sky + jdate + time ~BA + Evergreen5km, 
                 data=nobo.umf,K=105) 
summary(nobo.1)

nobo.coef <- coef(nobo.1) #coef() extracts coefficients from the nmix model
nobo.coef

#We can now plot the predicted response between the covariates and 
#abundance
BA.plot=seq(min(scale(nobo[,6])),max(scale(nobo[,6])),length=100)

plot(BA.plot,exp(nobo.coef[1] + nobo.coef[2]*BA.plot),type="l",
     xlab="Basal Area", ylab="Expected Abundance")

sky.plot=seq(min(scale(nobo[,11:13])),max(scale(nobo[,11:13])),length=100)

plot(sky.plot,(1/(1+exp(-(nobo.coef[3]+nobo.coef[4]*sky.plot)))),type="l",
     xlab="Sky Cover",ylab="Detection Probability")
