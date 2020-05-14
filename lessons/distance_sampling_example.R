#Count data workshop distance sampling example

#we want to survey the different ways that count data become biased and this is an example that occurs frequently with bird survey data
#when you survey an area you are often better at survey the area closest to you rather than further away
#distance sampling was created to try and understand how the chances of detecting a bird vary with their distance from the observer
#So we've added an additional source of information to the count data, for every count we have a distance to the observer
#We use the distance data to make another model that we combine with a simple Poisson regression model so that we can describe the patterns of abundance after accounting for detection bias
#once we estimate this source of bias, then we can account for the birds that we didn't see and correct our estimates

#we are using three packages today, unmarked is the package that does the distance sampling, AICcmodavg gives us some additional help with unmarked, reshape2 and dplyr help format data
library(unmarked)
library(AICcmodavg)
library(reshape2)
library(dplyr)

#let's pull in some data on Chipping Sparrow abundance from the Maine Bird Atlas
dat <- dget('data/mebba_spp_Chipping Sparrow.dat')

#these data are in a bit of weird format so I'm going to reorganize them into something that unmarked will recognize

#we create a list of sites that were visited
y <- data.frame(site = 1:dat$nsites)

#then we pull all the distance classess of our CHSP observation and the site that they below to
ds <- data.frame(site = dat$site.idx, dc = dat$distclass)

#here we have 5 distance classes: 20m intervals from 0-100m

#then we reorganize the data so I can count the number of individuals in each distance category in each site
ds <- dcast(ds, site ~ dc, fun.aggregate = length)

#finally we join back with the complete site list so that we know the sites that were surveyed
y <- left_join(y, ds, 'site')

#turn NAs in the merge into true zeroes as NA means that no birds from this species was detected
y[is.na(y)] <- 0

#I have a couple site covariates that I want to grab to describe either abundance or detection probability

sc <- data.frame(dat$habcovs)

#then we need to convert it to a format that unmarked understands
#unmarked has a function that helps us do that specifically for distance sampling data

umf <- unmarkedFrameDS(y = as.matrix(y[1:100 ,2:6]), siteCovs = sc[1:100,], dist.breaks = c(0, 20, 40, 60, 80, 100), survey = 'point', unitsIn = 'm')   #parse the data into the unmarked data frame and subset the data to (1) have models converge faster and (2) remove birds with unknown distances

#now we want to run the distance model using 'distsamp'
#the first element of the function needs a model formula that is similar to what we were using only we need multiple formulas for the distance model and the ecological model
#here we are running a null model with covariates on either submodel
#we also designate how we want to model the detection probability over distance (a half-normal distribution in this case)
fm <- distsamp(~1 ~1, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm)


#let's quantify the bias we saw in the counts
#to do this we'll look at the detection portion of the model to see how detectability varied with distance to the observer
backTransform(fm, type = 'det')   #backtransform the model estimates
hist(fm, xlab = 'Distance (m)')   #look at how well the model fits our distance data

#to determine how many birds we miss in a 100m radius point count circle, first we calculate the survey half-width

ehw <- integrate(gxhn, 0, 100, sigma = 39.1)$value  #calculate the survey half-width
ehw/100   #100 is the maximum survey distance and used to scale the half-width and calculate the proportion of birds detected




#look at our estimates of density
backTransform(fm, type = 'state')  #average across all sites
ranef(fm, K = 50)   #Empirical Bayes estimates for each site, K should be a value so high that the number of individuals at that site couldn't reach it

plot(ranef(fm, K = 50))

#we can also add covariates to either the abundance or detection models to test hypotheses about the importance of those variables

fm <- distsamp(~eff.domhabA ~eff.domhabF, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm)

#here we find that forest habitat is not strongly correlated with CHSP abundance but that it's easier to detect CHSP in developed habitat


#CHALLENGE

#make a few models and compare different detection functions

#note use 'aictab' in the AICcmodavg library to compare multiple models



#answer

fm1 <- fm <- distsamp(~eff.domhabA ~eff.domhabF, data = umf, keyfun = 'halfnorm', output = 'density')
fm2 <- fm <- distsamp(~eff.domhabA ~eff.domhabF, data = umf, keyfun = 'exp', output = 'density')
fm3 <- fm <- distsamp(~eff.domhabA ~eff.domhabF, data = umf, keyfun = 'hazard', output = 'density')
fm4 <- fm <- distsamp(~eff.domhabA ~eff.domhabF, data = umf, keyfun = 'uniform', output = 'density')


aictab(list(fm1, fm2, fm3, fm4), second.ord = FALSE)  #create an AIC table to compare the results
