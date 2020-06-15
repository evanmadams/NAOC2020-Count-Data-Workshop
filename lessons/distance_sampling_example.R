#NAOC 2020 count data workshop distance sampling example
#By Evan Adams and Beth Ross
###########################################################


#we want to survey the different ways that count data become biased and this is an example that occurs frequently with bird survey data
#when you survey an area you are often better at survey the area closest to you rather than further away
#distance sampling was created to try and understand how the chances of detecting a bird vary with their distance from the observer
#So we've added an additional source of information to the count data, for every count we have a distance to the observer
#We use the distance data to make another model that we combine with a simple Poisson regression model so that we can describe the patterns of abundance after accounting for detection bias
#once we estimate this source of bias, then we can account for the birds that we didn't see and correct our estimates

#we are using three packages today, unmarked is the package that does the distance sampling, AICcmodavg gives us some additional help with unmarked, reshape2 and dplyr help format data, and ggplot2 is essentially useful for all R projects
library(unmarked)
library(AICcmodavg)
library(reshape2)
library(dplyr)
library(ggplot2)

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

#pull in covariates that describe either abundance or detection probability
#these include habitat data and conditions under which the survey was conducted 
#Note that the wind speed, sky condition and background noise are categorical and the date/time data are scaled continuous

sc <- data.frame(dat$habcovs, dat$distcovs, dat$avcovs)


#then we need to convert it to a format that unmarked understands
#unmarked has a function that helps us do that specifically for distance sampling data
#I'll just be grabbing the first 100 sites so that the models finish quickly but you can explore the full data set if you wish

umf <- unmarkedFrameDS(y = as.matrix(y[1:100 ,2:6]), siteCovs = sc[1:100,], dist.breaks = c(0, 20, 40, 60, 80, 100), survey = 'point', unitsIn = 'm')   #parse the data into the unmarked data frame and subset the data to (1) have models converge faster and (2) remove birds with unknown distances

#let's have a look at these data to get a sense of what kinds of models will be useful to build with them

#first, we should have a look at the count data as that could give us some insight into what kind of error distribution to select

mean(rowSums(umf@y))

#changing the plot settings to compare two figures
par(mfrow=c(1,2))

#plotting the histogram of the actual count data
hist(rowSums(umf@y))

#does a poisson distribution with the same mean look similar?
hist(rpois(nrow(umf@y), mean(rowSums(umf@y))))


#now we want to run the distance model using 'distsamp'
#the first element of the function needs a model formula that is similar to what we were using only we need multiple formulas for the distance model and the ecological model
#here we are running a null model with covariates on either submodel
#we also designate how we want to model the detection probability over distance (a half-normal distribution in this case)
fm <- distsamp(~1 ~1, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm)
plot(fm)
plot(ranef(fm))

#LOOK AT OVERALL MODEL FIT
#one of the best ways to do this is to use bootstrapping to determine how our model fit compares to all the boostrapped fits
fit <- Nmix.gof.test(fm, nsim = 100)

#in our case the observed model (the red line) is not so different from all the bootstrapped cases (hence the P > 0.05) and that overall model fit is likely good

#also can look at something called c-hat to see if there is evidence that the count data are overdispersed to the modeled error estiamtes
print(fit)

#c-hat values close to 1 indicate that overdispersion is fairly minimal that our model fits the data reasonable well


#DISTANCE MODEL FIT

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

#what if detection probability was influenced by wind speed at the site?

fm <- distsamp(~Wind.Speed.Code ~1, data = umf, keyfun = 'halfnorm', output = 'density')

summary(fm)

#not much evidence for that in this data set, but we perhaps it's an effect that we need more data to estimate properly

#what if the presence of developed habitat at the site?

fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'halfnorm', output = 'density', unitsOut = 'ha')

summary(fm)

#Well that's pretty strong evidence that CHSP densities are increasing in developed habitat

#let's try a more complicated model that suggests that detections are a function of wind speed, background noise, and time of year; and abundance is forest and developed habitat

fm <- distsamp(~Wind.Speed.Code + Background.Noise.Code + ydate ~eff.domhabA + eff.domhabF, data = umf, keyfun = 'halfnorm', output = 'density', unitsOut = 'ha')

summary(fm)

#let's check the fit of this model

fit <- Nmix.gof.test(fm, nsim = 100)

print(fit)

#even better than before, so that suggests that the covariates we added were perhaps useful

#now what do we need to visualize these effects?
#to do this, we want to use the predict function in unmarked to show these effects
#so first let's create some new data that we want to predict to

#for abundance we want to predict for different combinations of the dominant habitat types

newdat <- data.frame(eff.domhabF = c(0, 0, 1, 1), eff.domhabA = c(0, 1, 0, 1))

ps <- predict(fm, newdat, type = 'state')

#so we have predicted densities and 95% CIs for each combination of mixed forest and developed habitat
#let's add a couple thigns to the data frame and then plot the differences

ps <- data.frame(ps, newdat)
ps$Habitat <- c('None', 'Developed', 'Mixed Forest', 'Both')


p <- ggplot(ps, aes(x = Habitat, y = Predicted, ymin = lower, ymax = upper))
p <- p + geom_point() + geom_errorbar(aes(width = 0.2))
p <- p + theme_bw() + ylab('Predicted Density (birds/ha)')

plot(p)

#so now we can see how both developed and mixed forest habitat contributes to the density of CHSP


#we can also visualize how covariates influence detection probability

#effect plots for detection
#This time we will make predictions for each factor while holding the others constant. So we vary one and use the mean for the others

#wind speed

newdat <- data.frame(Wind.Speed.Code = c(0, 1, 2, 3, 4), Background.Noise.Code = mean(sc$Background.Noise.Code), ydate = 0)

ps <- predict(fm, newdat, type = 'det')

ps <- data.frame(ps, newdat)

#We estimated sigma but we still need to conver this to detection probability, so we do that by using a similar process to above
#only I'll use lapply to speed up the process across multiple predictions

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = Wind.Speed.Code, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Wind Speed (Beaufort scale)')

plot(p)

#Background noise

newdat <- data.frame(Wind.Speed.Code = mean(sc$Wind.Speed.Code), Background.Noise.Code = c(0, 1, 2), ydate = 0)

ps <- predict(fm, newdat, type = 'det')

ps <- data.frame(ps, newdat)

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = Background.Noise.Code, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Background Noise (categorical)')

plot(p)

#Date

newdat <- data.frame(Wind.Speed.Code = mean(sc$Wind.Speed.Code), Background.Noise.Code = mean(sc$Background.Noise.Code), ydate = seq(from = min(sc$ydate), to = max(sc$ydate), by = .1))

ps <- predict(fm, newdat, type = 'det')

ps <- data.frame(ps, newdat)

ps$det <- unlist(lapply(ps$Predicted, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.lower <- unlist(lapply(ps$lower, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))
ps$det.upper <- unlist(lapply(ps$upper, FUN = function(x) integrate(gxhn, 0, 100, sigma = x)$value/100))

p <- ggplot(ps, aes(x = ydate, y = det, ymin = det.lower, ymax = det.upper))
p <- p + geom_point() + geom_line() + geom_ribbon(alpha = 0.2)
p <- p + theme_bw() + ylab('Detection Probability') + xlab('Date (scaled)')

plot(p)

#Okay, now you've seen the basics of what you can do with distance sampling in unmarked, let's test some of your new skills

####CHALLENGE####

#Often we are interested in finding the detection function that fits the data that we have the best
#Try your hand at doing this. Make a few models to compare different detection functions

#note please use 'aictab' in the AICcmodavg library to compare multiple models

#check out ?distsamp for the other types of detectin functions in can fit


##answer

fm1 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'halfnorm', output = 'density')
fm2 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'exp', output = 'density')
fm3 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'hazard', output = 'density')
fm4 <- fm <- distsamp(~1 ~eff.domhabA, data = umf, keyfun = 'uniform', output = 'density')


aictab(list(fm1, fm2, fm3, fm4), second.ord = FALSE)  #create an AIC table to compare the results

#compare the fits of the various models

par(mfrow=c(2,2))
hist(fm1, xlab = 'Distance (m)')
hist(fm2, xlab = 'Distance (m)')
hist(fm3, xlab = 'Distance (m)')
hist(fm4, xlab = 'Distance (m)')

                            
###END LESSON###
