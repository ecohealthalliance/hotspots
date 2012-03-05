# Summary stats

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())      # Clear all variables
graphics.off()       # Close graphics windows

# Load libraries
library(geosphere) # spherical trigonometry functions for geographic applications
library(foreign) # description here

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

dr<-read.table(file="data/eid08_drivers_19OCT11.csv",sep=",",header=TRUE)

# This was the weight used in the original analysis - a better one would be %land area in each grid cell
dr$weight=dr$landarea/mean(dr$landarea)
dr$abslat=abs(dr$lat)

# Total (eventual model goal)
# mt<-glm(anytotr1~lndensity+lnpubs+high_pop_g+mamdiv+precip+abslat, family=binomial,weights=weight, data=dr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic Visualization (linear model)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# view eid data distribution
plot(dr$anytotr1)

# driver column names
drivers <- c("lndensity","lnpubs","high_pop_g","mamdiv","precip","abslat")

# choose driver (1-6)
d <- 1
driver <- drivers[d]

# select response
response <- "anytotr1"

# summary plot with linear model
plot(dr[[driver]], dr[[response]],  
	 ylab = response, xlab = driver)
eid_model <- lm(dr[[response]] ~ dr[[driver]])
abline(eid_model)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Basic Visualization (log model)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eid_logist <- glm(dr[[response]] ~ dr[[driver]], family = binomial)
summary(eid_logist)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plog GLM (log model)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

plot(dr[[driver]], dr[[response]],  
	 ylab = response, xlab = driver)
eid_logist <- glm(dr[[response]] ~ dr[[driver]], family = binomial)
new_driver <- seq(0,8,1)
lines(new_driver,
  predict(eid_logist, data.frame(dr[[driver]] = new_driver), 
  type = "response"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Summary Stats
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eid_logist <- glm(dr[[response]] ~ dr[[driver]], family = binomial)
par(mfrow = c(2,2), mar = c(5,5,1.5,0.5))
plot(eid_logist, which = c(1:4), add.smooth = F, cex.id = 1)
par(mfrow = c(1,1), mar = c(4.5,4.5,0.5,0.5))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Add a variable
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eid_logist.2 <- glm(dr[[response]] ~ dr[[drivers[1]]] + dr[[drivers[2]]], family = binomial)
summary(eid_logist.2)


