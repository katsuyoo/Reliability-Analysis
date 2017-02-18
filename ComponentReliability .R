# Reliability Analysis - BETA

# For a good tutorial in survival (reliability) analysis using R, go to:
# http://www.openintro.org/stat/down/Survival-Analysis-in-R.pdf.

library(survival)

# PARAMETERS
n <- 35 # sample size
r <- 35  # stopping number of failures

confidence <- 0.95

# READ THE DATA
data <- read.csv("K25-Reliability.csv", header=TRUE, sep=",")
print(data)

# NONPARAMETRIC ANALYSIS

# Create a survival object where 
# 1. time = time of failure (if failed) or censoring (if survived)
# for right-censored data
# 2. event = 1 if failed, 0 if survived
#
# For right censoring, set the time variable as the 
survObj <- Surv(time=data$Days, event=data$Censored, type='right')
print(survObj)

# Estimate Reliability (R) using Kaplan-Meier nonparametric analysis.
# 
# The formula takes the data versus an intercept of 1
# (i.e., 100% reliabile at t=0.
#
# For the conf.type, use "log" if exponential and "log-log" if Weibull. 
kp <- survfit(formula=survObj ~ 1, 
                 conf.int=confidence, 
                 conf.type="log-log")
print(kp)
fit <- summary(kp)
print(fit)
t <- fit$time # time
R.hat <- fit$surv # reliability
Lamda.hat <- -log(R.hat) # cumulative hazard
print(Lamda.hat)
lamda.hat <- (Lamda.hat[2:r] - Lamda.hat[1:(r-1)]) / (t[2:r] - t[1:(r-1)]) # hazard rate
print(lamda.hat)
f.hat <- lamda.hat * R.hat[1:(r-1)] # probability density
print(f.hat)

# Plot the Kaplan-Meier nonparametric curves:  
# Reliability, P(Failure), Hazard Rate, Cumulative Hazard Rate
par(mfrow=c(2,2), oma=c(0,0,2,0)) # put all four plots in a 2x2 layout

plot(kp, 
     main="Reliability (95% Confidence)", 
     xlab="Days", 
     ylab="Reliability")

plot(x=t[1:(r-1)], 
     y=f.hat, 
     main="Probability Density", 
     xlab="Days", 
     ylab="Probability Density",
     xlim=c(0,t[r]), 
     cex=0.25)

# Calculate and plot the cumulative hazard function.
plot(kp, 
     main="Cumulative Hazard (95% Confidence)", 
     xlab="Days", 
     ylab="Cumulative Hazard",
     fun="cumhaz")

# Calculate and plot the hazard rate
plot(x=t[1:(r-1)], 
     y=lamda.hat,
     main="Hazard Rate", 
     xlab="Days", 
     ylab="Hazard Rate",
     xlim=c(0,t[r]), 
     cex=0.25)

# Print the mean time to failure (rmean=MTTF), the stanard error of the mean, 
# the median time to failure, and the confidence interval of the median.
print(R.hat, print.rmean=TRUE)

title("Kaplan-Meier", outer=TRUE)

# PARAMETRIC ANALYSIS

# Fit the exponential distribution.
eFit <- survreg(survObj~1, dist="exponential")
summary(eFit)

lamda <- exp(-eFit$coeff)
print(lamda)

stderr <- summary(eFit)$table[,2] # standard error
lamda.lcl <- exp(-(eFit$coeff+2*stderr))
print(lamda.lcl)
lamda.ucl <- exp(-(eFit$coeff-2*stderr))
print(lamda.ucl)

MTTF <- 1 / lamda
print(MTTF)

# Plot the exponential distribution's curves:  
# Reliability, P(Failure), Hazard Rate, Cumulative Hazard Rate
par(mfrow=c(2,2), oma=c(0,0,2,0)) # put all four plots in a 2x2 layout

x<-c(0, 1:100/100) # quantiles (centiles) to use in curve function

curve(exp(-lamda*x), 
      from=0, 
      to=max(data$Days), 
      ylim=c(0,1),
      col="green", 
      main="Reliability",
      ylab="Reliability", 
      xlab='Days')
points(x=t, y=R.hat, cex=0.25)

curve((lamda*exp(-lamda*x)), 
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="P(Failure)",
      ylab="P(Failure)", 
      xlab='Days')
points(x=t[1:(r-1)], y=f.hat, cex=0.25)

curve((lamda*x), 
      from=0, 
      to=max(data$Days), 
      col="darkred", 
      main="Cumulative Hazard Rate",
      ylab="Cumulative Hazard Rate", 
      xlab='Days')
points(x=t, y=Lamda.hat, cex=0.25)

curve(lamda*x/x, # need x/x for function to work
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="Hazard Rate",
      ylab="Hazard Rate", 
      xlab='Days')
points(x=t[1:(r-1)], y=lamda.hat, cex=0.25)

title("Exponential Distribution", outer=TRUE)

# Fit the Weibull distribution.
# 
# Unfortunately, the survreg function provides non-intuitive results:
# -- slope = 1 / 'scale' = 1 / exp('log(scale)')
# -- the intercept is stored as 'coefficients[1]'
wFit <- survreg(survObj~1, dist="weibull")
summary(wFit)
summaryTable <- summary(wFit)$table

slope <- 1 / wFit$scale
intercept <- wFit$coefficients[1]

delta.hat <- slope # shape
print(delta.hat)
delta.hat.lcl <- 1 / exp(summaryTable[2,1] + 2*summaryTable[2,2])
print(delta.hat.lcl)
delta.hat.ucl <- 1 / exp(summaryTable[2,1] - 2*summaryTable[2,2])
print(delta.hat.ucl)

theta.hat <- exp(intercept) # scale
print(theta.hat)
theta.hat.lcl <- exp(intercept - 2*summaryTable[1,2])
print(theta.hat.lcl)
theta.hat.ucl <- exp(intercept + 2*summaryTable[1,2])
print(theta.hat.ucl)

MTTF <- theta.hat * gamma(1 + (1/delta.hat))
print(MTTF)

# Plot the Weibull distribution's curves:  
# Reliability, P(Failure), Hazard Rate, Cumulative Hazard Rate
x<-c(0, 1:100/100) # quantiles (centiles) to use

par(mfrow=c(2,2), oma=c(0,0,2,0)) # put all four plots in a 2x2 layout

curve(exp(-((x/theta.hat)^(delta.hat))), 
      from=0, 
      to=max(data$Days), 
      ylim=c(0,1),
      col="green", 
      main="Reliability",
      ylab="Reliability", 
      xlab='Days')
points(x=t, y=R.hat, cex=0.25)

curve((delta.hat/theta.hat)*(x/theta.hat)^(delta.hat-1)*exp(-((x/theta.hat)^(delta.hat))), 
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="P(Failure)",
      ylab="P(Failure)", 
      xlab='Days')
points(x=t[1:(r-1)], y=f.hat, cex=0.25)

curve((x/theta.hat)^delta.hat, 
      from=0, 
      to=max(data$Days), 
      col="darkred", 
      main="Cumulative Hazard Rate",
      ylab="Cumulative Hazard Rate", 
      xlab='Days')
points(x=t, y=Lamda.hat, cex=0.25)

curve((delta.hat/theta.hat)*(x/theta.hat)^(delta.hat-1), 
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="Hazard Rate",
      ylab="Hazard Rate", 
      xlab='Days')
points(x=t[1:(r-1)], y=lamda.hat, cex=0.25)

title("Weibull Distribution", outer=TRUE)

# Fit the lognormal distribution.
lnFit <- survreg(survObj~1, dist="lognormal")
summary(lnFit)
summaryTable <- summary(lnFit)$table
print(summaryTable)
slope <- lnFit$scale
intercept <- summaryTable[1,1]

mu.lnt <- intercept
print(mu.lnt)
mu.lnt.lcl <- intercept - 2*summaryTable[1,2]
print(mu.lnt.lcl)
mu.lnt.ucl <- intercept + 2*summaryTable[1,2]
print(mu.lnt.ucl)

#sigma.lnt <- exp(summaryTable[2,1])
sigma.lnt <- slope
print(sigma.lnt)
sigma.lnt.lcl <- exp(summaryTable[2,1] - 2*summaryTable[2,2])
print(sigma.lnt.lcl)
sigma.lnt.ucl <- exp(summaryTable[2,1] + 2*summaryTable[2,2])
print(sigma.lnt.ucl)

MTTF <- exp(mu.lnt + 0.5*(sigma.lnt^2))
print(MTTF)

# Plot the lognormal distribution's curves:  
# Reliability, P(Failure), Hazard Rate, Cumulative Hazard Rate
x<-c(0, 1:100/100) # quantiles (centiles) to use

par(mfrow=c(2,2), oma=c(0,0,2,0)) # put all four plots in a 2x2 layout

curve(1-pnorm((log(x)-mu.lnt)/sigma.lnt), 
      from=0, 
      to=max(data$Days), 
      ylim=c(0,1),
      col="green", 
      main="Reliability",
      ylab="Reliability", 
      xlab='Days')
points(x=t, y=R.hat, cex=0.25)

curve(exp(-0.5*(((log(x)-mu.lnt)/sigma.lnt)^2)) / (sqrt(2*pi) * sigma.lnt * x), 
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="P(Failure)",
      ylab="P(Failure)", 
      xlab='Days')
points(x=t[1:(r-1)], y=f.hat, cex=0.25)

curve(-log(1-pnorm((log(x)-mu.lnt)/sigma.lnt)), 
      from=0, 
      to=max(data$Days), 
      col="darkred", 
      main="Cumulative Hazard Rate",
      ylab="Cumulative Hazard Rate", 
      xlab='Days')
points(x=t, y=Lamda.hat, cex=0.25)

curve((exp(-.5*(((log(x)-mu.lnt)/sigma.lnt)^2)) / (sqrt(2*pi) * sigma.lnt * x)) / (1-pnorm((log(x)-mu.lnt)/sigma.lnt)), 
      from=0, 
      to=max(data$Days), 
      col="red", 
      main="Hazard Rate",
      ylab="Hazard Rate", 
      xlab='Days')
points(x=t[1:(r-1)], y=lamda.hat, cex=0.25)

title("Lognormal Distribution", outer=TRUE)

# Reset plotting back to one chart at a time.
par(mfrow=c(1,1))

