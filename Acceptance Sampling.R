# Script for Acceptance Sampling.

# R has a lot of base capabilities, but we also need some 
# special libraries.
library(AcceptanceSampling) # Acceptance Sampling

# Given
N = 1000 # lot size
p1 = 0.01 # acceptable quality number (AOQ)
p2 = 0.05 # lot tolerance percent defective (LTPD) > p1
alpha = 0.05 # false positive ratio
beta = 0.1 # false negative ratio

# Find a single sampling plan (equivalent to a nomograph lookup).
#
# PRP is the 'producer risk plan' with a probability of accepting
# at least 1-alpha (i.e., rejecting alpha) with an AOQ of p1.
# CRP is the 'consumer risk plan' with a probability of accepting 
# at leat beta with a LTPD of p2.
samplingPlan <- find.plan(PRP=c(p1, (1-alpha)), CRP=c(p2, beta), type="binom")
n <- samplingPlan$n # sample size
print(n)
c <- samplingPlan$c # acceptance number
print(c)
r <- samplingPlan$r # rejection number
print(r) # should be c + 1 for single sampling

# Get the OC curve for the single sampling plan.
oc <- OC2c(n=n, c=c, type="b") # binomial
plot(oc, xlim=c(0.0, 0.1)) # plot with x-axis p=[0.0, 0.1]

# Get the OC curve for a double sampling plan.  We can experiment with
# different n1, n2, c1, c2, r1 and r2.
n1 <- n # stage 1 sample size
n2 <- n # stage 2
c1 <- c / 2 # stage 1 acceptance number
c2 <- c # stage 2
r1 <- r # stage 1 rejection number
r2 <- r # stage 2
oc2 <- OC2c(n=c(n, n2), c=c(c, c2), r=c(r2, r2), type="b")
plot(oc2, xlim=c(0.0, 0.1)) # plot with x-axis p=[0.0, 0.1]
