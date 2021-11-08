# ASSIGNMENT 4
# APPLIED BAYESIAN DATA ANALYSIS
# AUTHOR: MARC VILA
# DATE: 04 NOVEMBER 2021

library(rjags)
library(statip)
library(bayestestR)
library(dplyr)
library(ggplot2)

### 2. 
## a. GIVEN THE FOLLOWING MEASUREMENTS:  y=[1,0,1,1,0,1,1,1,0,1,1,1,1,1]

## i. Expected probability of getting a head?

y <- c(1,0,1,1,0,1,1,1,0,1,1,1,1,1) # coin flips
l <- seq(0.0, 0.0, length.out=length(theta)) # Initializing a vector of the same length as seq of theta values

z <- sum(y)
N <- length(y)
a = 1 # Uniform distribution
b = 1 # Uniform distribution

# DEFINE the model

flip_model <- "model{
    # Likelihood model for z
    z ~ dbin(p, N)
    
    # Prior model for p
    p ~ dbeta(a, b)
}"

# COMPILE the model    
flip_jags <- jags.model(textConnection(flip_model),
                        data = list("a" = a, "b" = b, "z" = z, "N" = N),
                        inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100))

list.samplers(flip_jags) # Sampler method used by jags is slicing

# SIMULATE the posterior
update(flip_jags, n.iter=1000)
flip_sim <- coda.samples(model = flip_jags, variable.names = c("p"), n.iter = 10000)
summary(flip_sim)

#i. What is the expected probability of getting a head?

flip_sim_mat <- as.matrix(flip_sim) # Change MCMC.list into a matrix
flip_sim_dat <- as.data.frame(flip_sim_mat) # Change matrix into a data.frame

mean(flip_sim_mat) # Probability of getting a head

## Give a 95% credible interval of this probability 
## (use either equal-tailed or HDI or both)

ci_hdi <- ci(flip_sim, method = "HDI")
print(ci_hdi)

## What is the probability of theta > 0.5?

sum(flip_sim_mat>0.5)/length(flip_sim_mat)

# PLOT the posterior
mcmc_hist(flip_sim,binwidth = 0.01)
plot(flip_sim, type="l",col="green", trace=FALSE)

## b.Given an additional set of measurements , w = [1,0,0,0,0,0,0,1,1,0]
## are y and w measurements from the same coin? You may answer this in different ways.

w <- c(1,0,0,0,0,0,0,1,1,0) # coin flips
l <- seq(0.0, 0.0, length.out=length(theta))
zw <- sum(w)
Nw <- length(w)

# DEFINE the model

flip_model_w <- "model{
    # Likelihood model for z
    z ~ dbin(p, N)
    
    # Prior model for p
    p ~ dbeta(a, b)
}"

# COMPILE the model    
flip_jags_w <- jags.model(textConnection(flip_model_w),
                        data = list("a" = a, "b" = b, "z" = zw, "N" = Nw),
                        inits = list(.RNG.name = "base::Wichmann-Hill", .RNG.seed = 100))

list.samplers(flip_jags_w) # Sampler method used by jags is slicing

# SIMULATE the posterior
update(flip_jags_w, n.iter=1000)
flip_sim_w <- coda.samples(model = flip_jags_w, variable.names = c("p"), n.iter = 10000)
summary(flip_sim_w)

flip_sim_w_mat <- as.matrix(flip_sim_w) # Change MCMC.list into a matrix
flip_sim_w_dat <- as.data.frame(flip_sim_w_mat) # Change matrix into a data.frame

# b.1. Answer this by calculating the probability that θy > θw given the measurements y,w
# count the number of samples satisfying the condition in relation to the total number of samples

length(which(flip_sim_mat>flip_sim_w_mat))/length(flip_sim_mat)

# b.2 Answer this by creating a new variable θy - θw and calculating a 95% credible interval.

difference <- flip_sim_mat-flip_sim_w_mat
ci_hdi_diff <- ci(difference, method = "HDI")
print(ci_hdi_diff)

# b.3 Plot a histogram representing p(dθ|y,z) . 
# Is it beta distributed? Motivate your answer.


## PLOT 1

set.seed(42)
p1 <- hist(flip_sim_mat,breaks=80)                    
p2 <- hist(flip_sim_w_mat,breaks=100)      
plot( p1, col=rgb(0,0,1,1/4), xlim=c(0,1), freq=FALSE, main="Sampling of w and y",
      xlab="Theta")  # first histogram
abline(v = mean(flip_sim_w_mat), col="red", lwd=3, lty=2)
abline(v = mean(flip_sim_mat), col="blue", lwd=3, lty=2)
abline(v = 0.54, col="light blue", lwd=3, lty=2)
abline(v = 0.93, col="light blue", lwd=3, lty=2)
text(0.538,4.0,cex = 0.88, "0.538",pos = 2, col = "lightblue")
text(0.932,4.0,cex = 0.88, "0.932",pos = 2, col = "lightblue")
text(0.331,4.0,cex = 0.88, "0.331",pos = 2, col = "red")
text(0.749,4.0,cex = 0.88, "0.749",pos = 2, col = "blue")
plot( p2, col=rgb(1,0,0,1/4), xlim=c(0,1),freq=FALSE, add=TRUE)  # second


## PLOT 2

p3 <- hist(difference,breaks=100)      
plot( p3, col=rgb(0,0,1,1/4), xlim=c(0,1), freq=FALSE, main="dTheta sampling",
      xlab="dTheta(thetay - thetaw)")
abline(v = mean(difference), col="red", lwd=3, lty=2)
abline(v = 0.083, col="blue", lwd=3, lty=2)
abline(v = 0.727, col="blue", lwd=3, lty=2)
text(0.083,2.5,cex = 0.88, "0.083",pos = 2, col = "blue")
text(0.727,2.5,cex = 0.88, "0.727",pos = 2, col = "blue")
text(0.417,2.5,cex = 0.88, "0.417",pos = 2, col = "red")

