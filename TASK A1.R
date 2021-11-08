# ASSIGNMENT 4
# APPLIED BAYESIAN DATA ANALYSIS
# AUTHOR: MARC VILA
# DATE: 04 NOVEMBER 2021

library(rjags)
#library(R2jags) #  interface between R and JAGS
#library(runjags)
library(coda) #general tools for analyzing and graphing MCMC algorithms
library(statip)
library(bayestestR)
library(dplyr)
library(bayesplot) #a number of useful plots using ggplot2
library(ggplot2)

# VARIABLES FROM FIGURE 6.4 OF THE BOOK

theta <- seq(0, 1, length.out = 100000)  # Generate sequence of theta values

z = 17
N = 20
a <- c(100) #(18.25),(100)
b <- c(100) #(6.75),(100)

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
list.samplers(flip_jags)

# SIMULATE the posterior
update(flip_jags, n.iter=1000)
flip_sim <- coda.samples(model = flip_jags, variable.names = c("p"), n.iter = 10000)
summary(flip_sim)

# PLOT the posterior

plot_posterior <- hist(as.matrix(flip_sim),breaks=100)      
plot( plot_posterior, col=rgb(0,0,1,1/4), xlim=c(0,1), freq=FALSE, main="Posterior sampling with a=100 and b=100",
      xlab="Theta")

