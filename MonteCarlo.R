#Set up ----
# Clear environment
rm(list=ls())

# load necessary libraries
library("lmtest")
library("tictoc")

# Simulated LR model ----
# Set up parameters etc.
N <- 200

e <- rnorm(N)
x <- rnorm(N)
a <- 1.3
b <- 0.8
y <- a + b*x + e

par(mfrow=c(1,2))
plot(y)
plot(x)

## Estimation ----
out <- lm(y~x)
summary(out)

## Asymptotic inference ----
# CI for model parameters
confint(out, level=0.95)

# Monte Carlo simulations ----
tic("Simulating a AR(1) model with R = 10000 replications")
## AR(p) models ----
### Example AR(1)----
## The setting:
# We want to build the model
# y(t) = a + b * x(t) + e(t)
#
# However, don't know the distribution of alpha and beta
# but we know that e~N(0,1) [1]
# and we also know that x~N(0,1) [2]
#
# What can you we about alpha and beta?
# 
# Derive the distributional properties
# MC to check your theory
# 
# a_hat~Normal variable
# b_hat~Normal variable
# 
# check if your theoretical derivations
# are correct

# Set the number of Replications
R <- 10000
N <- 200
ar_order <- c(1,0,0)

# Set your sparse matrix
MCc <- matrix(NA, R, 4)
colnames(MCc) <- c("mu_hat", "phi_hat", "Mean_y", "Var_y")

for(i in 1:R)
{
  # First generate your data
  y <- arima.sim(list(order = ar_order, ar=0.5), n=N)
  
  # Estimation
  # AR(1)
  out <- arima(y, order = ar_order)
  coeftest(out)
  out_coeffs <- out$coef
  mu_hat <- out_coeffs[2]
  phi_hat <- out_coeffs[1]
  
  # Saving estimates
  MCc[i,1] <- mu_hat
  MCc[i,2] <- phi_hat
  MCc[i,3] <- mean(y)
  MCc[i,4] <- sd(y)^2

  # For all loops, also add a loop tracker
  cat("Now doing replication ", i, " of ", R, "\n")
}

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MCc[,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MCc[,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(phi), " values")), 
     main=expression(paste(hat(phi), " (coefficient)")), lty=1, lwd=0.5)
plot(MCc[,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MCc[,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)
par(mfrow=c(1,1))

# Check the density plots
par(mfrow=c(2,2))
plot(density(MCc[,1]), main="mu_hat")
plot(density(MCc[,2]), main="phi_hat")
plot(density(MCc[,3]), main="mean of y")
plot(density(MCc[,4]), main="variance of y")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(density(MCc[,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1)
curve(dnorm(x, mean=mean(MCc[,1]), sd=sd(MCc[,1])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1)
curve(dnorm(x, mean=mean(MCc[,2]), sd=sd(MCc[,2])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1)
curve(dnorm(x, mean=mean(MCc[,3]), sd=sd(MCc[,3])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1)
curve(dnorm(x, mean=mean(MCc[,4]), sd=sd(MCc[,4])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")
par(mfrow=c(1,1))

## Inference - confidence bands and t-test
# if you have a linear model of Y(t)=a+b*x(t)+e(t)
# 90%: (5%, 95%)
# 95%: (2.5%, 97.5%)
# 99%: (0.5%, 99.5%)
quantile(MCc[,1], probs=c(2.5/100, 97.5/100))
quantile(MCc[,2], probs=c(2.5/100, 97.5/100))
quantile(MCc[,3], probs=c(2.5/100, 97.5/100))

t.test(y)
toc(log = TRUE)

### Changing the sample size ----
tic("Simulating a AR(1) model and different sample sizes, with R = 10000 replications")
# Set the number of Replications
R <- 10000
N_all <- c(50,100,200,500,700,1000)
ar_order <- c(1,0,0)

# Set your sparse matrix
MC_smplsize <- list()

for(j in 1:length(N_all))
{
  MCc <- matrix(NA, R, 4)
  colnames(MCc) <- c("mu_hat", "phi_hat", "Mean_y", "Var_y")
  
for(i in 1:R)
{
  #Sampling
  y <- arima.sim(list(order = ar_order, ar=0.5), n=N_all[j])
  
  # Estimation
  # AR(1)
  out <- arima(y, order = ar_order)
  coeftest(out)
  out_coeffs <- out$coef
  mu_hat <- out_coeffs[2]
  phi_hat <- out_coeffs[1]
  
  # Saving estimates
  MCc[i,1] <- mu_hat
  MCc[i,2] <- phi_hat
  MCc[i,3] <- mean(y)
  MCc[i,4] <- sd(y)^2

  
  # For all loops, also add a loop tracker
  cat("Now doing replication ", i, " of ", R, " using sample size = ", N_all[j],
      " (", j, "/", length(N_all),")", "\n", sep = "")
}
  
MC_smplsize[[j]] <- MCc
names(MC_smplsize)[j] <- paste("sample_", N_all[j], sep = "")
}

### Plots
color <- c("#005824", "#238b45", "#41ae76", "#66c2a4", "#99d8c9", "#ccece6")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_smplsize[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(phi), " values")), 
     main=expression(paste(hat(phi), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_AR1_densities_diff_smpl.pdf")

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different sample sizes (N)",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     xlim=c(-1,1), ylim=c(0,7))
curve(dnorm(x, mean=mean(MC_smplsize[[1]][,1]), sd=sd(MC_smplsize[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,1]), 
      main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
      add=TRUE)
curve(dnorm(x, mean=mean(MC_smplsize[[i]][,1]), sd=sd(MC_smplsize[[i]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     ylim=c(0,14))
curve(dnorm(x, mean=mean(MC_smplsize[[1]][,2]), sd=sd(MC_smplsize[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,2]), sd=sd(MC_smplsize[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")   
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-1,1), ylim=c(0,7))
curve(dnorm(x, mean=mean(MC_smplsize[[1]][,3]), sd=sd(MC_smplsize[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,3]), sd=sd(MC_smplsize[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     ylim=c(0,6))
curve(dnorm(x, mean=mean(MC_smplsize[[1]][,4]), sd=sd(MC_smplsize[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,4]), sd=sd(MC_smplsize[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off() 
toc(log = TRUE)


### Changing the autoregressive parameters ----
tic("Simulating a AR(1) model wih different AR parameters, with R = 10000 replications")
# Set the number of Replications
R <- 10000
N_all <- c(50,100,200,500,700,1000)
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
ar_order <- c(1,0,0)

# Set your sparse matrix
MC_arcoef <- list()

for(j in 1:length(ar_all))
{
  MCc <- matrix(NA, R, 4)
  colnames(MCc) <- c("mu_hat", "phi_hat", "Mean_y", "Var_y")
  
  for(i in 1:R)
  {
    
    #Sampling
    # N <- 200
    # e <- rnorm(N)
    # x <- rnorm(N)
    # x_lag <- c(NA, head(x, -1))
    y <- arima.sim(list(order = ar_order, ar=ar_all[j]), n=N_all[3])
    
    # Estimation
    # AR(1)
    out <- arima(y, order = ar_order)
    coeftest(out)
    out_coeffs <- out$coef
    mu_hat <- out_coeffs[2]
    phi_hat <- out_coeffs[1]
    
    # Saving estimates
    MCc[i,1] <- mu_hat
    MCc[i,2] <- phi_hat
    MCc[i,3] <- mean(y)
    MCc[i,4] <- sd(y)^2
    MC_arcoef[[j]] <- MCc
    names(MC_arcoef)[j] <- paste("p = ", ar_all[j], sep = "")
    
    # For all loops, also add a loop tracker
    cat("Now doing replication ", i, " of ", R, " AR(p), where p = ", ar_all[j],
        " (", j, "/", length(ar_all),")", "\n", sep = "")
  }
}

### Plots
color <- c("#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824",
           "#238b45", "#41ae76", "#66c2a4", "#99d8c9", "#ccece6")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_arcoef[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(phi), " values")), 
     main=expression(paste(hat(phi), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_AR1_densities_diff_ARq.pdf")
m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different AR(p) parameters",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,12), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,1]), sd=sd(MC_arcoef[[1]][,1])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,1]), sd=sd(MC_arcoef[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#phi
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(-1,1))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,2]), sd=sd(MC_arcoef[[1]][,2])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,2]), 
        main=expression(paste("Density Estimate using MC_arcoef, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,2]), sd=sd(MC_arcoef[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topleft",        
       legend = paste("p = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")   

#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-0.5,0.5), ylim=c(0, 12))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,3]), sd=sd(MC_arcoef[[1]][,3])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,3]), sd=sd(MC_arcoef[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0,8), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,4]), sd=sd(MC_arcoef[[1]][,4])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,4]), sd=sd(MC_arcoef[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

dev.off()
toc(log = TRUE)

### Changing the number of replications ----
tic("Simulating a AR(1) model with different number of replications (R)")
# Set the number of Replications
R_all <- c(50, 100, 200, 500, 700, 1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
ar_order <- c(1,0,0)

# Set your sparse matrix
MC_norepl <- list()

for(j in 1:length(R_all))
{
  MCc <- matrix(NA, R_all[j], 4)
  colnames(MCc) <- c("mu_hat", "phi_hat", "Mean_y", "Var_y")
  
  for(i in 1:R_all[j])
  {
    
    #Sampling
    # N <- 200
    # e <- rnorm(N)
    # x <- rnorm(N)
    # x_lag <- c(NA, head(x, -1))
    y <- arima.sim(list(order = ar_order, ar=ar_all[9]), n=N_all[3])
    
    # Estimation
    # AR(1)
    out <- arima(y, order = ar_order)
    coeftest(out)
    out_coeffs <- out$coef
    mu_hat <- out_coeffs[2]
    phi_hat <- out_coeffs[1]
    
    # Saving estimates
    MCc[i,1] <- mu_hat
    MCc[i,2] <- phi_hat
    MCc[i,3] <- mean(y)
    MCc[i,4] <- sd(y)^2
    MC_norepl[[j]] <- MCc
    names(MC_norepl)[j] <- paste("R = ", R_all[j], sep = "")
    
    # For all loops, also add a loop tracker
    cat("Now doing replication ", i, " of ", R_all[j],
        " (", j, "/", length(R_all), ")",
        "\n", sep = "")
  }
}

### Plots
color <- c("#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824",
           "#238b45")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_norepl[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(phi), " values")), 
     main=expression(paste(hat(phi), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_AR1_densities_diff_norepl.pdf")
m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different number or replications (R)",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,4), xlim=c(-0.4,0.4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,1]), sd=sd(MC_norepl[[1]][,1])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,1]), sd=sd(MC_norepl[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0.3,0.7))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,2]), sd=sd(MC_norepl[[1]][,2])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,2]), sd=sd(MC_norepl[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topleft",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")   

#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-0.3,0.3), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,3]), sd=sd(MC_norepl[[1]][,3])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,3]), sd=sd(MC_norepl[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(1,1.7), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,4]), sd=sd(MC_norepl[[1]][,4])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,4]), sd=sd(MC_norepl[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

dev.off()

## MA(q) models ----
### Example MA(1)----
## The setting:
# We want to build the model
# y(t) = a + b * x(t) + e(t)
#
# However, don't know the distribution of alpha and beta
# but we know that e~N(0,1) [1]
# and we also know that x~N(0,1) [2]
#
# What can you we about alpha and beta?
# 
# Derive the distributional properties
# MC to check your theory
# 
# a_hat~Normal variable
# b_hat~Normal variable
# 
# check if your theoretical derivations
# are correct

# Set the number of Replications
R <- 10000
N <- 200
ma_order <- c(0,0,1)


# Set your sparse matrix
MCc <- matrix(NA, R, 4)
colnames(MCc) <- c("mu_hat", "theta_hat", "Mean_y", "Var_y")

for(i in 1:R)
{
  
  #Sampling
  # N <- 200
  # e <- rnorm(N)
  # x <- rnorm(N)
  # x_lag <- c(NA, head(x, -1))
  
  #Or use arima.sim to simulate an ARIMA model
  y <- arima.sim(list(order = ma_order, ma=0.5), n=N)
  
  # Estimation
  # AR(1)
  out <- arima(y, order=c(0,0,1))
  coeftest(out)
  out_coeffs <- out$coef
  mu_hat <- out_coeffs[2]
  theta_hat <- out_coeffs[1]
  
  # Saving estimates
  MCc[i,1] <- mu_hat
  MCc[i,2] <- theta_hat
  MCc[i,3] <- mean(y)
  MCc[i,4] <- sd(y)^2
  
  # For all loops, also add a loop tracker
  cat("Now doing replication ", i, " of ", R, "\n")
}

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MCc[,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MCc[,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(theta), " values")), 
     main=expression(paste(hat(theta), " (coefficient)")), lty=1, lwd=0.5)
plot(MCc[,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MCc[,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)
par(mfrow=c(1,1))

# Check the density plots
par(mfrow=c(2,2))
plot(density(MCc[,1]), main="mu_hat")
plot(density(MCc[,2]), main="theta_hat")
plot(density(MCc[,3]), main="mean of y")
plot(density(MCc[,4]), main="variance of y")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(density(MCc[,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1)
curve(dnorm(x, mean=mean(MCc[,1]), sd=sd(MCc[,1])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1)
curve(dnorm(x, mean=mean(MCc[,2]), sd=sd(MCc[,2])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1)
curve(dnorm(x, mean=mean(MCc[,3]), sd=sd(MCc[,3])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")

plot(density(MCc[,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1)
curve(dnorm(x, mean=mean(MCc[,4]), sd=sd(MCc[,4])), col="#9ecae1",
      lty=2, lwd=1, add=TRUE, yaxt="n")
par(mfrow=c(1,1))

## Inference - confidence bands and t-test
# if you have a linear model of Y(t)=a+b*x(t)+e(t)
# 90%: (5%, 95%)
# 95%: (2.5%, 97.5%)
# 99%: (0.5%, 99.5%)
quantile(MCc[,1], probs=c(2.5/100, 97.5/100))
quantile(MCc[,2], probs=c(2.5/100, 97.5/100))
quantile(MCc[,3], probs=c(2.5/100, 97.5/100))

t.test(y)
toc(log = TRUE)

### Changing the sample size ----
tic("Simulating a MA(1) model and different sample sizes, with R = 10000 replications")
# Set the number of Replications
R <- 10000
N_all <- c(50,100,200,500,700,1000)
ma_order <- c(0,0,1)

# Set your sparse matrix
MC_smplsize <- list()

for(j in 1:length(N_all))
{
  MCc <- matrix(NA, R, 4)
  colnames(MCc) <- c("mu_hat", "theta_hat", "Mean_y", "Var_y")
  
  for(i in 1:R)
  {
    
    #Sampling
    # N <- 200
    # e <- rnorm(N)
    # x <- rnorm(N)
    # x_lag <- c(NA, head(x, -1))
    y <- arima.sim(list(order = ma_order, ma=0.5), n=N_all[j])
    
    # Estimation
    # AR(1)
    out <- arima(y, order=c(0,0,1))
    coeftest(out)
    out_coeffs <- out$coef
    mu_hat <- out_coeffs[2]
    theta_hat <- out_coeffs[1]
    
    # Saving estimates
    MCc[i,1] <- mu_hat
    MCc[i,2] <- theta_hat
    MCc[i,3] <- mean(y)
    MCc[i,4] <- sd(y)^2
    MC_smplsize[[j]] <- MCc
    names(MC_smplsize)[j] <- paste("sample_", N, sep = "")
    
    # For all loops, also add a loop tracker
    cat("Now doing replication ", i, " of ", R, " using sample size = ", N,
        " (", j, "/", length(N_all),")", "\n", sep = "")
  }
}

### Plots
color <- c("#005824", "#238b45", "#41ae76", "#66c2a4", "#99d8c9", "#ccece6")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_smplsize[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(theta), " values")), 
     main=expression(paste(hat(theta), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_smplsize[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_MA1_densities_diff_smpl.pdf")
m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different sample sizes (N)",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,7), xlim=c(-1,1))
curve(dnorm(x, mean=mean(MC_smplsize[[1]][,1]), sd=sd(MC_smplsize[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,1]), sd=sd(MC_smplsize[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

#theta
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[i]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1,
     xlim=c(0,1))
curve(dnorm(x, mean=mean(MC_smplsize[[i]][,2]), sd=sd(MC_smplsize[[i]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,2]), sd=sd(MC_smplsize[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")   

#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[i]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(MC_smplsize[[i]][,3]), sd=sd(MC_smplsize[[i]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,3]), sd=sd(MC_smplsize[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_smplsize[[i]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0.5,2))
curve(dnorm(x, mean=mean(MC_smplsize[[i]][,4]), sd=sd(MC_smplsize[[i]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_smplsize))
{
  lines(density(MC_smplsize[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_smplsize[[i]][,4]), sd=sd(MC_smplsize[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("N = ", N_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off() 
toc(log = TRUE)

### Changing the autoregressive parameters ----
tic("Simulating a MA(1) model wih different AR parameters, with R = 10000 replications")
# Set the number of Replications
R <- 10000
N_all <- c(50,100,200,500,700,1000)
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
ma_order <- c(0,0,1)

# Set your sparse matrix
MC_arcoef <- list()

for(j in 1:length(ma_all))
{
  MCc <- matrix(NA, R, 4)
  colnames(MCc) <- c("mu_hat", "theta_hat", "Mean_y", "Var_y")
  
  for(i in 1:R)
  {
    
    #Sampling
    # N <- 200
    # e <- rnorm(N)
    # x <- rnorm(N)
    # x_lag <- c(NA, head(x, -1))
    y <- arima.sim(list(order = ma_order, ma=ma_all[j]), n=N_all[3])
    
    # Estimation
    # AR(1)
    out <- arima(y, order=c(0,0,1))
    coeftest(out)
    out_coeffs <- out$coef
    mu_hat <- out_coeffs[2]
    theta_hat <- out_coeffs[1]
    
    # Saving estimates
    MCc[i,1] <- mu_hat
    MCc[i,2] <- theta_hat
    MCc[i,3] <- mean(y)
    MCc[i,4] <- sd(y)^2
    MC_arcoef[[j]] <- MCc
    names(MC_arcoef)[j] <- paste("p = ", ma_all[j], sep = "")
    
    # For all loops, also add a loop tracker
    cat("Now doing replication ", i, " of ", R, " AR(p), where p = ", ma_all[j],
        " (", j, "/", length(ma_all),")", "\n", sep = "")
  }
}

### Plots
color <- c("#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824",
           "#238b45", "#41ae76", "#66c2a4", "#99d8c9", "#ccece6")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_arcoef[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(theta), " values")), 
     main=expression(paste(hat(theta), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_arcoef[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_MA1_densities_diff_ARq.pdf")
m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different MA(q) parameters",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,12), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,1]), sd=sd(MC_arcoef[[1]][,1])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,1]), sd=sd(MC_arcoef[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#theta
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1,
     xlim=c(-1,1))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,2]), sd=sd(MC_arcoef[[1]][,2])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,2]), 
        main=expression(paste("Density Estimate using MC_arcoef, ", hat(theta))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,2]), sd=sd(MC_arcoef[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topleft",        
       legend = paste("p = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")   

#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-0.5,0.5), ylim=c(0, 12))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,3]), sd=sd(MC_arcoef[[1]][,3])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,3]), sd=sd(MC_arcoef[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_arcoef[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0,8), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_arcoef[[1]][,4]), sd=sd(MC_arcoef[[1]][,4])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_arcoef))
{
  lines(density(MC_arcoef[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_arcoef[[i]][,4]), sd=sd(MC_arcoef[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("p = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

dev.off()
toc(log = TRUE)

### Changing the number of replications ----
tic("Simulating a MA(1) model with different number of replications (R)")
# Set the number of Replications
R_all <- c(50, 100, 200, 500, 700, 1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
ma_order <- c(0,0,1)

# Set your sparse matrix
MC_norepl <- list()

for(j in 1:length(R_all))
{
  MCc <- matrix(NA, R_all[j], 4)
  colnames(MCc) <- c("mu_hat", "theta_hat", "Mean_y", "Var_y")
  
  for(i in 1:R_all[j])
  {
    
    #Sampling
    # N <- 200
    # e <- rnorm(N)
    # x <- rnorm(N)
    # x_lag <- c(NA, head(x, -1))
    y <- arima.sim(list(order = ma_order, ma=ma_all[9]), n=N_all[3])
    
    # Estimation
    # AR(1)
    out <- arima(y, order=c(0,0,1))
    coeftest(out)
    out_coeffs <- out$coef
    mu_hat <- out_coeffs[2]
    theta_hat <- out_coeffs[1]
    
    # Saving estimates
    MCc[i,1] <- mu_hat
    MCc[i,2] <- theta_hat
    MCc[i,3] <- mean(y)
    MCc[i,4] <- sd(y)^2
    MC_norepl[[j]] <- MCc
    names(MC_norepl)[j] <- paste("R = ", R_all[j], sep = "")
    
    # For all loops, also add a loop tracker
    cat("Now doing replication ", i, " of ", R_all[j],
        " (", j, "/", length(R_all), ")",
        "\n", sep = "")
  }
}

### Plots
color <- c("#ccece6", "#99d8c9", "#66c2a4", "#41ae76", "#238b45", "#005824",
           "#238b45")

# So, using Monte Carlo
# we are able to "approximate" the distribution of our statistic#
# of interest
par(mfrow=c(2,2), oma=c(0,0,2,0))
plot(MC_norepl[[1]][,1], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(mu), " values")), 
     main=expression(paste(hat(mu), " (intercept)")), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,2], type="l", 
     xlab="Monte Carlo Replication",
     ylab=expression(paste(hat(theta), " values")), 
     main=expression(paste(hat(theta), " (coefficient)")), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,3], type="l", 
     xlab="Monte Carlo Replication",
     ylab="mean values", 
     main=expression("mean of y"), lty=1, lwd=0.5)
plot(MC_norepl[[1]][,4], type="l", 
     xlab="Monte Carlo Replication",
     ylab="variance values", 
     main=expression("variance of y"), lty=1, lwd=0.5)
mtext("Parameters", outer = TRUE, cex = 1.5)

# Check the density plots
pdf("plots/MC_MA1_densities_diff_norepl.pdf")
m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Monte Carlo simulation - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different number or replications (R)",cex=1,font=2, pos = 1)
#mu
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,4), xlim=c(-0.4,0.4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,1]), sd=sd(MC_norepl[[1]][,1])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,1]), sd=sd(MC_norepl[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#theta
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1,
     xlim=c(0.3,0.7))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,2]), sd=sd(MC_norepl[[1]][,2])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(theta))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,2]), sd=sd(MC_norepl[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topleft",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")   

#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     xlim=c(-0.3,0.3), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,3]), sd=sd(MC_norepl[[1]][,3])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,3]), sd=sd(MC_norepl[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(MC_norepl[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(1,1.7), ylim=c(0, 4))
curve(dnorm(x, mean=mean(MC_norepl[[1]][,4]), sd=sd(MC_norepl[[1]][,4])), col=color[1],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(MC_norepl))
{
  lines(density(MC_norepl[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(MC_norepl[[i]][,4]), sd=sd(MC_norepl[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("R = ", R_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.4,
       bg = "transparent",
       bty = "n")

dev.off()
toc(log = TRUE)
log.lst <- tic.log(format = FALSE)
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
total_time <- sum(timings)/60
cat("Simulations finished after", round(total_time, 1), "minutes", "\n")
tic.clearlog()