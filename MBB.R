#Set up ----
# Clear environment
rm(list=ls())

## Load necessary libraries
library("tseries")
library("blocklength")
library("lmtest")
library("tictoc")

# Load a custom-made MBB function
source("functions/MovingBlockBoostrap.R")

#time it
tic("Counting running time for script")

# AR(1) example with US CPI ----
## Data ----
### Raw ----
# Load the CPI data
df <- read.csv("data/CPIAUCSL_full.csv", header=TRUE)

#Cut the sample
df <- df[400:NROW(df),]
d <- as.Date(df[,1])
data <- as.numeric(df[,2])

# Plot the data
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(d, data, type="l", main="US CPI")
acf(data, main = "ACF - US CPI")
pacf(data, main = "PACF - US CPI")

### Stationary ----
# Make data stationary
# period-to-period inflation (log growth of CPI)
y <- diff(log(data))*100
d <- d[2:NROW(d)]

#Save raw series
y_raw <- y
d_raw <- d

# Plot again the data
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
plot(d, y, type="l", main="US CPI")
acf(y, main = "ACF - US CPI")
pacf(y, main = "PACF - US CPI")

t.test(y) #this suggests that the true mean is equal to 0 
# (we don't reject the null H_0: true mean is equal to 0) 
jarque.bera.test(y) #this suggests that our stationary data is not normal
# (we reject the null hypothesis of normality)

# The below sections use y as the US CPI (transformed as above)
## MBB ----

# Set the bootstrap draws
B <- 10000

#Dataset and sample size
y <- y_raw
N <- NROW(y)

### Model selection for y ----
# Here y is the log growth of US CPI
# Say that we model our inflation as an AR(1) 
#(see PACF of see at first lag)
out <- arima(y, order=c(1,0,0))
coeftest(out)

### Selection of block size ----
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
b1 <- round(N^(1/3))
b2 <- round(N^(1/4))
b3 <- round(N^(1/5))
b4 <- pwsd(y, correlogram = FALSE)

# Initialize matrix of resamples
boot_est <- matrix(NA, B, 3)
colnames(boot_est) <- c("mean", "ar", "intercept")

### Bootstrapping ----
for(i in 1:B)
{
  #Resample via MBB
  yb_mbb <- MBB(y, b=b1)
  
  #Re-calculate mean of resample
  boot_est[i,1] <- mean(yb_mbb)
  
  #Re-estimate and AR(1) for y
  out_boot <- arima(yb_mbb, order=c(1,0,0))
  coef_boot <- out_boot$coef
  boot_est[i,2] <- coef_boot[1]
  boot_est[i,3] <- coef_boot[2]
  
  cat("Now doing", i, "of", B, "bootstrap draws", "\n")
}

# what about the properties of dependence?
# par(mfrow=c(1,3))
# acf(y)
# acf(yb_iid)
# acf(yb_mbb)
# par(mfrow=c(1,1))

### Plots ----
m <- matrix(c(1, 1, 2, 3, 4, 5, 6, 7), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.10,0.3,0.3,0.3))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - AR(1) model for US CPI",cex=1.5,font=2, pos = 1)
# u <- par("usr")
# text(1,u[4],labels = "Here",col = "black",pos = 1)

par(mar = c(3,1,1,1)+1) 
plot(boot_est[,1], type="l", main = "Mean across draws", 
     ylab = "Value", xlab = "No. of bootstrap draw")
plot(density(boot_est[,1]), main = "Density of mean")
     curve(dnorm(x, mean=mean(boot_est[,1]),
            sd=sd(boot_est[,1])),
      col="#9ecae1", lty=2, lwd=2,
      add=TRUE, yaxt="n")

par(mar = c(3,1,1,1)+1) 
plot(boot_est[,2], type="l", main = "AR coefficient across draws", 
    ylab = "Value", xlab = "No. of bootstrap draw")
plot(density(boot_est[,2]), main = "Density of AR coefficient")
curve(dnorm(x, mean=mean(boot_est[,2]),
           sd=sd(boot_est[,2])),
     col="#9ecae1", lty=2, lwd=2,
     add=TRUE, yaxt="n")

par(mar = c(3,1,1,1)+1) 
plot(boot_est[,3], type="l", main = "AR intercept across draws", 
     ylab = "Value", xlab = "No. of bootstrap draw")
plot(density(boot_est[,3]), main = "Density of AR intercept")
curve(dnorm(x, mean=mean(boot_est[,3]),
            sd=sd(boot_est[,3])),
      col="#9ecae1", lty=2, lwd=2,
      add=TRUE, yaxt="n")


### Confidence Intervals ----
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))

# AR(p) models ----
## AR(1): changing the sample size of y ----
cat("Now doing changing sample size of y in AR(1) model", "\n")
Sys.sleep(2)
# Settings
ar_order <- c(1,0,0)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in N_all){
jj <- jj + 1

#Simulate y series 
y <- arima.sim(list(order = ar_order, ar=ar_all[9]), n=N_all[jj])

#Sample size
N <- NROW(y) 

#Block size
lambda1 <- round(N^(1/3))

# Set the bootstrap draws
B <- 10000

# Initialize matrix of resamples
boot_est <- matrix(NA, B, 4)
colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
### Bootstrapping
for(i in 1:B)
{
  #Resample via MBB
  yb_mbb <- MBB(y, b=lambda1)

  #Re-estimate and AR(1) for y
  out_boot <- arima(yb_mbb, order=c(1,0,0))
  coef_boot <- out_boot$coef
  mu_hat <- coef_boot[2]
  phi_hat <- coef_boot[1]
  
  # Saving estimates
  boot_est[i,1] <- mu_hat
  boot_est[i,2] <- phi_hat
  boot_est[i,3] <- mean(yb_mbb)
  boot_est[i,4] <- sd(yb_mbb)^2
  
  # For all loops, also add a loop tracker
  cat("Now doing ", i, " of ", B, " bootstrap draws, using sample size N = ", N,
      " (", jj, "/", NROW(N_all), ")", sep = "", "\n")
}
BB[[jj]] <- boot_est
names(BB)[jj] <- paste("sample_", N, sep = "")
}

### Plots
color <- c("#99d8c9","#66c2a4","#41ae76",
           "#238b45","#006d2c","#00441b")

# Check the density plots
pdf("plots/MBB_AR1_densities_diff_smpl.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different sample sizes (N)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
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
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0, 0.8), ylim=c(0,14))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
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
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
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
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0, 3), ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
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

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))

## AR(1): changing the autoregressive coefficient ----
cat("Now doing changing the autoregressive coefficient in AR(1) model", "\n")
Sys.sleep(2)

# Settings
ar_order <- c(1,0,0)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in ar_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ar_order, ar=ar_all[jj]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- round(N^(1/3))
  
  # Set the bootstrap draws
  B <- 10000
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and AR(1) for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws, using AR coefficient phi = ", ar_all[jj],
        " (", jj, "/", NROW(ar_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("phi =", ar_all[jj], sep = "")
}

### Plots
color <- c("#f7fcf5","#e5f5e0","#c7e9c0","#a1d99b",
  "#74c476","#41ab5d","#238b45","#006d2c",
  "#00441b", "#252525", "#000000")

# Check the density plots
pdf("plots/MBB_AR1_densities_diff_ARq.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different AR coefficients(phi)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,9), xlim=c(-0.5, 0.5))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(-1, 1), ylim=c(0,12))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,9), xlim=c(-0.5, 0.5))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0, 6), ylim=c(0,6))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ar_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))


## AR(1): changing the number of draws B ----
cat("Now doing changing the number of draws in AR(1) model", "\n")
Sys.sleep(2)

# Settings
ar_order <- c(1,0,0)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in B_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ar_order, ar=ar_all[9]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- round(N^(1/3))
  
  # Set the bootstrap draws
  B <- B_all[jj]
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and AR(1) for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws",
        " (", jj, "/", NROW(B_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("B = ", B_all[jj], sep = "")
}

### Plots
color <- c("#edf8fb","#ccece6","#99d8c9",
           "#66c2a4","#41ae76","#238b45","#005824")

# Check the density plots
pdf("plots/MBB_AR1_densities_diff_drawsB.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different number of draws (B)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,6), xlim=c(-0.5,1))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0.1, 0.7), ylim=c(0,8))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,6), xlim=c(-0.5,1))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0.5, 2), ylim=c(0,4))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))


## AR(1): changing the block length (lambda) ----
cat("Now doing changing the block length in AR(1) model", "\n")
Sys.sleep(2)

# Settings
ar_order <- c(1,0,0)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
N_lambda <- N_all[3]
ar_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_lambda^(1/3)), round(N_lambda^(1/4)), round(N_lambda^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in lambda_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ar_order, ar=ar_all[9]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- lambda_all[jj]
  
  # Set the bootstrap draws
  B <- B_all[7]
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and AR(1) for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws, using block size l = ",
        lambda_all[jj], " (", jj, "/", NROW(lambda_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("B = ", B_all[jj], sep = "")
}

### Plots
color <- c("#41ae76","#238b45","#005824")

# Check the density plots
pdf("plots/MBB_AR1_densities_diff_blocklength.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated AR(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different block size (lambda)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,4), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0.1, 0.7), ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,4), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0.5, 2), ylim=c(0,4))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))

# MA(q) models ----
## MA(1): changing the sample size of y ----
cat("Now doing changing sample size of y in MA(1) model", "\n")
Sys.sleep(2)

# Settings
ma_order <- c(0,0,1)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in N_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ma_order, ma=ma_all[9]), n=N_all[jj])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- round(N^(1/3))
  
  # Set the bootstrap draws
  B <- 10000
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and  for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws, using sample size N = ", N,
        " (", jj, "/", NROW(N_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("sample_", N, sep = "")
}

### Plots
color <- c("#99d8c9","#66c2a4","#41ae76",
           "#238b45","#006d2c","#00441b")

# Check the density plots
pdf("plots/MBB_MA1_densities_diff_smpl.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different sample sizes (N)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,9))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
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
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0, 0.8), ylim=c(0,15))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
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
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,9))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
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
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0, 3), ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
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

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))

## MA(1): changing the autoregressive coefficient ----
cat("Now doing changing the autoregressive coefficient in MA(1) model", "\n")
Sys.sleep(2)

# Settings
ma_order <- c(0,0,1)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in ma_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ma_order, ma=ma_all[jj]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- round(N^(1/3))
  
  # Set the bootstrap draws
  B <- 10000
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and  for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws, using AR coefficient phi = ", ma_all[jj],
        " (", jj, "/", NROW(ma_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("phi =", ma_all[jj], sep = "")
}

### Plots
color <- c("#f7fcf5","#e5f5e0","#c7e9c0","#a1d99b",
           "#74c476","#41ab5d","#238b45","#006d2c",
           "#00441b", "#252525", "#000000")

# Check the density plots
pdf("plots/MBB_MA1_densities_diff_ARq.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different AR coefficients(phi)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,10), xlim=c(-0.5, 0.5))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(-1, 1), ylim=c(0,9))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,10), xlim=c(-0.5, 0.5))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0, 6), ylim=c(0,6))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("phi = ", ma_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))


## MA(1): changing the number of draws B ----
cat("Now doing changing the number of draws in MA(1) model", "\n")
Sys.sleep(2)

# Settings
ma_order <- c(0,0,1)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_all^(1/3)), round(N_all^(1/4)), round(N_all^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in B_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ma_order, ma=ma_all[9]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- round(N^(1/3))
  
  # Set the bootstrap draws
  B <- B_all[jj]
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and MA(1) for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws",
        " (", jj, "/", NROW(B_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("B = ", B_all[jj], sep = "")
}

### Plots
color <- c("#edf8fb","#ccece6","#99d8c9",
           "#66c2a4","#41ae76","#238b45","#005824")

# Check the density plots
pdf("plots/MBB_MA1_densities_diff_drawsB.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different number of draws (B)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,5), xlim=c(-0.5,1))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0.1, 0.7), ylim=c(0,8))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,5), xlim=c(-0.5,1))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0.5, 2), ylim=c(0,5))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("B = ", B_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))


## MA(1): changing the block length (lambda) ----
cat("Now doing changing the block length in MA(1) model", "\n")
Sys.sleep(2)

# Settings
ma_order <- c(0,0,1)
B_all <- c(50,100,200,500,700,1000, 10000)
N_all <- c(50,100,200,500,700,1000)
N_lambda <- N_all[3]
ma_all <- c(-0.9,-0.7,-0.5,- 0.3,-0.1,0,0.1,0.3,0.5,0.7,0.9)
lambda_all <- c(round(N_lambda^(1/3)), round(N_lambda^(1/4)), round(N_lambda^(1/5)))
# b=T^(1/3), t^(1/4), t(1^5), or Politis and White (2004)
# ad-hoc or pwsd(y, correlogram = FALSE)

# Set an identifier
jj <- 0
# Set an empty list for each set of bootstrap
BB <- list()

#Loop over sample size/AR coefficients/no. of draws/lambdas 
for (j in lambda_all){
  jj <- jj + 1
  
  #Simulate y series 
  y <- arima.sim(list(order = ma_order, ma=ma_all[9]), n=N_all[3])
  
  #Sample size
  N <- NROW(y) 
  
  #Block size
  lambda1 <- lambda_all[jj]
  
  # Set the bootstrap draws
  B <- B_all[7]
  
  # Initialize matrix of resamples
  boot_est <- matrix(NA, B, 4)
  colnames(boot_est) <- c("mu_hat", "phi_hat", "mean", "variance")
  
  ### Bootstrapping
  for(i in 1:B)
  {
    #Resample via MBB
    yb_mbb <- MBB(y, b=lambda1)
    
    #Re-estimate and MA(1) for y
    out_boot <- arima(yb_mbb, order=c(1,0,0))
    coef_boot <- out_boot$coef
    mu_hat <- coef_boot[2]
    phi_hat <- coef_boot[1]
    
    # Saving estimates
    boot_est[i,1] <- mu_hat
    boot_est[i,2] <- phi_hat
    boot_est[i,3] <- mean(yb_mbb)
    boot_est[i,4] <- sd(yb_mbb)^2
    
    # For all loops, also add a loop tracker
    cat("Now doing ", i, " of ", B, " bootstrap draws, using block size l = ",
        lambda_all[jj], " (", jj, "/", NROW(lambda_all), ")", sep = "", "\n")
  }
  BB[[jj]] <- boot_est
  names(BB)[jj] <- paste("B = ", B_all[jj], sep = "")
}

### Plots
color <- c("#41ae76","#238b45","#005824")

# Check the density plots
pdf("plots/MBB_MA1_densities_diff_blocklength.pdf", width=11.27, height=8.69)

m <- matrix(c(1, 1, 2, 2, 3, 4, 5, 6), 
            nrow = 4, ncol = 2, byrow = TRUE)

layout(m, heights = c(0.06,0.06,0.44,0.44))

par(mar = c(0,0,0,0)) 
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Moving Block Bootstrap - Simulated MA(1) model",cex=1.5,font=2, pos = 1)
par(mar = c(0,0,0,0))
plot(1,1,type = "n",frame.plot = FALSE,axes = FALSE)
text(1,"Estimated distribution across different block size (lambda)",cex=1,font=2, pos = 1)

#mu
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,1]), 
     main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
     ylim=c(0,5), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(BB[[1]][,1]), sd=sd(BB[[1]][,1])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,1]), 
        main=expression(paste("Density Estimate using MC, ", hat(mu))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,1]), sd=sd(BB[[i]][,1])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#phi
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,2]), 
     main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
     xlim=c(0.1, 0.7), ylim=c(0,7))
curve(dnorm(x, mean=mean(BB[[1]][,2]), sd=sd(BB[[1]][,2])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,2]), 
        main=expression(paste("Density Estimate using MC, ", hat(phi))), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,2]), sd=sd(BB[[i]][,2])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n") 
#mean of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,3]), 
     main=expression("Density Estimate using MC, mean of y"), lwd=1,
     ylim=c(0,5), xlim=c(-0.5,0.5))
curve(dnorm(x, mean=mean(BB[[1]][,3]), sd=sd(BB[[1]][,3])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,3]), 
        main=expression("Density Estimate using MC, mean of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,3]), sd=sd(BB[[i]][,3])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")
#var of y
par(mar = c(5,1,1,1)+1) 
plot(density(BB[[1]][,4]), 
     main=expression("Density Estimate using MC, variance of y"), lwd=1,
     xlim=c(0.5, 2), ylim=c(0,3))
curve(dnorm(x, mean=mean(BB[[1]][,4]), sd=sd(BB[[1]][,4])), col=color[i],
      lty=2, lwd=1, add=TRUE, yaxt="n")
for (i in 2:length(BB))
{
  lines(density(BB[[i]][,4]), 
        main=expression("Density Estimate using MC, variance of y"), lwd=1,
        add=TRUE)
  curve(dnorm(x, mean=mean(BB[[i]][,4]), sd=sd(BB[[i]][,4])), col=color[i],
        lty=2, lwd=1, add=TRUE, yaxt="n")
}
legend(x = "topright",        
       legend = paste("lambda = ", lambda_all, sep=""),
       lty = 2,           
       col = color, 
       lwd = 2,
       cex = 0.5,
       bg = "transparent",
       bty = "n")

dev.off()

### Confidence Intervals
# bootstrap confidence intervals for the mean are
quantile(boot_est[,1], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,2], probs=c(2.5/100, 97.5/100))
quantile(boot_est[,3], probs=c(2.5/100, 97.5/100))

# Running time
toc(log = TRUE)
log.lst <- tic.log(format = FALSE)
timings <- unlist(lapply(log.lst, function(x) x$toc - x$tic))
total_time <- sum(timings)/60
cat("Simulations finished after", round(total_time, 1), "minutes", "\n")
tic.clearlog()
