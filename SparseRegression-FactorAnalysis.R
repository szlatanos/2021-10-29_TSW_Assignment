#Set up ----
# Clear environment
rm(list=ls())

# Load necessary libraries
library("tidyverse")
library("lubridate")
library("boeCharts")
library("zoo")
library("BVAR")
library("tseries")
library("fredr")
library("glmnet")
library("pls")
library("changepoint")
library("fbi")

# Load custom functions
source("functions/insample-functions.R")

# Data ----
# Long format
df_q <- fred_qd %>% 
  mutate(date = as.Date(rownames(.)),
         date_qrt = as.yearqtr(date)) %>%
  select(date, date_qrt, everything()) %>%
  pivot_longer(cols = 3:ncol(.))

#Get series names
series_id <- df_q %>% distinct(name) %>% pull(name)
fredqd_df <- fredqd_description %>% select(fred_mnemonic, description, group)

## Plots ----
# Raw data

p_q <- ggplot(df_q, aes(x=date, y=value, colour = name)) +
  geom_line(show.legend = FALSE) +
  labs(
    title = "Raw quarterly indicators", 
    subtitle = "FRED-QD",
    y = "Value", x = NULL
  ) +
  theme_overground() +
  # apply custom axis settings
  scale_y_continuous(
    expand = c(0, 0), breaks = boe_breaks_numeric(), 
    limits = boe_limits_numeric(), position = "left"
  )
p_q

# Filtering ----
## Sample size
smpl_start <- as.Date("1980-01-01")
smpl_end <- as.Date("2019-12-01")


df_q <- df_q %>%
  filter(date >= smpl_start) %>%
  filter(date <= smpl_end)

## NAs ----
# Find series which don't start with the smpl_start date but later
df_q_na <- df_q %>% 
  filter(is.na(value)) %>%
  group_by(name) %>%
  slice(which(is.na(value))) %>%
  filter(date <= as.Date(today()-years(2))) %>%
  slice(which.max(date)) %>%
  ungroup() %>%
  select(date, date_qrt, name) %>%
  arrange(date) %>%
  pull(name)

#Filter out these short series
# (i.e. series with starting date after smpl_start)

df_q <- df_q %>%
  filter(!name %in% df_q_na)

## Clean dataset ----
df <- df_q %>%
  select(-date) %>%
  arrange(name) %>%
  pivot_wider(1:3, names_from = name, values_from = value)

data <- as.matrix(df[,2:ncol(df)])
dates <- as.matrix(df[,1])

# Stationarity ----
# Remove the NAs, if any
data <- na.omit(data)

## ADF test for stationarity on raw series ----
# Check the all series are stationary
ADF <- NULL
for(j in 1:NCOL(data))
{
  adf.out <- adf.test(data[,j])
  adf.pval <- adf.out$p.value
  ADF <- c(ADF, adf.pval)
}

non_stat <- which(ADF>0.1) #clearly most, if not all of the series, are non-stationary
non_stat_name <- colnames(data[,non_stat])

df_f <- df %>% select(-date_qrt)

## Transformation ----
# Save the transformation codes
codes_t <- fred_code(vars = colnames(data), table = T)

#Transform the series
df_t <- fred_transform(df_f, type = "fred_qd")

#Get data in a matrix form
data_t <- as.matrix(df_t)
dates_t <- dates[3:NROW(dates)]

## ADF test for stationarity on transformed series ----
# ADF test for stationarity
# Check the all series are stationary
ADF <- NULL

for(j in 1:NCOL(data_t))
{
  adf.out <- adf.test(data_t[,j])
  adf.pval <- adf.out$p.value
  ADF <- c(ADF, adf.pval)
}
non_stat <- which(ADF>0.1)
non_stat_name <- colnames(data[,non_stat])

cat("The following variables are non-stationary, even after the suggested transformation:", 
    str_c(non_stat_name, collapse = ", "), "\n")

non_stat_comp <- left_join(codes_t, 
                           data.frame(variable = non_stat_name, non_stat = "1"), 
                           by = "variable") 

# Plot non-stationary series after transformation
non_stat_plots <- list()

for(j in 1:NROW(non_stat))
{
  p <- ggplot() +
    geom_line(aes_string(x=seq_along(data_t[,j]), y=data_t[,j]),
              show.legend = FALSE) +
    labs(
      title = colnames(data_t)[j],
      subtitle = "",
      y = "Value", x = NULL
    ) +
    theme_overground() +
    # apply custom axis settings
    scale_y_continuous(
      expand = c(0, 0), breaks = boe_breaks_numeric(),
      limits = boe_limits_numeric(), position = "left"
    )
  non_stat_plots[[j]] <- p
 
  names(non_stat_plots)[j] <- colnames(data_t)[j]
}

# Drop non-stationary series after transformation
df_t <- df_t[,-c(non_stat)]

transf_series_plots <- list()

for(j in 1:NCOL(df_t))
{
  p <- ggplot() +
    geom_line(aes_string(x=as.yearqtr(dates_t), y=df_t[,j]),
              show.legend = FALSE) +
    labs(
      title = colnames(df_t)[j],
      subtitle = "",
      y = "Value", x = NULL
    ) +
    theme_overground() +
    scale_y_continuous(
      expand = c(0, 0), breaks = boe_breaks_numeric(),
      limits = boe_limits_numeric(), position = "left"
    )
  transf_series_plots[[j]] <- p
  
  names(transf_series_plots)[j] <- colnames(df_t)[j]
}

# Clean data ----
## Final (transformed) dataset
# df_t -> transformed data
# dates_t -> corresponding dates

# Xs
X <- as.matrix(df_t[!names(df_t) %in% c("GDPC1")]) #X variables
rownames(X) <- dates_t
# Y Target variable
Y <- as.matrix(df_t["GDPC1"])
rownames(Y) <- dates_t

d <- as.Date(as.yearqtr(dates_t, 
  format = "%Y Q%q"))

## Plots ----
#Transformed series with change points (BIC criterion)
pdf("plots/transformed_series.pdf", width= 8.69, height=11.27)
xx <- df_t
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  #Title
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 40)
  
  #Plot
  plot(d, xx[,i], 
       ylab = "value", xlab = "time",
       col="blue", lwd=2, type="l")
  
  title(main = title, cex.main = 0.8)
}
dev.off()


#Transformed series with change points (BIC criterion)
pdf("plots/transformed_series_breaks.pdf", width= 8.69, height=11.27)
xx <- df_t
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  #Title
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 40)

  #Change points in mean
  xx_br <- cpt.mean(data = xx[,i],
           penalty = "BIC", 
           method = "BinSeg", Q = 20, 
           class = T)
  #Plot
  plot(xx_br, 
       ylab = "value", xlab = "index",
       col="blue", lwd=2, type="l")

  title(main = title, cex.main = 0.8)
}
dev.off()


# Ridge ----
#Discuss how the beta coefficients change 
#for each variable during the last 20 years

## Background info:
# Y, X 
# 158 obs -> quarterly
# how my betas change in the past 20 years * 4 quarters = 80 periods

## Create the matrix to store your results
Nout <- 20*4   # period to check your betas across time
ridge_betas <- matrix(NA, NROW(Y), NCOL(X)+1) # +1 for your intercept
rownames(ridge_betas) <- rownames(Y)
colnames(ridge_betas) <- c("Intercept", colnames(X))
j <- 1

# The loop below gives us the betas across time in a recursive estimation
for(i in (NROW(Y)-Nout+1):NROW(Y))
{
  i_in <- 1:i
  X_in <- X[i_in,]
  Y_in <- Y[i_in,]
  
  # Estimation
  out <- cv.glmnet(X_in, Y_in, type.measure="mse", alpha=0, #set = 0 for Ridge
                   family="gaussian", standardize=TRUE,
                   nfolds=NROW(Y_in))
  ll <- out$lambda.min
  bb <- as.matrix(coef(out, ll))
  ridge_betas[i,] <- bb
  
  cat("Now estimating with sample size ", i, " (", j , " of ", Nout, 
      " different sample sizes)", "\n", sep = "")
  j <- j + 1
}

# Remove NAs
ridge_betas <- na.omit(ridge_betas)
d <- as.Date(as.yearqtr(
  dates_t[(NROW(Y)-Nout+1):NROW(Y)], 
  format = "%Y Q%q"))

#Long format
ridge_betas_long <- ridge_betas %>% 
  as.data.frame(.) %>% 
  mutate(date = as.Date(as.yearqtr(rownames(.), format = "%Y Q%q"))) %>% 
           select(date, everything()) %>% 
           pivot_longer(., cols = 2:NCOL(.), names_to = "variable", values_to = "value")
  
## Plots
pdf("plots/ridge_betas.pdf", width= 8.69, height=11.27)
xx <- ridge_betas
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 60)
  
  plot(d, xx[,i], 
       ylab = "coefficient", xlab = "time",
       col="blue", lwd=2, type="l")
  title(main = title, cex.main = 0.8)
  
}
dev.off()

# Lasso ----
#Discuss how the beta coefficients change 
#for each variable during the last 20 years

## Background info:
# Y, X 
# 158 obs -> quarterly
# how my betas change in the past 20 years * 4 quarters = 80 periods

## Create the matrix to store your results
Nout <- 20*4   # period to check your betas across time
lasso_betas <- matrix(NA, NROW(Y), NCOL(X)+1) # +1 for your intercept
rownames(lasso_betas) <- rownames(Y)
colnames(lasso_betas) <- c("Intercept", colnames(X))
j <- 1

# The loop below gives us the betas across time in a recursive estimation
for(i in (NROW(Y)-Nout+1):NROW(Y))
{
  i_in <- 1:i
  X_in <- X[i_in,]
  Y_in <- Y[i_in,]
  
  # Estimation
  out <- cv.glmnet(X_in, Y_in, type.measure="mse", alpha=1, # set = 1 for Lasso
                   family="gaussian", standardize=TRUE,
                   nfolds=NROW(Y_in))
  ll <- out$lambda.min
  bb <- as.matrix(coef(out, ll))
  lasso_betas[i,] <- bb
  
  cat("Now estimating with sample size ", i, " (", j , " of ", Nout, 
      " different sample sizes)", "\n", sep = "")
  j <- j + 1
}

# Remove NAs
lasso_betas <- na.omit(lasso_betas)
d <- as.Date(as.yearqtr(
  dates_t[(NROW(Y)-Nout+1):NROW(Y)], 
  format = "%Y Q%q"))

#Long format
lasso_betas_long <- lasso_betas %>% 
  as.data.frame(.) %>% 
  mutate(date = as.Date(as.yearqtr(rownames(.), format = "%Y Q%q"))) %>% 
  select(date, everything()) %>% 
  pivot_longer(., cols = 2:NCOL(.), names_to = "variable", values_to = "value")

## Plots
pdf("plots/lasso_betas.pdf", width= 8.69, height=11.27)
xx <- lasso_betas
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 60)
  
  plot(d, xx[,i], 
       ylab = "coefficient", xlab = "time",
       col="blue", lwd=2, type="l")
  title(main = title, cex.main = 0.8)
  
}
dev.off()

# PCA ----
#Discuss how the factor loadings change 
#for each variable during the last 20 years

## Background info:
# Y, X 
# 158 obs -> quarterly
# how my betas change in the past 20 years * 4 quarters = 80 periods

## Create the matrix to store your results
Nout <- 20*4   # period to check your betas across time
pca_loads <- matrix(NA, NROW(Y), NCOL(X))
rownames(pca_loads) <- rownames(Y)
colnames(pca_loads) <- colnames(X)
j <- 1

# The loop below gives us the betas across time in a recursive estimation
for(i in (NROW(Y)-Nout+1):NROW(Y))
{
  i_in <- 1:i
  X_in <- X[i_in,]
  Y_in <- Y[i_in,]
  
  # first PCA factor
  out <- pcr(Y_in~X_in, ncomp=1, scale=TRUE, center=TRUE)
  floads <- as.numeric(out$loadings)
  pca_loads[i,] <- floads
  
  cat("Now estimating with sample size ", i, " (", j , " of ", Nout, 
      " different sample sizes)", "\n", sep = "")
  j <- j + 1
}

# Remove NAs
pca_loads <- na.omit(pca_loads)
d <- as.Date(as.yearqtr(
  dates_t[(NROW(Y)-Nout+1):NROW(Y)], 
  format = "%Y Q%q"))

#Long format
pca_loads_long <- pca_loads %>% 
  as.data.frame(.) %>% 
  mutate(date = as.Date(as.yearqtr(rownames(.), format = "%Y Q%q"))) %>% 
  select(date, everything()) %>% 
  pivot_longer(., cols = 2:NCOL(.), names_to = "variable", values_to = "value")

## Plots
pdf("plots/pca_loads.pdf", width= 8.69, height=11.27)
xx <- pca_loads
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 60)
  
  plot(d, xx[,i], 
       ylab = "coefficient", xlab = "time",
       col="blue", lwd=2, type="l")
  title(main = title, cex.main = 0.8)
  
}
dev.off()

# PLS ----
#Discuss how the factor loadings change 
#for each variable during the last 20 years

## Background info:
# Y, X 
# 158 obs -> quarterly
# how my betas change in the past 20 years * 4 quarters = 80 periods

## Create the matrix to store your results
Nout <- 20*4   # period to check your betas across time
pls_loads <- matrix(NA, NROW(Y), NCOL(X))
rownames(pls_loads) <- rownames(Y)
colnames(pls_loads) <- colnames(X)
j <- 1

# The loop below gives us the betas across time in a recursive estimation
for(i in (NROW(Y)-Nout+1):NROW(Y))
{
  i_in <- 1:i
  X_in <- X[i_in,]
  Y_in <- Y[i_in,]
  
  # first pls factor
  out <- plsr(Y_in~X_in, ncomp=1, scale=TRUE, center=TRUE)
  floads <- as.numeric(out$loadings)
  pls_loads[i,] <- floads
  
  cat("Now estimating with sample size ", i, " (", j , " of ", Nout, 
      " different sample sizes)", "\n", sep = "")
  j <- j + 1
}

# Remove NAs
pls_loads <- na.omit(pls_loads)
d <- as.Date(as.yearqtr(
  dates_t[(NROW(Y)-Nout+1):NROW(Y)], 
  format = "%Y Q%q"))

#Long format
pls_loads_long <- pls_loads %>% 
  as.data.frame(.) %>% 
  mutate(date = as.Date(as.yearqtr(rownames(.), format = "%Y Q%q"))) %>% 
  select(date, everything()) %>% 
  pivot_longer(., cols = 2:NCOL(.), names_to = "variable", values_to = "value")

## Plots
pdf("plots/pls_loads.pdf", width= 8.69, height=11.27)
xx <- pls_loads
par(mfrow=c(4,3))
for(i in 1:NCOL(xx))
{
  title <- fredqd_df %>% 
    filter(fred_mnemonic == colnames(xx)[i]) %>% 
    pull(description) %>% str_wrap(., width = 60)
  
  plot(d, xx[,i], 
       ylab = "coefficient", xlab = "time",
       col="blue", lwd=2, type="l")
  title(main = title, cex.main = 0.8)
}
dev.off()