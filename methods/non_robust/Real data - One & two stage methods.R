#####################################################################################
##### One-stage non-robust method for real dataset
#####################################################################################

### A: Loading the necessary packages and dataset ###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Desktop")
getwd()

rm(list = ls())

# Installing necessary libraries (if not already installed)
install.packages("nlme")
install.packages("ie2misc")
install.packages("caTools")
install.packages("foreach")
install.packages("doParallel")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)  # For splitting the dataset
library(caret)    # Another option for splitting
library(glmnet)
library(ie2misc)
library(MASS)
library(foreach)
library(doParallel)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
X_ww <- X
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
X_ss <- X
X_combined <- rbind(X_ww, X_ss)
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Set up parallel backend
cl <- makeCluster(detectCores() - 2)  # Leave one core free
registerDoParallel(cl)

####################################################################

num_iterations <- 100  # Number of iterations

# Initialize vectors to store results
cor_RMLA <- numeric(num_iterations)
mse_RMLA <- numeric(num_iterations)
mad_RMLA <- numeric(num_iterations)
cor_RMLV <- numeric(num_iterations)
mse_RMLV <- numeric(num_iterations)
mad_RMLV <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)
heritability_difference_values_RMLA <- numeric(num_iterations)
heritability_difference_values_RMLV <- numeric(num_iterations)
heritability_RMLA_values <- numeric(num_iterations)
heritability_RMLV_values <- numeric(num_iterations)

### B: Split the dataset into train and testing and fitting the model ###
########################################################################

results <- for(t in 1:num_iterations){
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    dataset= training_dataset
    colSums(dataset==2)
    dim(dataset)
    
    X_combined_train <- X_combined[training_indices, ]
    X_combined_test <- X_combined[-training_indices, ]
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  lmer(y ~ -1 +  as.factor(geno) + (1|rep)+(1|rep:sets), training_dataset)
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    r <- 2
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    X_combined_train[X_combined_train == 0] <- 0  # aa 
    X_combined_train[X_combined_train == 1] <- 1   # Aa 
    X_combined_train[X_combined_train == -1] <- 0   # AA 
    monomorphic_cols <- apply(X_combined_train, 2, function(col) length(unique(col)) == 1)
    X_combined__train_filtered <- X_combined_train[, !monomorphic_cols]
    h1 <- colSums(X_combined__train_filtered == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(X_combined__train_filtered == 1, na.rm = TRUE)  # Count of allele Aa
    
    # Use the number of columns in the filtered dataset
    num_snp_cols <- ncol(X_combined__train_filtered)
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(num_snp_cols)
    k <- numeric(num_snp_cols)
    l <- numeric(num_snp_cols)
    s <- numeric(num_snp_cols)
    su <- numeric(num_snp_cols)
    sum <- numeric(num_snp_cols+1)
    lamda_rmla <- numeric(num_snp_cols)
    lamda_rr <- numeric(num_snp_cols)
    dd_new <- data.frame(training_dataset$y, X_combined__train_filtered)
    
    # Use foreach loop to replace the inner for loops
    for(i in 1:num_snp_cols)  {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (370 - ((1 / 370) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:num_snp_cols)
    {
      lamda_rmla[j]=lam*(sum[num_snp_cols]/su[j])
      lamda_rr[j]=(1/0.16-1)*num_snp_cols*(sum[num_snp_cols+1]/su[j])
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(X_combined__train_filtered) %*% X_combined__train_filtered
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(X_combined__train_filtered) %*% training_dataset$y
    beta2 <- ginv(S3) %*% t(X_combined__train_filtered) %*% training_dataset$y
    
    g_p1 <- X_combined__train_filtered %*% beta1  # For RRWA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(training_dataset$y, g_p2)
    g_p2_vector <- as.vector(g_p2[, 1])
    M2 <- mse(training_dataset$y, g_p2_vector)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- training_dataset$y - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    
    ###  Step 9: RMLV  ###
    D <- (t(X_combined__train_filtered) %*% X_combined__train_filtered) + ((diag(lam, num_snp_cols)))
    cvv <- ginv(D)
    tr <- tr(cvv)
    sigma <- numeric(num_snp_cols)
    lamda_rmlv <- numeric(num_snp_cols)
    
    U <- beta1
    
    for(i in 1:num_snp_cols){
      sigma[i] <- ((U[i] * U[i]) - (Ve * tr)) / i
      lamda_rmlv[i] <- (Ve / sigma[i])
    }
    
    lamda_RMLV <- diag(lamda_rmlv)
    A <- t(X_combined__train_filtered) %*% training_dataset$y
    S1 <- ((t(X_combined__train_filtered) %*% X_combined__train_filtered) + lamda_RMLV)
    beta_RMLV <- solve(S1) %*% A
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    cr <- cor(training_dataset$y, g_pr)
    g_pr_vector <- as.vector(g_pr[, 1])
    M <- mse(training_dataset$y, g_pr_vector)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- training_dataset$y - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
    
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", t, e$message))
  })
}

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV)
average_MSE_RMLV <- mean(mse_RMLV)
average_MAD_RMLV <- mean(mad_RMLV)

average_cr_RMLA <- mean(cor_RMLA)
average_MSE_RMLA <- mean(mse_RMLA)
average_MAD_RMLA <- mean(mad_RMLA)

average_MSD_heritability <- mean(heritability_difference_values)
average_heritability <- mean(heritability_values)

average_heritability_RMLA <- mean(heritability_RMLA_values, na.rm = T)
average_heritability_RMLV <- mean(heritability_RMLV_values, na.rm = T)
heritability_MSD_RMLA <- mean(heritability_difference_values_RMLA, na.rm = T)
heritability_MSD_RMLV <- mean(heritability_difference_values_RMLV, na.rm = T)

cat("Average Heritability RMLA:", average_heritability_RMLA, "\n")
cat("Average Heritability RMLV:", average_heritability_RMLV, "\n")
cat("heritability_MSD_RMLA:", heritability_MSD_RMLA, "\n")
cat("heritability_MSD_RMLV:", heritability_MSD_RMLV, "\n")

# Print results
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", average_MSD_heritability, "\n")
cat("RMLA - Average Correlation:", average_cr_RMLA, "MSE:", average_MSE_RMLA, "MAD:", average_MAD_RMLA, "\n")
cat("RMLV - Average Correlation:", average_cr_RMLV, "MSE:", average_MSE_RMLV, "MAD:", average_MAD_RMLV, "\n")

################################## End of One-stage non-robust method for real dataset ################################################


#####################################################################################
##### Two-stage non-robust method for real dataset
#####################################################################################

### A: Loading the necessary packages and dataset ###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Desktop")
getwd()

rm(list = ls())

# Installing necessary libraries (if not already installed)
install.packages("nlme")
install.packages("ie2misc")
install.packages("caTools")
install.packages("foreach")
install.packages("doParallel")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)  # For splitting the dataset
library(caret)    # Another option for splitting
library(glmnet)
library(ie2misc)
library(foreach)
library(doParallel)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
X_ww <- X
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
X_ss <- X
X_combined <- rbind(X_ww, X_ss)
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Set up parallel backend
cl <- makeCluster(detectCores() - 2)  # Leave one core free
registerDoParallel(cl)

####################################################################

num_iterations <- 100  # Number of iterations

# Initialize vectors to store results
cor_RMLA <- numeric(num_iterations)
mse_RMLA <- numeric(num_iterations)
mad_RMLA <- numeric(num_iterations)
cor_RMLV <- numeric(num_iterations)
mse_RMLV <- numeric(num_iterations)
mad_RMLV <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)
heritability_RMLA_values <- numeric(num_iterations)
heritability_RMLV_values <- numeric(num_iterations)
heritability_difference_values_RMLA <- numeric(num_iterations)
heritability_difference_values_RMLV <- numeric(num_iterations)

### B: Split the dataset into train and testing and fitting the model ###
########################################################################

results <- for(t in 1:num_iterations){
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    dataset= training_dataset
    colSums(dataset==2)
    dim(dataset)
    
    X_combined_train <- X_combined[training_indices, ]
    X_combined_test <- X_combined[-training_indices, ]
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  rlmer(y ~ -1 +  as.factor(geno) + (1|rep)+(1|rep:sets), training_dataset,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    r <- 2
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    X_combined_train[X_combined_train == 0] <- 0  # aa 
    X_combined_train[X_combined_train == 1] <- 1   # Aa 
    X_combined_train[X_combined_train == -1] <- 0   # AA 
    monomorphic_cols <- apply(X_combined_train, 2, function(col) length(unique(col)) == 1)
    X_combined__train_filtered <- X_combined_train[, !monomorphic_cols]
    h1 <- colSums(X_combined__train_filtered == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(X_combined__train_filtered == 1, na.rm = TRUE)  # Count of allele Aa
    
    # Use the number of columns in the filtered dataset
    num_snp_cols <- ncol(X_combined__train_filtered)
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(num_snp_cols)
    k <- numeric(num_snp_cols)
    l <- numeric(num_snp_cols)
    s <- numeric(num_snp_cols)
    su <- numeric(num_snp_cols)
    sum <- numeric(num_snp_cols+1)
    lamda_rmla <- numeric(num_snp_cols)
    lamda_rr <- numeric(num_snp_cols)
    dd_new <- data.frame(training_dataset$y, X_combined__train_filtered)
    
    # Use foreach loop to replace the inner for loops
    for(i in 1:num_snp_cols)  {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (370 - ((1 / 370) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:num_snp_cols)
    {
      lamda_rmla[j]=lam*(sum[num_snp_cols]/su[j])
      lamda_rr[j]=(1/0.16-1)*num_snp_cols*(sum[num_snp_cols+1]/su[j])
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(X_combined__train_filtered) %*% X_combined__train_filtered
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(X_combined__train_filtered) %*% training_dataset$y
    beta2 <- ginv(S3) %*% t(X_combined__train_filtered) %*% training_dataset$y
    
    g_p1 <- X_combined__train_filtered %*% beta1  # For RRWA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(training_dataset$y, g_p2)
    g_p2_vector <- as.vector(g_p2[, 1])
    M2 <- mse(training_dataset$y, g_p2_vector)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- training_dataset$y - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    D <- (t(X_combined__train_filtered) %*% X_combined__train_filtered) + ((diag(lam, num_snp_cols)))
    cvv <- solve(D)
    tr <- tr(cvv)
    sigma <- numeric(num_snp_cols)
    lamda_rmlv <- numeric(num_snp_cols)
    
    U <- beta1
    
    for(i in 1:num_snp_cols){
      sigma[i] <- ((U[i] * U[i]) - (Ve * tr)) / i
      lamda_rmlv[i] <- (Ve / sigma[i])
    }
    
    lamda_RMLV <- diag(lamda_rmlv)
    A <- t(X_combined__train_filtered) %*% training_dataset$y
    S1 <- ((t(X_combined__train_filtered) %*% X_combined__train_filtered) + lamda_RMLV)
    beta_RMLV <- ginv(S1) %*% A
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    cr <- cor(training_dataset$y, g_pr)
    g_pr_vector <- as.vector(g_pr[, 1])
    M <- mse(training_dataset$y, g_pr_vector)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- training_dataset$y - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
    
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", t, e$message))
  })
}

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV,na.rm = T)
average_MSE_RMLV <- mean(mse_RMLV,na.rm = T)
average_MAD_RMLV <- mean(mad_RMLV,na.rm = T)

average_cr_RMLA <- mean(cor_RMLA,na.rm = T)
average_MSE_RMLA <- mean(mse_RMLA,na.rm = T)
average_MAD_RMLA <- mean(mad_RMLA,na.rm = T)

average_MSD_heritability <- mean(heritability_difference_values,na.rm = T)
average_heritability <- mean(heritability_values,na.rm = T)

average_heritability_RMLA <- mean(heritability_RMLA_values, na.rm = T)
average_heritability_RMLV <- mean(heritability_RMLV_values, na.rm = T)
heritability_MSD_RMLA <- mean(heritability_difference_values_RMLA, na.rm = T)
heritability_MSD_RMLV <- mean(heritability_difference_values_RMLV, na.rm = T)

cat("Average Heritability RMLA:", average_heritability_RMLA, "\n")
cat("Average Heritability RMLV:", average_heritability_RMLV, "\n")
cat("heritability_MSD_RMLA:", heritability_MSD_RMLA, "\n")
cat("heritability_MSD_RMLV:", heritability_MSD_RMLV, "\n")


# Print results
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", average_MSD_heritability, "\n")
cat("RMLA - Average Correlation:", average_cr_RMLA, "MSE:", average_MSE_RMLA, "MAD:", average_MAD_RMLA, "\n")
cat("RMLV - Average Correlation:", average_cr_RMLV, "MSE:", average_MSE_RMLV, "MAD:", average_MAD_RMLV, "\n")

############################### End of One-stage robust method for real dataset #######################################



#####################################################################################
##### Two-stage non-robust method for real dataset
#####################################################################################

### A: Loading the necessary packages and dataset ###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Desktop")
getwd()

rm(list = ls())

# Installing necessary libraries (if not already installed)
install.packages("nlme")
install.packages("ie2misc")
install.packages("caTools")
install.packages("foreach")
install.packages("doParallel")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)  # For splitting the dataset
library(caret)    # Another option for splitting
library(glmnet)
library(ie2misc)
library(MASS)
library(foreach)
library(doParallel)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
X_ww <- X
weights_ww <- as.matrix(weights)
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
X_ss <- X
weights_ss <- as.matrix(weights)
X_combined <- rbind(X_ww, X_ss)
weights_combined <- rbind(weights_ww,weights_ss)
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
registerDoParallel(cl)

####################################################################

num_iterations <- 10  # Number of iterations

# Initialize vectors to store results
cor_RMLA <- numeric(num_iterations)
mse_RMLA <- numeric(num_iterations)
mad_RMLA <- numeric(num_iterations)
cor_RMLV <- numeric(num_iterations)
mse_RMLV <- numeric(num_iterations)
mad_RMLV <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)
heritability_RMLA_values <- numeric(num_iterations)
heritability_RMLV_values <- numeric(num_iterations)
heritability_difference_values_RMLA <- numeric(num_iterations)
heritability_difference_values_RMLV <- numeric(num_iterations)

### B: Split the dataset into train and testing and fitting the model ###
########################################################################

results <- list()
mu_values <- list()

for(t in 1:num_iterations){
  
  # Print progress
  print(paste("Iteration:", t, "out of", num_iterations))
  
  # Assuming as[[t]] is your data, use caTools to split into train and test sets
  set.seed(123*t)  # For reproducibility
  
  #Data splitting and generation of train & testing data
  training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
  training_dataset  <- df[training_indices, ]
  testing_dataset <- df[-training_indices, ]
  
  dataset= training_dataset
  colSums(dataset==2)
  dim(dataset)
  
  ## Step 1: Fit a linear mixed model on the training data
  # (Use train_data instead of as[[t]])
  fit<-  lmer(y ~ -1 +  as.factor(geno) + (1|rep)+(1|rep:sets), training_dataset)
  
  mu = summary(fit)$coefficients[,1]
  mu_values <- fitted(fit)
  
  results[[t]] <- list(mu = mu_values, fit = fit)
  
}

# Results from the second loop
results_2 <- vector("list", num_iterations)  # Initialize results list for second loop

for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    set.seed(123*t)  # For reproducibility
    
    t=1
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    X_combined_train <- X_combined[training_indices, ]
    X_combined_test <- X_combined[-training_indices, ]
    
    weights_combined_train <- weights_combined[training_indices, ]
    
    # Retrieve fitted values for the current iteration
    mu_values <- results[[t]]$mu
    
    fit_2 = mixed.solve(mu_values,X_combined_train, method="REML") ##2nd stage with REML estimate
    beta=fit_2$u
    g_p=X_combined_train%*%beta
    
    c=cor(training_dataset$y,g_p)
    g_p_vector <- as.vector(g_p[, 1])
    M=mse(training_dataset$y,g_p_vector) ##RM11
    
    ### Step 2: Calculate variance components on the training data ###
    fit<-  lmer(y ~ -1 +  as.factor(geno) + (1|rep)+(1|rep:sets), training_dataset)
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    r=2
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    X_combined_train[X_combined_train == 0] <- 0  # aa 
    X_combined_train[X_combined_train == 1] <- 1   # Aa 
    X_combined_train[X_combined_train == -1] <- 0   # AA 
    monomorphic_cols <- apply(X_combined_train, 2, function(col) length(unique(col)) == 1)
    X_combined__train_filtered <- X_combined_train[, !monomorphic_cols]
    h1 <- colSums(X_combined__train_filtered == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(X_combined__train_filtered == 1, na.rm = TRUE)  # Count of allele Aa
    
    # Use the number of columns in the filtered dataset
    num_snp_cols <- ncol(X_combined__train_filtered)
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(num_snp_cols)
    k <- numeric(num_snp_cols)
    l <- numeric(num_snp_cols)
    s <- numeric(num_snp_cols)
    su <- numeric(num_snp_cols)
    sum <- numeric(num_snp_cols+1)
    lamda_rmla <- numeric(num_snp_cols)
    lamda_rr <- numeric(num_snp_cols)
    dd_new <- data.frame(training_dataset$y, X_combined__train_filtered)
    
    # Use foreach loop to replace the inner for loops
    for(i in 1:num_snp_cols)  {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (370 - ((1 / 370) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:num_snp_cols)
    {
      lamda_rmla[j]=lam*(sum[num_snp_cols]/su[j])
      lamda_rr[j]=(1/0.16-1)*num_snp_cols*(sum[num_snp_cols+1]/su[j])
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(X_combined__train_filtered) %*% X_combined__train_filtered
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(X_combined__train_filtered) %*% training_dataset$y
    beta2 <- ginv(S3) %*% t(X_combined__train_filtered) %*% training_dataset$y
    
    g_p1 <- X_combined__train_filtered %*% beta1  # For RRWA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(training_dataset$y, g_p2)
    g_p2_vector <- as.vector(g_p2[, 1])
    M2 <- mse(training_dataset$y, g_p2_vector)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- training_dataset$y - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    u=as.matrix(fit_2$u)
    lam=(fit_2$Ve/fit_2$Vu)
    D=(t(X_combined__train_filtered)%*% X_combined__train_filtered) + ((diag(rep(lam,num_snp_cols))))
    cvv=solve(D)
    tr=tr(cvv)
    sigma=0
    lamda_rmlv=0
    for(i in 1:num_snp_cols)
    {
      sigma[i]=((u[i,]*u[i,])-(fit_2$Ve*t))/i
      lamda_rmlv[i]=(fit_2$Ve/sigma[i])
    }
    lamda_RMLV=diag(lamda_rmlv)
    A=t(X_combined__train_filtered)%*% training_dataset$y
    S1=((t(X_combined__train_filtered)%*%X_combined__train_filtered) + lamda_RMLV)
    beta_RMLV=ginv(S1)%*%A
    g_pr=X_combined__train_filtered%*%beta_RMLV
    
    cr <- cor(training_dataset$y, g_pr)
    g_pr_vector <- as.vector(g_pr[, 1])
    M <- mse(training_dataset$y, g_pr_vector)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- training_dataset$y - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Combine results for this iteration
    c(cr, M, c2, M2, mad_RMLA[t], mad_RMLV[t],  heritability_values, heritability_difference_values)
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", t, e$message))
  })
}


# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV,na.rm=T)
average_MSE_RMLV <- mean(mse_RMLV,na.rm=T)
average_MAD_RMLV <- mean(mad_RMLV,na.rm=T)

average_cr_RMLA <- mean(cor_RMLA,na.rm=T)
average_MSE_RMLA <- mean(mse_RMLA,na.rm=T)
average_MAD_RMLA <- mean(mad_RMLA,na.rm=T)

average_MSD_heritability <- mean(heritability_difference_values,na.rm=T)
average_heritability <- mean(heritability_values,na.rm=T)

average_heritability_RMLA <- mean(heritability_RMLA_values, na.rm = T)
average_heritability_RMLV <- mean(heritability_RMLV_values, na.rm = T)
heritability_MSD_RMLA <- mean(heritability_difference_values_RMLA, na.rm = T)
heritability_MSD_RMLV <- mean(heritability_difference_values_RMLV, na.rm = T)

cat("Average Heritability RMLA:", average_heritability_RMLA, "\n")
cat("Average Heritability RMLV:", average_heritability_RMLV, "\n")
cat("heritability_MSD_RMLA:", heritability_MSD_RMLA, "\n")
cat("heritability_MSD_RMLV:", heritability_MSD_RMLV, "\n")

# Print results
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", average_MSD_heritability, "\n")
cat("RMLA - Average Correlation:", average_cr_RMLA, "MSE:", average_MSE_RMLA, "MAD:", average_MAD_RMLA, "\n")
cat("RMLV - Average Correlation:", average_cr_RMLV, "MSE:", average_MSE_RMLV, "MAD:", average_MAD_RMLV, "\n")

######################################## End of Two-stage non-robust method for real dataset ###########################################



#####################################################################################
##### Two-stage robust method for real dataset
#####################################################################################

### A: Loading the necessary packages and dataset ###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Desktop")
getwd()

rm(list = ls())

# Installing necessary libraries (if not already installed)
install.packages("nlme")
install.packages("ie2misc")
install.packages("caTools")
install.packages("foreach")
install.packages("doParallel")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)  # For splitting the dataset
library(caret)    # Another option for splitting
library(glmnet)
library(ie2misc)
library(MASS)
library(foreach)
library(doParallel)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
X_ww <- X
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
X_ss <- X
X_combined <- rbind(X_ww, X_ss)
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
registerDoParallel(cl)

####################################################################

num_iterations <- 10  # Number of iterations

# Initialize vectors to store results
cor_RMLA <- numeric(num_iterations)
mse_RMLA <- numeric(num_iterations)
mad_RMLA <- numeric(num_iterations)
cor_RMLV <- numeric(num_iterations)
mse_RMLV <- numeric(num_iterations)
mad_RMLV <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)
heritability_RMLA_values <- numeric(num_iterations)
heritability_RMLV_values <- numeric(num_iterations)
heritability_difference_values_RMLA <- numeric(num_iterations)
heritability_difference_values_RMLV <- numeric(num_iterations)

### B: Split the dataset into train and testing and fitting the model ###
########################################################################

results <- list()

for(t in 1:num_iterations){
  
  # Print progress
  print(paste("Iteration:", t, "out of", num_iterations))
  
  # Assuming as[[t]] is your data, use caTools to split into train and test sets
  set.seed(123)  # For reproducibility
  
  #Data splitting and generation of train & testing data
  training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
  training_dataset  <- df[training_indices, ]
  testing_dataset <- df[-training_indices, ]
  
  dataset= training_dataset
  colSums(dataset==2)
  dim(dataset)
  
  ## Step 1: Fit a linear mixed model on the training data
  # (Use train_data instead of as[[t]])
  fit<-  rlmer(y ~ -1 + as.factor(geno) + (1|rep)+(1|rep:sets),dataset ,
               rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
  
  mu = summary(fit)$coefficients[,1]
  mu_values <- fitted(fit)
  
  results[[t]] <- list(mu = mu_values, fit = fit)
  
}

# Results from the second loop
results_2 <- vector("list", num_iterations)  # Initialize results list for second loop

for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    set.seed(123*t)  # For reproducibility
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    X_combined_train <- X_combined[training_indices, ]
    X_combined_test <- X_combined[-training_indices, ]
    
    # Retrieve fitted values for the current iteration
    mu_values <- results[[t]]$mu
    
    fit_2 = mixed.solve(mu_values,X_combined_train, method="REML") ##2nd stage with REML estimate
    beta=fit_2$u
    g_p=X_combined_train%*%beta
    
    c=cor(training_dataset$y,g_p)
    g_p_vector <- as.vector(g_p[, 1])
    M=mse(training_dataset$y,g_p_vector) ##RM11
    
    ### Step 2: Calculate variance components on the training data ###
    fit<-  rlmer(y ~ 1 + (1|rep)+(1|rep:sets),dataset ,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    r=2
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    X_combined_train[X_combined_train == 0] <- 0  # aa 
    X_combined_train[X_combined_train == 1] <- 1   # Aa 
    X_combined_train[X_combined_train == -1] <- 0   # AA 
    monomorphic_cols <- apply(X_combined_train, 2, function(col) length(unique(col)) == 1)
    X_combined__train_filtered <- X_combined_train[, !monomorphic_cols]
    h1 <- colSums(X_combined__train_filtered == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(X_combined__train_filtered == 1, na.rm = TRUE)  # Count of allele Aa
    
    # Use the number of columns in the filtered dataset
    num_snp_cols <- ncol(X_combined__train_filtered)
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(num_snp_cols)
    k <- numeric(num_snp_cols)
    l <- numeric(num_snp_cols)
    s <- numeric(num_snp_cols)
    su <- numeric(num_snp_cols)
    sum <- numeric(num_snp_cols+1)
    lamda_rmla <- numeric(num_snp_cols)
    lamda_rr <- numeric(num_snp_cols)
    dd_new <- data.frame(training_dataset$y, X_combined__train_filtered)
    
    # Use foreach loop to replace the inner for loops
    for(i in 1:num_snp_cols)  {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (370 - ((1 / 370) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:num_snp_cols)
    {
      lamda_rmla[j]=lam*(sum[num_snp_cols]/su[j])
      lamda_rr[j]=(1/0.16-1)*num_snp_cols*(sum[num_snp_cols+1]/su[j])
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(X_combined__train_filtered) %*% X_combined__train_filtered
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(X_combined__train_filtered) %*% training_dataset$y
    beta2 <- ginv(S3) %*% t(X_combined__train_filtered) %*% training_dataset$y
    
    g_p1 <- X_combined__train_filtered %*% beta1  # For RRWA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(training_dataset$y, g_p2)
    g_p2_vector <- as.vector(g_p2[, 1])
    M2 <- mse(training_dataset$y, g_p2_vector)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- X_combined__train_filtered %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- training_dataset$y - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    
    ###  Step 9: RMLV  ###
    u=as.matrix(fit_2$u)
    lam=(fit_2$Ve/fit_2$Vu)
    D=(t(X_combined__train_filtered)%*% X_combined__train_filtered) + ((diag(rep(lam,num_snp_cols))))
    cvv=solve(D)
    tr=tr(cvv)
    sigma=0
    lamda_rmlv=0
    for(i in 1:num_snp_cols)
    {
      sigma[i]=((u[i,]*u[i,])-(fit_2$Ve*t))/i
      lamda_rmlv[i]=(fit_2$Ve/sigma[i])
    }
    lamda_RMLV=diag(lamda_rmlv)
    A=t(X_combined__train_filtered)%*% training_dataset$y
    S1=((t(X_combined__train_filtered)%*%X_combined__train_filtered) + lamda_RMLV)
    beta_RMLV=ginv(S1)%*%A
    g_pr=X_combined__train_filtered%*%beta_RMLV
    
    cr <- cor(training_dataset$y, g_pr)
    g_pr_vector <- as.vector(g_pr[, 1])
    M <- mse(training_dataset$y, g_pr_vector)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- X_combined__train_filtered %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- training_dataset$y - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Combine results for this iteration
    c(cr, M, c2, M2, mad_RMLA[t], mad_RMLV[t],  heritability_values, heritability_difference_values)
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", t, e$message))
  })
}



# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV,na.rm=T)
average_MSE_RMLV <- mean(mse_RMLV,na.rm=T)
average_MAD_RMLV <- mean(mad_RMLV,na.rm=T)

average_cr_RMLA <- mean(cor_RMLA,na.rm=T)
average_MSE_RMLA <- mean(mse_RMLA,na.rm=T)
average_MAD_RMLA <- mean(mad_RMLA,na.rm=T)

average_MSD_heritability <- mean(heritability_difference_values,na.rm=T)
average_heritability <- mean(heritability_values,na.rm=T)

average_heritability_RMLA <- mean(heritability_RMLA_values, na.rm = T)
average_heritability_RMLV <- mean(heritability_RMLV_values, na.rm = T)
heritability_MSD_RMLA <- mean(heritability_difference_values_RMLA, na.rm = T)
heritability_MSD_RMLV <- mean(heritability_difference_values_RMLV, na.rm = T)

cat("Average Heritability RMLA:", average_heritability_RMLA, "\n")
cat("Average Heritability RMLV:", average_heritability_RMLV, "\n")
cat("heritability_MSD_RMLA:", heritability_MSD_RMLA, "\n")
cat("heritability_MSD_RMLV:", heritability_MSD_RMLV, "\n")

# Print results
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", average_MSD_heritability, "\n")
cat("RMLA - Average Correlation:", average_cr_RMLA, "MSE:", average_MSE_RMLA, "MAD:", average_MAD_RMLA, "\n")
cat("RMLV - Average Correlation:", average_cr_RMLV, "MSE:", average_MSE_RMLV, "MAD:", average_MAD_RMLV, "\n")

######################################## End of Two-stage non-robust method for real dataset ###########################################


