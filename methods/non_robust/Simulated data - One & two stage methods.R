######################################################################################
######## One stage non-robust method
######################################################################################


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
library(lme4)
library(insight)
library(caTools)  # For splitting the dataset
library(caret)    # Another option for splitting
library(glmnet)
library(ie2misc)
library(MASS)
library(foreach)
library(doParallel)

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
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

results <- foreach(t = 1:num_iterations, .combine = 'rbind', .packages = c('caTools', 'lme4', 'Metrics', 'robustlmm', 'insight','MASS','ie2misc','lme4')) %dopar% {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Create a split: 70% training, 30% testing
    split <- sample.split(as[[t]]$y, SplitRatio = 0.7)
    
    # Subset the data into train and test
    train_data <- subset(as[[t]], split == TRUE)
    test_data <- subset(as[[t]], split == FALSE)
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  lmer(y ~ 1 + (1 | `rep(bvalues, 2)`) + (1|rep) + (1|rep:block_effects), train_data)
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    h1 <- colSums(Z_replicated == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(Z_replicated == 1, na.rm = TRUE)  # Count of allele Aa
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(11640)
    k <- numeric(11640)
    l <- numeric(11640)
    s <- numeric(11640)
    su <- numeric(11640)
    sum <- numeric(11641)
    lamda_rmla <- numeric(11640)
    lamda_rr <- numeric(11640)
    dd_new <- data.frame(res[,t], a[[t]])
    
    for(i in 1:11640) {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (1430 - ((1 / 1430) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:11640) {
      lamda_rmla[j] <- lam * (sum[11641] / su[j])
      lamda_rr[j] <- (1 / 0.16 - 1) * 11640 * (sum[11641] / su[j])
      lamda_rmla[j]
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(a[[t]]) %*% a[[t]]
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(a[[t]]) %*% as[[t]]$y
    beta2 <- ginv(S3) %*% t(a[[t]]) %*% as[[t]]$y
    str(beta1)
    str(beta2)
    
    g_p1 <- a[[t]] %*% beta1  # For RRWA
    g_p2 <- a[[t]] %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(res[,t], g_p2)
    M2 <- mse(res[,t], g_p2)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- a[[t]] %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- res[,t] - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    D <- (t(a[[t]]) %*% a[[t]]) + ((diag(rep(lam, 11640))))
    cvv <- solve(D)
    tr <- tr(cvv)
    sigma <- numeric(11640)
    lamda_rmlv <- numeric(11640)
    
    U <- beta1
    
    for(j in 1:11640) {
      sigma[i] <- ((U[i] * U[i]) - (Ve * tr)) / i
      lamda_rmlv[i] <- (Ve / sigma[i])
    }
    
    lamda_RMLV <- diag(lamda_rmlv)
    A <- t(a[[t]]) %*% res[, t]
    S1 <- ((t(a[[t]]) %*% a[[t]]) + lamda_RMLV)
    beta_RMLV <- ginv(S1) %*% A
    g_pr <- a[[t]] %*% beta_RMLV
    
    cr <- cor(res[,t], g_pr)
    M <- mse(res[,t], g_pr)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- a[[t]] %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- res[,t] - g_pr
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
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

# Stop the parallel cluster
stopCluster(cl)

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV)
average_MSE_RMLV <- mean(mse_RMLV)
average_MAD_RMLV <- mean(mad_RMLV)

average_cr_RMLA <- mean(cor_RMLA)
average_MSE_RMLA <- mean(mse_RMLA)
average_MAD_RMLA <- mean(mad_RMLA)

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

############################# END OF ONE-STAGE NON-ROBUST METHOD ##########################################


##############################################################################
##### One-stage robust method 
##############################################################################

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

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
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

results <- foreach(t = 1:num_iterations, .combine = 'rbind', .packages = c('caTools', 'lme4', 'Metrics', 'robustlmm', 'insight','MASS','ie2misc','lme4')) %dopar% {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Create a split: 70% training, 30% testing
    split <- sample.split(as[[t]]$y, SplitRatio = 0.7)
    
    # Subset the data into train and test
    train_data <- subset(as[[t]], split == TRUE)
    test_data <- subset(as[[t]], split == FALSE)
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  rlmer(y ~ 1 + (1 | `rep(genotype_values)`) + (1|rep) + (1|rep:block_effects), train_data, rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    h1 <- colSums(Z_replicated == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(Z_replicated == 1, na.rm = TRUE)  # Count of allele Aa
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(11640)
    k <- numeric(11640)
    l <- numeric(11640)
    s <- numeric(11640)
    su <- numeric(11640)
    sum <- numeric(11641)
    lamda_rmla <- numeric(11640)
    lamda_rr <- numeric(11640)
    dd_new <- data.frame(res[,t], a[[t]])
    
    for(i in 1:11640) {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (1430 - ((1 / 1430) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:11640) {
      lamda_rmla[j] <- lam * (sum[11641] / su[j])
      lamda_rr[j] <- (1 / 0.16 - 1) * 11640 * (sum[11641] / su[j])
      lamda_rmla[j]
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(a[[t]]) %*% a[[t]]
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- ginv(S2) %*% t(a[[t]]) %*% as[[t]]$y
    beta2 <- ginv(S3) %*% t(a[[t]]) %*% as[[t]]$y
    str(beta1)
    str(beta2)
    
    g_p1 <- a[[t]] %*% beta1  # For RRWA
    g_p2 <- a[[t]] %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(res[,t], g_p2)
    M2 <- mse(res[,t], g_p2)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- a[[t]] %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- res[,t] - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    D <- (t(a[[t]]) %*% a[[t]]) + ((diag(rep(lam, 11640))))
    cvv <- solve(D)
    tr <- tr(cvv)
    sigma <- numeric(11640)
    lamda_rmlv <- numeric(11640)
    
    U <- beta1
    
    for(j in 1:11640) {
      sigma[i] <- ((U[i] * U[i]) - (Ve * tr)) / i
      lamda_rmlv[i] <- (Ve / sigma[i])
    }
    
    lamda_RMLV <- diag(lamda_rmlv)
    A <- t(a[[t]]) %*% res[, t]
    S1 <- ((t(a[[t]]) %*% a[[t]]) + lamda_RMLV)
    beta_RMLV <- ginv(S1) %*% A
    g_pr <- a[[t]] %*% beta_RMLV
    
    cr <- cor(res[,t], g_pr)
    M <- mse(res[,t], g_pr)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- a[[t]] %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- res[,t] - g_pr
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
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

# Stop the parallel cluster
stopCluster(cl)

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV)
average_MSE_RMLV <- mean(mse_RMLV)
average_MAD_RMLV <- mean(mad_RMLV)

average_cr_RMLA <- mean(cor_RMLA)
average_MSE_RMLA <- mean(mse_RMLA)
average_MAD_RMLA <- mean(mad_RMLA)

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
################################### END OF ONE-STAGE ROBUST METHOD #########################################################


#################################################################################################
######### Two-stage non-robust method 
#################################################################################################

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
library(MASS)
library(ie2misc)
library(foreach)
library(doParallel)

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
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

results <- list()
train_data_list <- list()
test_data_list <- list()
train_a <- list()
res_train <- list()

for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Create a split: 70% training, 30% testing
    split <- sample.split(as[[t]]$y, SplitRatio = 0.7)
    
    # Subset the data into train and test, storing in lists
    train_data <- subset(as[[t]], split == TRUE)
    test_data <- subset(as[[t]], split == FALSE)
    
    res_train[[t]] <- subset(res[, t], split == TRUE)
    
    # Now split a[[t]] accordingly
    train_a[[t]] <- subset(a[[t]], split == TRUE)
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  lmer(y ~ 1 + (1 | `rep(genotype_values)`) + (1|rep)+(1|rep:block), train_data)
    
    mu = summary(fit)$coefficients[,1]
    
    # Extract mu from fit
    mu_values <- fitted(fit)
    
    results[[t]] <- list(mu = mu_values, fit=fit)
    
  }, error = function(e) {
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Now apply `mixed.solve` with the aligned `beta` and `a[[t]]`
    fit_2 <- mixed.solve(results[[t]]$mu, train_a[[t]], method="REML")
    
    
    beta=fit_2$u
    g_p=a[[t]]%*%beta
    
    c=cor(res[,t],g_p)
    M=mse(res[,t],g_p) ##RM11
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    h1 <- colSums(Z_replicated == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(Z_replicated == 1, na.rm = TRUE)  # Count of allele Aa
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(11640)
    k <- numeric(11640)
    l <- numeric(11640)
    s <- numeric(11640)
    su <- numeric(11640)
    sum <- numeric(11641)
    lamda_rmla <- numeric(11640)
    lamda_rr <- numeric(11640)
    dd_new <- data.frame(res[,t], a[[t]])
    
    for(i in 1:11640) {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (1430 - ((1 / 1430) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:11640) {
      lamda_rmla[j] <- lam * (sum[11641] / su[j])
      lamda_rr[j] <- (1 / 0.16 - 1) * 11640 * (sum[11641] / su[j])
      lamda_rmla[j]
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(a[[t]]) %*% a[[t]]
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- solve(S2) %*% t(a[[t]]) %*% res[, t]
    beta2 <- solve(S3) %*% t(a[[t]]) %*% res[, t]
    str(beta1)
    str(beta2)
    
    g_p1 <- a[[t]] %*% beta1  # For RRWA
    g_p2 <- a[[t]] %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(res[, t], g_p2)
    M2 <- mse(res[, t], g_p2)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- a[[t]] %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- res[,t] - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    u=as.matrix(fit_2$u)
    lam=(fit_2$Ve/fit_2$Vu)
    D=(t(a[[t]])%*%a[[t]]) + ((diag(rep(lam,11640))))
    cvv=solve(D)
    tr=tr(cvv)
    sigma=0
    lamda_rmlv=0
    for(i in 1:11640)
    {
      sigma[i]=((u[i,]*u[i,])-(fit_2$Ve*t))/i
      lamda_rmlv[i]=(fit_2$Ve/sigma[i])
    }
    lamda_RMLV=diag(lamda_rmlv)
    A=t(train_a[[t]])%*%results[[t]]$mu
    S1=((t(a[[t]])%*%a[[t]]) + lamda_RMLV)
    beta_RMLV=solve(S1)%*%A
    g_pr=a[[t]]%*%beta_RMLV
    
    cr <- cor(res[, t], g_pr)
    M <- mse(res[, t], g_pr)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- a[[t]] %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- res[,t] - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Combine results for this iteration
    c(cr, M, c2, M2, mad_RMLA[t], mad_RMLV[t], heritability_RMLA_values[t], heritability_RMLV_values[t])
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
    
  }, error = function(e) {
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV)
average_MSE_RMLV <- mean(mse_RMLV)
average_MAD_RMLV <- mean(mad_RMLV)

average_cr_RMLA <- mean(cor_RMLA)
average_MSE_RMLA <- mean(mse_RMLA)
average_MAD_RMLA <- mean(mad_RMLA)

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

########################## END OF TWO-STAGE NON-ROBUST METHOD ################################################


#################################################################################################
######### Two-stage robust method 
#################################################################################################

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
library(MASS)
library(ie2misc)
library(foreach)
library(doParallel)

# Set up parallel backend
cl <- makeCluster(detectCores() - 1)  # Leave one core free
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

results <- list()
train_data_list <- list()
test_data_list <- list()
train_a <- list()
res_train <- list()

for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Create a split: 70% training, 30% testing
    split <- sample.split(as[[t]]$y, SplitRatio = 0.7)
    
    # Subset the data into train and test, storing in lists
    train_data <- subset(as[[t]], split == TRUE)
    test_data <- subset(as[[t]], split == FALSE)
    
    res_train[[t]] <- subset(res[, t], split == TRUE)
    
    # Now split a[[t]] accordingly
    train_a[[t]] <- subset(a[[t]], split == TRUE)
    
    ## Step 1: Fit a linear mixed model on the training data
    # (Use train_data instead of as[[t]])
    fit<-  lmer(y ~ 1 + (1 | `rep(genotype_values)`) + (1|rep)+(1|rep:block), train_data)
    
    mu = summary(fit)$coefficients[,1]
    
    # Extract mu from fit
    mu_values <- fitted(fit)
    
    results[[t]] <- list(mu = mu_values, fit=fit)
    
  }, error = function(e) {
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 


for(t in 1:num_iterations) {
  tryCatch({
    
    # Print progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    # Assuming as[[t]] is your data, use caTools to split into train and test sets
    set.seed(123*t)  # For reproducibility
    
    # Now apply `mixed.solve` with the aligned `beta` and `a[[t]]`
    fit_2 <- mixed.solve(results[[t]]$mu, train_a[[t]], method="REML")
    
    beta=fit_2$u
    g_p=a[[t]]%*%beta
    
    c=cor(res[,t],g_p)
    M=mse(res[,t],g_p) ##RM11
    
    ### Step 2: Calculate variance components on the training data ###
    Ve <- get_variance_residual(fit)
    Vu <- get_variance(fit)$var.intercept[[1]][[1]]
    lam <- Ve / Vu
    print(lam)
    
    # Step 3: Calculate heritability
    heritability <- Vu / (Vu + (Ve / r))  # Heritability estimation
    heritability_values[t] <- heritability
    heritability_difference_values[t] <- (heritability - 0.7)^2
    
    ### Step 3: Count allele occurrences in the genotype matrix ###
    h1 <- colSums(Z_replicated == 0, na.rm = TRUE)  # Count of allele aa and AA
    h2 <- colSums(Z_replicated == 1, na.rm = TRUE)  # Count of allele Aa
    
    ### Step 4: Initialize variables for the loop ###
    n <- numeric(11640)
    k <- numeric(11640)
    l <- numeric(11640)
    s <- numeric(11640)
    su <- numeric(11640)
    sum <- numeric(11641)
    lamda_rmla <- numeric(11640)
    lamda_rr <- numeric(11640)
    dd_new <- data.frame(res[,t], a[[t]])
    
    for(i in 1:11640) {
      n[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[1]]
      k[i] <- summary((aov(dd_new[, 1] ~ as.factor(dd_new[, i + 1]), data = dd_new)))[[1]][[3]][[2]]
      l[i] <- n[i] - k[i]
      s[i] <- (0.5) * (1430 - ((1 / 1430) * ((h1[i] * h1[i]) + (h2[i] * h2[i]))))
      su[i] <- l[i] / s[i]
      sum[i + 1] <- sum[i] + su[i]
    }
    
    for(j in 1:11640) {
      lamda_rmla[j] <- lam * (sum[11641] / su[j])
      lamda_rr[j] <- (1 / 0.16 - 1) * 11640 * (sum[11641] / su[j])
      lamda_rmla[j]
    }
    
    # Continue with the rest of the code
    lamda_RR <- diag(lamda_rr)
    lamda_RRLA <- diag(lamda_rmla)
    str(lamda_RR)
    str(lamda_RRLA)
    
    ### Step 7: Solve for SNP effects (beta1 and beta2)###
    S <- t(a[[t]]) %*% a[[t]]
    S2 <- (S + lamda_RR)
    S3 <- (S + lamda_RRLA)
    
    beta1 <- solve(S2) %*% t(a[[t]]) %*% res[, t]
    beta2 <- solve(S3) %*% t(a[[t]]) %*% res[, t]
    str(beta1)
    str(beta2)
    
    g_p1 <- a[[t]] %*% beta1  # For RRWA
    g_p2 <- a[[t]] %*% beta2
    
    ### Step 8: Calculate performance metrics (Correlation and MSE) ###
    c2 <- cor(res[, t], g_p2)
    M2 <- mse(res[, t], g_p2)
    
    # Calculate MAD for RMLA
    mad_RMLA[t] <- madstat(g_p2)
    
    # Calculate predicted values for RMLA
    g_p2 <- a[[t]] %*% beta2
    
    # Estimate residuals and variance components
    residuals_RMLA <- res[,t] - g_p2
    Ve_RMLA <- var(residuals_RMLA)  # Residual variance
    Vu_RMLA <- var(g_p2)  # Genetic variance for RMLA
    
    # Calculate heritability for RMLA
    heritability_RMLA <- Vu_RMLA / (Vu_RMLA + (Ve_RMLA / r))
    heritability_difference_values_RMLA[t] <- (heritability_RMLA - 0.7)^2
    
    ###  Step 9: RMLV  ###
    u=as.matrix(fit_2$u)
    lam=(fit_2$Ve/fit_2$Vu)
    D=(t(a[[t]])%*%a[[t]]) + ((diag(rep(lam,11640))))
    cvv=solve(D)
    tr=tr(cvv)
    sigma=0
    lamda_rmlv=0
    for(i in 1:11640)
    {
      sigma[i]=((u[i,]*u[i,])-(fit_2$Ve*t))/i
      lamda_rmlv[i]=(fit_2$Ve/sigma[i])
    }
    lamda_RMLV=diag(lamda_rmlv)
    A=t(train_a[[t]])%*%results[[t]]$mu
    S1=((t(a[[t]])%*%a[[t]]) + lamda_RMLV)
    beta_RMLV=solve(S1)%*%A
    g_pr=a[[t]]%*%beta_RMLV
    
    cr <- cor(res[, t], g_pr)
    M <- mse(res[, t], g_pr)
    
    # Calculate MAD for RMLV
    mad_RMLV[t] <- madstat(g_pr)
    
    # Calculate predicted values for RMLV
    g_pr <- a[[t]] %*% beta_RMLV
    
    # Estimate residuals and variance components
    residuals_RMLV <- res[,t] - g_pr
    Ve_RMLV <- var(residuals_RMLV)  # Residual variance
    Vu_RMLV <- var(g_pr)  # Genetic variance for RMLV
    
    # Calculate heritability for RMLV
    heritability_RMLV <- Vu_RMLV / (Vu_RMLV + (Ve_RMLV / r))
    heritability_difference_values_RMLV[t] <- (heritability_RMLV - 0.7)^2
    
    # Combine results for this iteration
    c(cr, M, c2, M2, mad_RMLA[t], mad_RMLV[t], heritability_RMLA_values[t], heritability_RMLV_values[t])
    
    # Store results in vectors
    cor_RMLA[t] <- c2
    mse_RMLA[t] <- M2
    cor_RMLV[t] <- cr
    mse_RMLV[t] <- M
    heritability_RMLA_values[t] <- heritability_RMLA
    heritability_RMLV_values[t] <- heritability_RMLV
    
  }, error = function(e) {
    cat("Error in iteration", t, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

# After the loop ends, process the results further
average_cr_RMLV <- mean(cor_RMLV)
average_MSE_RMLV <- mean(mse_RMLV)
average_MAD_RMLV <- mean(mad_RMLV)

average_cr_RMLA <- mean(cor_RMLA)
average_MSE_RMLA <- mean(mse_RMLA)
average_MAD_RMLA <- mean(mad_RMLA)

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
################################### END OF TWO-STAGE ROBUST METHOD ################################################