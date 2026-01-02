############################################################################################################
####### Proposed MDPDE one-stage method
############################################################################################################

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

# Set up parallel backend
cl <- makeCluster(detectCores() - 2)  # Leave one core free
registerDoParallel(cl)

####################################################################

# Initialize vectors to store results for 100 iterations
num_iterations <- 100  # Number of iterations
# Initialize vectors to store results
cor_values <- numeric(num_iterations)
mse_values <- numeric(num_iterations)
spearman_corr <- numeric(num_iterations)
mad_values <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)

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
    
    #Fitting the model 
    names(train_data)[names(train_data) == "rep(bvalues, 2)"] <- "bvalues_replicated"
    fit <- rlmer(y ~ 1 + (1 | `rep(genotype_values)`) + (1|rep) + (1|rep:block_effects), train_data,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    
    mu <- summary(fit)$coefficients[,1]
    
    X=matrix(rep(1,1430),nrow=1430,ncol=1) ## as intercept model is taken
    alpha=1 ## alpha have 4 more choices (0.1,0.3,0.5,0.7))
    
    M=as[[t]][which(train_data$rep==1),]
    A=as[[t]][which(train_data$rep==2),]
    B=cbind(M$block,A$block)
    Y=cbind(as[[t]]$y,as[[t]]$y)
    L=as.vector(B)
    
    P=matrix(rep(0,715*2*11640),nrow=715*2,ncol=11640)
    for(i in 1:1430)
    {
      P[i,L[i]]=1
    }
    
    Z=a[[t]]
    
    ##Dpd function for one stage model analysis
    
    initial=c(5,0.0012,0.12)  
    alpha=0.1
    pi <- base::pi
    
    dpd_function <- function(X, Y, Z, P, alpha, initial) {
      
      n= nrow(Y)
      s= ncol(Z)
      u= (nrow(P)/2)
      r= ncol(Y)
      
      # Split P into two matrices for block effects
      row_odd <- seq_len(nrow(P)) %% 2
      P1 <- P[row_odd == 1, ]  
      P2 <- P[row_odd == 0, ]
      
      m_function=function(theta)
      {
        beta <- theta[1]     # Coefficient for the fixed effect
        sigma.g <- theta[2]  # Variance for genotype effects
        sigma.e <- theta[3]  # Residual variance
        sigma.b <- theta[4]  # Variance for block effects
        
        # Initialize vectors to store results
        r_values <- numeric(n)
        l_values <- numeric(n)
        k_values <- numeric(n)
        
        # Loop through each observation
        for (i in 1:n) {
          # Correctly map i to access corresponding rows in P1 and P2
          p_index <- (i + 1) / 2  # Index for P1 and P2
          
          w <- rbind(Z[i, ], Z[i, ])  # Genotype effect for observation i
          Q <- rbind(P1[p_index, ], P2[p_index, ])  # Block effect for observation i (from P1 and P2)
          
          # Step 2: Construct V matrix
          V <- sigma.e * diag(1, nrow = 2) + (w %*% t(w)) * sigma.g + (Q %*% t(Q)) * sigma.b
          
          # Calculate residual vector, ensuring it remains a matrix
          residual_vector <- matrix(Y[i, ] - (X[i, ] * beta), nrow = 2)  # 1 x r
          
          # Now calculate r_i
          r_values[i] <- exp((-alpha / 2) * (t(residual_vector) %*% solve(V) %*% residual_vector))
          
          # Calculate l_i
          l_values[i] <- 1 / (((2 * pi) ^ alpha) * ((det(V)) ^ (alpha / 2)) * (alpha + 1))
          
          # Calculate k_i
          k_values[i] <- l_values[i] - ((1 + (1 / alpha)) * r_values[i] * (1 / (((2 * pi) ^ alpha) * ((det(V)) ^ (alpha / 2)))))
        }
        
        # Objective function value to minimize
        f <- (1 / n) * sum(k_values)
        
        return(f)  # Return the value of the objective function
      }
      
      result <- optim(initial, m_function, gr = NULL,
                      method = "BFGS",   # Try a different method
                      control = list(maxit = 10000, reltol = 1e-6),  # Increase iterations
                      hessian = TRUE)  # Hessian for better convergence
      
      
      return(result$par)  # Return estimated parameters
    }
    ############################################################
    
    ##one-stage model analysis with DPD
    d=dpd_function(X,Y,Z,P,alpha,initial) ## estimate of parameters using DPD
    
    lamda=d[3]/d[2]
    l=diag(c(rep(lamda,11640)))
    Z_1=rbind(a[[t]])
    
    S=t(Z_1)%*%Z_1
    S2=(S+l)
    
    beta=solve(S2)%*%t(Z_1)%*%bv[,t]
    X_1=matrix(rep(1,1430),nrow=1430)
    
    Y_p=X_1*d[1]+Z_1%*%beta
    c=cor(as[[t]]$y,Y_p)       ## correlation coeffecient
    M=mse(as[[t]]$y,Y_p)         ## mean square error
    corr = cor.test(as[[t]]$`rep(bvalues, 2)`, Y_p, method = 'spearman')
    sp=corr$estimate     ##spearman's rank correlation 
    mad=median(abs(as[[t]]$`rep(bvalues, 2)`-Y_p))  ##median absolute deviation
    
    ##Val=me_function(S,d[3],c) ## two different prediction measure (Pa_1 , heredity)
    hr= d[2]/(d[2] + (d[3]/2))
    heritability_values[t] <- hr
    heritability_difference_values[t] <- (hr - 0.7)^2
    
    
    cor_values[t] <- c
    mse_values[t] <- M
    spearman_corr[t] <- sp
    mad_values[t] <- mad
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
  })
}

# After the loop, compute averages
average_correlation <- mean(cor_values, na.rm = TRUE)
average_mse <- mean(mse_values, na.rm = TRUE)
average_spearman_corr <- mean(spearman_corr, na.rm = TRUE)
average_mad <- mean(mad_values, na.rm = TRUE)
average_heritability <- mean(heritability_values, na.rm = TRUE)
msd_heritability <- mean(heritability_difference_values, na.rm = TRUE)

# Print results
cat("Average Correlation:", average_correlation, "\n")
cat("Average MSE:", average_mse, "\n")
cat("Average Spearman Correlation:", average_spearman_corr, "\n")
cat("Average MAD:", average_mad, "\n")
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", msd_heritability, "\n")

##################################### END OF ONE-STAGE MDPDE METHOD ########################################




############################################################################################################
######## Proposed MDPDE two-stage method
############################################################################################################

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

# Set up parallel backend
cl <- makeCluster(detectCores() - 2)  # Leave one core free
registerDoParallel(cl)

####################################################################

# Initialize vectors to store results for 100 iterations
num_iterations <- 100  # Number of iterations

# Initialize vectors to store results
cor_values <- numeric(num_iterations)
mse_values <- numeric(num_iterations)
spearman_corr <- numeric(num_iterations)
mad_values <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)

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
    
    #Fitting the model 
    names(train_data)[names(train_data) == "rep(bvalues, 2)"] <- "bvalues_replicated"
    fit <- rlmer(y ~ 1 + (1 | bvalues_replicated) + (1|rep) + (1|rep:block_effects), train_data,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    
    mu <- summary(fit)$coefficients[,1]
    
    X=matrix(rep(1,1430),nrow=1430,ncol=1) ## as intercept model is taken
    
    
    M=as[[t]][which(train_data$rep==1),]
    A=as[[t]][which(train_data$rep==2),]
    B=cbind(M$block,A$block)
    Y=cbind(as[[t]]$y,as[[t]]$y)
    L=as.vector(B)
    
    pi <- base::pi
    
    P=matrix(rep(0,715*2*11640),nrow=715*2,ncol=11640)
    for(i in 1:1430)
    {
      P[i,L[i]]=1
    }
    
    Z=a[[t]]
    
    # Set alpha (robustness parameter)
    alpha = 0.1  # You can explore 0.1, 0.3, 0.5, 0.7 later
    
    # Define initial parameter guesses
    initial=c(5,0.00012,0.12)
    
    # Define the dpd_function_2
    dpd_function_2 <- function(X, Y, Z, alpha, initial) {
      n <- nrow(Z)
      
      # Define the objective function for optimization
      m_function <- function(theta) {
        beta <- theta[1]   # These values will be optimized by 'optim'
        sigma.a <- theta[2]
        sigma <- theta[3]
        
        V <- sigma * diag(1, n) + Z %*% t(Z) * sigma.a
        
        r <- numeric(n)  # Initialize vectors for the loop
        l <- numeric(n)
        k <- numeric(n)
        
        for (i in 1:n) {
          r[i] <- exp((-alpha / 2) * (1 / V[i, i]) * ((Y[i, ] - X[i, ] * beta)^2))
          l[i] <- 1 / (((2 * pi)^(alpha / 2)) * (V[i, i]^(alpha / 2)) * ((alpha + 1)^(1/2)))
          k[i] <- l[i] - ((1 + (1 / alpha)) * r[i] * (1 / (((2 * pi)^(alpha / 2)) * (V[i, i]^(alpha / 2)))))
        }
        
        f <- (1 / n) * sum(k)
        return(f)  # Optim expects a scalar value as output
      }
      
      # Set lower bounds for parameters in optim
      Low <- c(-Inf, 0, 0)
      
      # Run optimization
      result <- optim(par = initial, fn = m_function, method = "L-BFGS-B",
                      lower = Low, upper = Inf, control = list(), hessian = FALSE)
      
      return(result$par)  # Return optimized parameters
    }
    
    ############################################################
    
    ##one-stage model analysis with DPD
    d=dpd_function_2(X,Y,Z,alpha,initial) ## estimate of parameters using DPD
    
    lamda=d[3]/d[2]
    l=diag(c(rep(lamda,11640)))
    Z_1=rbind(a[[t]])
    
    S=t(Z_1)%*%Z_1
    S2=(S+l)
    
    beta=solve(S2)%*%t(Z_1)%*%as[[t]]$y
    X_1=matrix(rep(1,1430),nrow=1430)
    
    Y_p=X_1*d[1]+Z_1%*%beta
    c=cor(as[[t]]$`rep(bvalues, 2)`,Y_p)       ## correlation coeffecient
    M=mse(as[[t]]$`rep(bvalues, 2)`,Y_p)         ## mean square error
    corr = cor.test(as[[t]]$`rep(bvalues, 2)`, Y_p, method = 'spearman')
    sp=corr$estimate     ##spearman's rank correlation 
    mad=median(abs(as[[t]]$`rep(bvalues, 2)`-Y_p))  ##median absolute deviation
    
    ##Val=me_function(S,d[3],c) ## two different prediction measure (Pa_1 , heredity)
    hr= d[2]/(d[2] + (d[3]/2))
    heritability_values[t] <- hr
    heritability_difference_values[t] <- (hr - 0.7)^2    
    
    cor_values[t] <- c
    mse_values[t] <- M
    spearman_corr[t] <- sp
    mad_values[t] <- mad
    
  }, error = function(e) {
    message(sprintf("Error in iteration %d: %s", i, e$message))
  })
}

# After the loop, compute averages
average_correlation <- mean(cor_values, na.rm = TRUE)
average_mse <- mean(mse_values, na.rm = TRUE)
average_spearman_corr <- mean(spearman_corr, na.rm = TRUE)
average_mad <- mean(mad_values, na.rm = TRUE)
average_heritability <- mean(heritability_values, na.rm = TRUE)
msd_heritability <- mean(heritability_difference_values, na.rm = TRUE)

# Print results
cat("Average Correlation:", average_correlation, "\n")
cat("Average MSE:", average_mse, "\n")
cat("Average Spearman Correlation:", average_spearman_corr, "\n")
cat("Average MAD:", average_mad, "\n")
cat("Average Heritability:", average_heritability, "\n")
cat("Average MSD of Heritability:", msd_heritability, "\n")

####################################### END OF MDPDE TWO-STAGE METHOD ##############################################

