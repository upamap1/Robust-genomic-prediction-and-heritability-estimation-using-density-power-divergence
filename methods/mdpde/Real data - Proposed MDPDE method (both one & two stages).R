#####################################################################################
##### Proposed one stage MDPDE method
#####################################################################################


### A: Loading the necessary packages and dataset###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Documents")
getwd()

rm(list = ls())

# Installing necessary libraries
install.packages("nlme")
install.packages("rrBLUP")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)
library(caret)
library(glmnet)
library(rrBLUP)
library(ie2misc)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Initialize vectors to store results for 100 iterations
num_iterations <- 100  # Number of iterations
cor_values <- numeric(num_iterations)
mse_values <- numeric(num_iterations)
spearman_corr <- numeric(num_iterations)
mad_values <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)

for (t in 1:num_iterations) {
  tryCatch({
    
    # Print the progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    dataset= training_dataset
    colSums(dataset==2)
    dim(dataset)
    
    #### FItting the model ###################################################
    library(robustlmm)
    fit<-  rlmer(y ~ -1 + as.factor(geno) + (1|rep)+(1|rep:sets),dataset ,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    mu <- summary(fit)$coefficients[,1]
    
    Z_new=rbind(X,X)
    Z_t= Z_new[training_indices,]
    X_t= Z_t[1:185,]
    
    x=matrix(rep(1,185),nrow=185) 
    Y=as.matrix(mu)
    initial=c(5,0.0012,0.12)
    
    x=matrix(rep(1,528),nrow=264,ncol=2)
    alpha=0.9
    
    M=df[which(df$rep==1),]
    A=df[which(df$rep==2),]
    B=cbind(M$block,A$block)
    Y=cbind(M$y,A$y)
    L=as.vector(B)
    
    P=matrix(rep(0,264*2*10),nrow=264*2,ncol=10)
    for(i in 1:528)
    {
      P[i,L[i]]=1
    }
    #########################################################################
    
    ##Dpd function for one stage model analysis###############################
    # Set alpha (robustness parameter)
    alpha = 1  # You can explore 0.1, 0.3, 0.5, 0.7 later
    
    # Define initial parameter guesses
    initial=c(1.14,0.0015,0.055,1.5)
    
    dpd_function <- function(X, Y, Z, P, alpha, initial) {
      
      n= nrow(Y)
      s= ncol(Z)
      u= (nrow(P)/2)
      r= ncol(Y)
      
      # Split P into two matrices for block effects
      row_odd <- seq_len(nrow(P)) %% 2
      P1 <- P[row_odd == 1, ]  
      P2 <- P[row_odd == 0, ]
      
      # Function to compute the objective function
      m_function <- function(theta) {
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
    #################################################################
    
    #### DPD one stage analysis #####################################
    initial=c(1.14,0.0015,0.055,1.5)
    
    d=dpd_function(x,Y,X,P,alpha,initial) ##dpd with one stage modelling
    
    lamda=d[3]/d[2]
    l=diag(c(rep(lamda,1135)))
    Z_1=rbind(X,X)
    
    S=t(Z_1)%*%Z_1
    S2=(S+l)
    
    beta=solve(S2)%*%t(Z_1)%*%df$y
    X_1=matrix(rep(1,528),nrow=528)
    Y_p=X_1*d[1]+Z_1%*%beta
    
    c=cor(df$y,Y_p)
    M=mse(df$y,Y_p)
    corr = cor.test(df$y, Y_p, method = 'spearman')
    sp=corr$estimate
    mad=median(abs(df$y-Y_p))
    hr= d[2]/(d[2] + (d[3]/2))
    heritability_values[t] <- hr
    heritability_difference_values[t] <- (hr - 0.7)^2
    
    #########################################################
    # Store results in vectors
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

################################ End of proposed one stage MDPDE method #########################################


#####################################################################################
##### Proposed two-stage MDPDE method
#####################################################################################

### A: Loading the necessary packages and dataset###
####################################################
# Setting directory
setwd("C:/Users/RIKU/Documents")
getwd()

rm(list = ls())

# Installing necessary libraries
install.packages("nlme")
install.packages("rrBLUP")
library(lme4)
library(nlme)
library(Metrics)
library(robustlmm)
library(rrBLUP)
library(psych)
library(matrixStats)
library(Metrics)
library(insight)
library(caTools)
library(caret)
library(glmnet)
library(rrBLUP)
library(ie2misc)

# Loading the datasets
load("C:/Users/RIKU/Desktop/dataCorn_SS.RData")
load("C:/Users/RIKU/Desktop/dataCorn_WW.RData")
df <- read.csv("C:/Users/RIKU/Desktop/data.csv")

# Initialize vectors to store results for 100 iterations
num_iterations <- 100  # Number of iterations
cor_values <- numeric(num_iterations)
mse_values <- numeric(num_iterations)
spearman_corr <- numeric(num_iterations)
mad_values <- numeric(num_iterations)
heritability_values <- numeric(num_iterations)
heritability_difference_values <- numeric(num_iterations)

for (t in 1:num_iterations) {
  tryCatch({
    
    # Print the progress
    print(paste("Iteration:", t, "out of", num_iterations))
    
    #Data splitting and generation of train & testing data
    training_indices <- createDataPartition(df$rep, p = 0.7, list = FALSE)
    training_dataset  <- df[training_indices, ]
    testing_dataset <- df[-training_indices, ]
    
    dataset= training_dataset
    colSums(dataset==2)
    dim(dataset)
    
    #### Fitting the model ###################################################
    library(robustlmm)
    fit<-  rlmer(y ~ -1 + as.factor(geno) + (1|rep)+(1|rep:sets),dataset ,
                 rho.sigma.e = psi2propII(smoothPsi, k = 2.28))
    mu <- summary(fit)$coefficients[,1]
    
    Z_new=rbind(X,X)
    Z_t= Z_new[training_indices,]
    X_t= Z_t[1:185,]
    Y=as.matrix(mu)
    x=matrix(rep(1,185),nrow=185) 
    initial=c(5,0.0012,0.12)
    
    M=df[which(df$rep==1),]
    A=df[which(df$rep==2),]
    B=cbind(M$block,A$block)
    Y=cbind(M$y,A$y)
    L=as.vector(B)
    
    P=matrix(rep(0,185*2*10),nrow=185*2,ncol=10)
    for(i in 1:370)
    {
      P[i,L[i]]=1
    }
    #########################################################################
    
    ##Dpd function for one stage model analysis###############################
    # Set alpha (robustness parameter)
    alpha = 0.1  # You can explore 0.1, 0.3, 0.5, 0.7 later
    
    # Define initial parameter guesses
    initial=c(5,0.00012,0.12)
    
    dpd_function_2=function(X,Y,Z_t,alpha,initial){
      n <- nrow(Z_t)
      m_function=function(theta)
      {
        beta <- theta[1]
        sigma.a <- theta[2]
        sigma <- theta[3]
        V <- sigma*diag(1,nrow = n) + Z_t%*%t(Z_t)*sigma.a
        r=0
        l=0
        k=0
        for (i in 1:n)
        {
          r[i]=exp((-alpha/2)*(1/V[i,i])*((Y[i,]-X[i,]*beta)^2))
          l[i]=1/(((2*pi)^(alpha/2))*(V[i,i]^(alpha/2))*((alpha+1)^(1/2)))
          k[i]=l[i]-((1+(1/alpha))*r[i]*(1/(((2*pi)^(alpha/2))*(V[i,i]^(alpha/2)))))
        }
        f=(1/n)*sum(k)
      }
      Low <- c(-Inf,0,0)
      result <- optim(initial, m_function , gr = NULL,
                      method = "L-BFGS-B",
                      lower = Low, upper = Inf,
                      control = list(), hessian = FALSE)
      
      e=result$par
      
      return(e)
      
    }
    #################################################################
    
    #### DPD one stage analysis #####################################
    initial=c(5,0.00012,0.12)
    
    d=dpd_function_2(x,Y,X_t,alpha,initial) ## second stage of two stage model using DPD
    
    lamda=(d[3]^2)/(d[2]^2)
    l=diag(c(rep(lamda,1135)))
    Z_1=rbind(X,X)
    
    S=t(Z_1)%*%Z_1
    S2=(S+l)
    
    beta=solve(S2)%*%t(Z_1)%*%df$y
    
    g_p=matrix(rep(d[1],528),nrow=528) + (Z_new%*%beta)
    
    c=cor(df$y,g_p)## correlation coeffecient
    M=mse(df$y,g_p)  ## Mean square error
    corr = cor.test(df$y, g_p, method = 'spearman')
    sp=corr$estimate  ## spearman rank correlation
    mad=median(abs(df$y-g_p)) ## median absolute deviation
    hr= d[2]/(d[2] + (d[3]/2))
    heritability_values[t] <- hr
    heritability_difference_values[t] <- (hr - 0.7)^2
    
    #########################################################
    # Store results in vectors
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


##################################  End of proposed two-stage MDPDE method ######################################
