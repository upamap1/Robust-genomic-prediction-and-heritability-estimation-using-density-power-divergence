##################################################################################
## Data simulation with no-contamination
##################################################################################

# Clear the workspace
rm(list = ls())

# Set working directory
setwd("C:/Users/RIKU/Documents")
getwd()

install.packages("LDlinkR")

# Load necessary packages
install.packages("devtools")  # Ensure devtools is installed
install.packages("AlphaSimR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
install.packages("simplePHENOTYPES")

# Replace the path below with the path to the downloaded package file on your desktop
package_path <- "C:/Users/RIKU/Desktop/hypred"

# Install the package
install.packages(package_path, repos = NULL, type = "source")
install.packages("predhy")

devtools::install_github("thierrygosselin/radiator")

# Ensure CJAMP is installed or remove if not available
# install.packages("CJAMP")  # Uncomment if CJAMP is available

# Load libraries
library(radiator)
library(AlphaSimR)
library(SNPRelate)
library(simplePHENOTYPES)
library(predhy)
library(lme4) 
library(robustlmm) 
library(psych) 
library(rrBLUP)
library(agricolae)
library(dplyr)
library(insight)
library(matrixStats)
library(Metrics)
library(LDlinkR)

as=list()
a = list()
bv=matrix(rep(0,16645200),nrow=1430,ncol=11640)
res=matrix(rep(0,16645200),nrow=1430,ncol=11640)

# Step 1: Simulate the Genetic Data with AlphaSimR############################
##############################################################################
# Loop for 100 iterations
for(m in 1:100) {
  tryCatch({
    
    set.seed(123 * m)  # Modify the seed with iteration number
    
    # Define Simulation Parameters
    nInd <- 715  # Number of individuals
    nChr <- 10    # Number of chromosomes
    totalSnps <- 11640  # Reduced total number of SNPs
    SNPsPerChr <- floor(totalSnps / nChr)  # SNP markers per chromosome
    
    ## Set number of replication
    r = 2
    
    # Overall mean
    phi <- 0.5
    
    # Define variance parameters
    var_s = 0.005892
    var_b = 6.3148
    var_e = 53.8715
    
    # Create founder population with genetic data
    founderPop <- runMacs(nInd = nInd, nChr = nChr, segSites = totalSnps)
    
    # Initialize simulation parameters
    SP <- SimParam$new(founderPop)
    SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)  # Add 1000 QTLs per chromosome
    SP$setVarE(H2 = 0.7)  # Set heritability
    SP$addSnpChip(nSnpPerChr = SNPsPerChr)  # Add reduced SNP markers
    
    # Generate Population Across Generations
    pop <- newPop(founderPop, simParam = SP)
    
    # Generate breeding values
    bvalues <- bv(pop, simParam = SP)
    
    
    # Step 2: Format Genotype Data & estimation of genotype values
    ###########################################################################################################
    
    # Extract genetic map
    gen_map <- getGenMap(SP)
    
    # Extract 
    Z <- pullSnpGeno(pop)
    str(Z)
    
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    Z[Z == 0] <- 0  # aa ->
    Z[Z == 1] <- 1   # Aa 
    Z[Z == 2] <- 0   # AA 
    
    # Replication of genotype matrix
    # **Important: Repeat genotype_values for replications**
    n <- nrow(Z)  # Get the current number of rows in Z
    num_replicates <- ceiling(1430 / n)
    Z_replicated <- Z[rep(seq_len(n), num_replicates), ]
    
    u = rnorm(totalSnps, 0, sqrt(var_s))  # Adjusted for total SNPs
    genotype_values = Z %*% u
    
  
    # Step 3: Randomized complete block design
    ###########################################################################################################
    
    # Create a randomized complete block design (RCBD)
    design <- design.rcbd(trt = 651:715, r = 2, seed = 123, serie = 2)  # Create the RCBD design
    book <- design$book
    bl <- c(rep(136, 13), rep(137, 13), rep(138, 13), rep(139, 13), rep(140, 13))
    data1 <- cbind(book, rep(1:13, 5))  # Bind the design with blocks
    str(data1)
    data1 <- data1[, c(1, 2, 3, 4)]  # Adjust columns as needed
    colnames(data1) <- c("plots", "rep", "geno", "block")  # Rename columns to match original
    
    trt_1 <- 1:650
    k1 = 10
    r1 = 2
    design_2 = design.alpha(trt_1, k1, r1, serie = 2, seed = 0, kinds = "Super-Duper", randomization = TRUE)
    data2 = design_2$book[, c(2, 3, 4, 5)]
    colnames(data2) <- c("plots", "block", "geno", "rep")
    
    # Convert the 'rep' column to numeric in both data frames
    data1$rep <- as.numeric(as.character(data1$rep))  # Convert to numeric
    data2$rep <- as.numeric(as.character(data2$rep))  # Convert to numeric
    
    # Convert the 'geno' column to character to avoid factor level issues
    data1$geno <- as.character(data1$geno)
    data2$geno <- as.character(data2$geno)
    
    # Now combine the two datasets
    firstdata <- rbind(data2, data1)  # Combine the two datasets
    
    # Remove the 'plots' column from the dataset
    firstdata <- firstdata[, -which(colnames(firstdata) == "plots")]
    
    
    # Step 4: creating the block effect
    ###########################################################################################################
    
    # **Important: Repeat genotype_values for replications**
    genotype_values <- rep(genotype_values, times = r)  # Now length is 1430
    
    # Generate block effects
    block_eff_1 <- rep(rnorm(5, 0, sqrt(var_b)), each = 13)  
    block_eff_2 <- rep(rnorm(65, 0, sqrt(var_b)), each = 10)   
    block_effects <- c(block_eff_1, block_eff_2)            
    
    # **Repeat block_effects for replications**
    block_effects <- rep(block_effects, times = r)  # Now length is 1430
    
    res_eff = rnorm(715 * 2, 0, sqrt(var_e))
    plot_eff = block_effects + res_eff
    
    finaldata = cbind(firstdata, block_effects, res_eff, plot_eff)
    
    dataset = (finaldata[order(finaldata$rep, finaldata$block), ])
    dataset = cbind(dataset, rep(genotype_values))
    dataset = cbind(dataset, rep(bvalues, 2))
    y = phi + dataset[, 6] + dataset[, 8]
    dataset = cbind(dataset, y)  ## pure dataset
    
    as[[m]] = dataset               # Store the current dataset in 'as'
    a[[m]] = Z_replicated           # Store the Z matrix for this iteration
    bv[, m] = dataset[1:1430, 8]    # Store breeding values in 'bv' for this iteration
    res[, m] = dataset[1:1430, 9]   # Store the response effect from the dataset
    
  }, error = function(e) {
    cat("Error in iteration", m, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

###################### END OF no-contamination code ################################################


####################################################################################################
## Data simulation with response-contamination
####################################################################################################

# Clear the workspace
rm(list = ls())

# Set working directory
setwd("C:/Users/RIKU/Documents")
getwd()

install.packages("LDlinkR")

# Load necessary packages
install.packages("devtools")  # Ensure devtools is installed
install.packages("AlphaSimR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
install.packages("simplePHENOTYPES")

# Replace the path below with the path to the downloaded package file on your desktop
package_path <- "C:/Users/RIKU/Desktop/hypred"

# Install the package
install.packages(package_path, repos = NULL, type = "source")
install.packages("predhy")

devtools::install_github("thierrygosselin/radiator")

# Ensure CJAMP is installed or remove if not available
# install.packages("CJAMP")  # Uncomment if CJAMP is available

# Load libraries
library(radiator)
library(AlphaSimR)
library(SNPRelate)
library(simplePHENOTYPES)
library(predhy)
library(lme4) 
library(robustlmm) 
library(psych) 
library(rrBLUP)
library(agricolae)
library(dplyr)
library(insight)
library(matrixStats)
library(Metrics)
library(LDlinkR)

as=list()
a = list()
bv=matrix(rep(0,16645200),nrow=1430,ncol=11640)
res=matrix(rep(0,16645200),nrow=1430,ncol=11640)

# Step 1: Simulate the Genetic Data with AlphaSimR############################
##############################################################################
# Loop for 100 iterations
for(m in 1:100) {
  tryCatch({
    
    set.seed(123 * m)  # Modify the seed with iteration number
    
    # Define Simulation Parameters
    nInd <- 715  # Number of individuals
    nChr <- 10    # Number of chromosomes
    totalSnps <- 11640  # Reduced total number of SNPs
    SNPsPerChr <- floor(totalSnps / nChr)  # SNP markers per chromosome
    
    ## Set number of replication
    r = 2
    
    # Overall mean
    phi <- 0.5
    
    # Define variance parameters
    var_s = 0.005892
    var_b = 6.3148
    var_e = 53.8715
    
    # Create founder population with genetic data
    founderPop <- runMacs(nInd = nInd, nChr = nChr, segSites = totalSnps)
    
    # Initialize simulation parameters
    SP <- SimParam$new(founderPop)
    SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)  # Add 1000 QTLs per chromosome
    SP$setVarE(H2 = 0.7)  # Set heritability
    SP$addSnpChip(nSnpPerChr = SNPsPerChr)  # Add reduced SNP markers
    
    # Generate Population Across Generations
    pop <- newPop(founderPop, simParam = SP)
    
    # Generate breeding values
    bvalues <- bv(pop, simParam = SP)
    
    # Step 2: Format Genotype Data & estimation of genotype values
    ###########################################################################################################
    
    # Extract genetic map
    gen_map <- getGenMap(SP)
    
    # Extract 
    Z <- pullSnpGeno(pop)
    str(Z)
    
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    Z[Z == 0] <- 0  # aa ->
    Z[Z == 1] <- 1   # Aa 
    Z[Z == 2] <- 0   # AA 
    
    # Replication of genotype matrix
    # **Important: Repeat genotype_values for replications**
    n <- nrow(Z)  # Get the current number of rows in Z
    num_replicates <- ceiling(1430 / n)
    Z_replicated <- Z[rep(seq_len(n), num_replicates), ]
    
    u = rnorm(totalSnps, 0, sqrt(var_s))  # Adjusted for total SNPs
    genotype_values = Z %*% u
    
   
    # Step 3: Randomized complete block design
    ###########################################################################################################
    
    # Create a randomized complete block design (RCBD)
    design <- design.rcbd(trt = 651:715, r = 2, seed = 123, serie = 2)  # Create the RCBD design
    book <- design$book
    bl <- c(rep(136, 13), rep(137, 13), rep(138, 13), rep(139, 13), rep(140, 13))
    data1 <- cbind(book, rep(1:13, 5))  # Bind the design with blocks
    str(data1)
    data1 <- data1[, c(1, 2, 3, 4)]  # Adjust columns as needed
    colnames(data1) <- c("plots", "rep", "geno", "block")  # Rename columns to match original
    
    trt_1 <- 1:650
    k1 = 10
    r1 = 2
    design_2 = design.alpha(trt_1, k1, r1, serie = 2, seed = 0, kinds = "Super-Duper", randomization = TRUE)
    data2 = design_2$book[, c(2, 3, 4, 5)]
    colnames(data2) <- c("plots", "block", "geno", "rep")
    
    # Convert the 'rep' column to numeric in both data frames
    data1$rep <- as.numeric(as.character(data1$rep))  # Convert to numeric
    data2$rep <- as.numeric(as.character(data2$rep))  # Convert to numeric
    
    # Convert the 'geno' column to character to avoid factor level issues
    data1$geno <- as.character(data1$geno)
    data2$geno <- as.character(data2$geno)
    
    # Now combine the two datasets
    firstdata <- rbind(data2, data1)  # Combine the two datasets
    
    # Remove the 'plots' column from the dataset
    firstdata <- firstdata[, -which(colnames(firstdata) == "plots")]
    
    
    # Step 4: creating the block effect
    ###########################################################################################################
    
    # **Important: Repeat genotype_values for replications**
    genotype_values <- rep(genotype_values, times = r)  # Now length is 1430
    
    # Generate block effects
    block_eff_1 <- rep(rnorm(5, 0, sqrt(var_b)), each = 13)  
    block_eff_2 <- rep(rnorm(65, 0, sqrt(var_b)), each = 10)   
    block_effects <- c(block_eff_1, block_eff_2)            
    
    # **Repeat block_effects for replications**
    block_effects <- rep(block_effects, times = r)  # Now length is 1430
    
    res_eff = rnorm(715 * 2, 0, sqrt(var_e))
    plot_eff = block_effects + res_eff
    
    finaldata = cbind(firstdata, block_effects, res_eff, plot_eff)
    
    dataset = (finaldata[order(finaldata$rep, finaldata$block), ])
    dataset = cbind(dataset, rep(genotype_values))
    dataset = cbind(dataset, rep(bvalues, 2))
    y = phi + dataset[, 6] + dataset[, 8]
    dataset = cbind(dataset, y)  ## pure dataset
    
    # Step 5: Introducing response contamination & creating the final dataset
    ###########################################################################################################
    
    # Add response contamination here
    perc.cont = 0.05
    n.outliers = round(1430 * perc.cont)  # Adjusted for the current dataset size
    outlier.positions = sample(1:1430, n.outliers)
    dataset[outlier.positions, 9] = dataset[outlier.positions, 9] + 5 * sqrt(var_e)
    
    as[[m]] = dataset               # Store the current dataset in 'as'
    a[[m]] = Z_replicated           # Store the Z matrix for this iteration
    bv[, m] = dataset[1:1430, 8]    # Store breeding values in 'bv' for this iteration
    res[, m] = dataset[1:1430, 9]   # Store the response effect from the dataset
    
  }, error = function(e) {
    cat("Error in iteration", m, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
}

############################### END of response contaminated code ##########################################


############################################################################################################
## Data simulation with block-contamination
############################################################################################################
# Clear the workspace
rm(list = ls())

# Set working directory
setwd("C:/Users/RIKU/Documents")
getwd()

install.packages("LDlinkR")

# Load necessary packages
install.packages("devtools")  # Ensure devtools is installed
install.packages("AlphaSimR")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SNPRelate")
install.packages("simplePHENOTYPES")

# Replace the path below with the path to the downloaded package file on your desktop
package_path <- "C:/Users/RIKU/Desktop/hypred"

# Install the package
install.packages(package_path, repos = NULL, type = "source")
install.packages("predhy")

devtools::install_github("thierrygosselin/radiator")

# Ensure CJAMP is installed or remove if not available
# install.packages("CJAMP")  # Uncomment if CJAMP is available

# Load libraries
library(radiator)
library(AlphaSimR)
library(SNPRelate)
library(simplePHENOTYPES)
library(predhy)
library(lme4) 
library(robustlmm) 
library(psych) 
library(rrBLUP)
library(agricolae)
library(dplyr)
library(insight)
library(matrixStats)
library(Metrics)
library(LDlinkR)

as=list()
a = list()
bv=matrix(rep(0,16645200),nrow=1430,ncol=11640)
res=matrix(rep(0,16645200),nrow=1430,ncol=11640)


# Step 1: Simulate the Genetic Data with AlphaSimR
######################################################################################
# Loop for 100 iterations
for(m in 1:100) {
  tryCatch({
    
    set.seed(123 * m)  # Modify the seed with iteration number
    
    # Define Simulation Parameters
    nInd <- 715  # Number of individuals
    nChr <- 10    # Number of chromosomes
    totalSnps <- 11640  # Reduced total number of SNPs
    SNPsPerChr <- floor(totalSnps / nChr)  # SNP markers per chromosome
    
    ## Set number of replication
    r = 2
    
    # Overall mean
    phi <- 0.5
    
    # Define variance parameters
    var_s = 0.005892
    var_b = 6.3148
    var_e = 53.8715
    
    # Create founder population with genetic data
    founderPop <- runMacs(nInd = nInd, nChr = nChr, segSites = totalSnps)
    
    # Initialize simulation parameters
    SP <- SimParam$new(founderPop)
    SP$addTraitA(nQtlPerChr = 1000, mean = 0, var = 1)  # Add 1000 QTLs per chromosome
    SP$setVarE(H2 = 0.7)  # Set heritability
    SP$addSnpChip(nSnpPerChr = SNPsPerChr)  # Add reduced SNP markers
    
    # Generate Population Across Generations
    pop <- newPop(founderPop, simParam = SP)
    
    # Generate breeding values
    bvalues <- bv(pop, simParam = SP)
    
    # Step 2: Format Genotype Data & estimation of genotype value
    ###########################################################################################################
    # Extract genetic map
    gen_map <- getGenMap(SP)
    
    # Extract 
    Z <- pullSnpGeno(pop)
    str(Z)
    
    # Convert the genotype coding from 0, 1, 2 to 1 and 0
    Z[Z == 0] <- 0  # aa ->
    Z[Z == 1] <- 1   # Aa 
    Z[Z == 2] <- 0   # AA 
    
    # Replication of genotype matrix
    # **Important: Repeat genotype_values for replications**
    n <- nrow(Z)  # Get the current number of rows in Z
    num_replicates <- ceiling(1430 / n)
    Z_replicated <- Z[rep(seq_len(n), num_replicates), ]
    
    u = rnorm(totalSnps, 0, sqrt(var_s))  # Adjusted for total SNPs
    genotype_values = Z %*% u
    
    
    # Step 3: Randomized complete block design
    ###########################################################################################################
    # Create a randomized complete block design (RCBD)
    design <- design.rcbd(trt = 651:715, r = 2, seed = 123, serie = 2)  # Create the RCBD design
    book <- design$book
    bl <- c(rep(136, 13), rep(137, 13), rep(138, 13), rep(139, 13), rep(140, 13))
    data1 <- cbind(book, rep(1:13, 5))  # Bind the design with blocks
    str(data1)
    data1 <- data1[, c(1, 2, 3, 4)]  # Adjust columns as needed
    colnames(data1) <- c("plots", "rep", "geno", "block")  # Rename columns to match original
    
    trt_1 <- 1:650
    k1 = 10
    r1 = 2
    design_2 = design.alpha(trt_1, k1, r1, serie = 2, seed = 0, kinds = "Super-Duper", randomization = TRUE)
    data2 = design_2$book[, c(2, 3, 4, 5)]
    colnames(data2) <- c("plots", "block", "geno", "rep")
    
    # Convert the 'rep' column to numeric in both data frames
    data1$rep <- as.numeric(as.character(data1$rep))  # Convert to numeric
    data2$rep <- as.numeric(as.character(data2$rep))  # Convert to numeric
    
    # Convert the 'geno' column to character to avoid factor level issues
    data1$geno <- as.character(data1$geno)
    data2$geno <- as.character(data2$geno)
    
    # Now combine the two datasets
    firstdata <- rbind(data2, data1)  # Combine the two datasets
    
    # Remove the 'plots' column from the dataset
    firstdata <- firstdata[, -which(colnames(firstdata) == "plots")]
    
    
    # Step 4: creating the block effect
    ###########################################################################################################
    
    # **Important: Repeat genotype_values for replications**
    genotype_values <- rep(genotype_values, times = r)  # Now length is 1430
    
    # Generate block effects
    block_eff_1 <- rep(rnorm(5, 0, sqrt(var_b)), each = 13)  
    block_eff_2 <- rep(rnorm(65, 0, sqrt(var_b)), each = 10)   
    block_effects <- c(block_eff_1, block_eff_2)       
    
    # Step 5: Introducing block contamination & creating the final dataset
    ###########################################################################################################
    
    # **Repeat block_effects for replications**
    block_effects <- rep(block_effects, times = r)  # Now length is 1430
    
    res_eff = rnorm(715 * 2, 0, sqrt(var_e))
    plot_eff = block_effects + res_eff
    
    finaldata = cbind(firstdata, block_effects, res_eff, plot_eff)
    
    dataset = (finaldata[order(finaldata$rep, finaldata$block), ])
    dataset = cbind(dataset, rep(genotype_values))
    dataset = cbind(dataset, rep(bvalues, 2))
    y = phi + dataset[, 6] + dataset[, 8]
    dataset = cbind(dataset, y)  ## pure dataset
    
    # Add block contamination after creating `dataset`
    ## contaminated block for replication 1
    n.blocks <- 3
    out.blocks <- sort(sample(1:135, n.blocks))
    for(i in 1:n.blocks) {
      dataset[dataset$block == out.blocks[i] & dataset$rep == 1, 9] <- 
        dataset[dataset$block == out.blocks[i] & dataset$rep == 1, 9] + 8 * sqrt(var_e)
    }
    
    # Additional contaminated blocks for replication 1 (other block range)
    n <- 2
    out.blocks2 <- sort(sample(136:139, n))
    for(i in 1:n) {
      dataset[dataset$block == out.blocks2[i] & dataset$rep == 1, 9] <- 
        dataset[dataset$block == out.blocks2[i] & dataset$rep == 1, 9] + 8 * sqrt(var_e)
    }
    
    ## contaminated block for replication 2
    n.blocks2 <- 3
    out.blocks3 <- sort(sample(1:135, n.blocks2))
    for(i in 1:n.blocks2) {
      dataset[dataset$block == out.blocks3[i] & dataset$rep == 2, 9] <- 
        dataset[dataset$block == out.blocks3[i] & dataset$rep == 2, 9] + 8 * sqrt(var_e)
    }
    
    # Additional contaminated blocks for replication 2 (other block range)
    n2 <- 2
    out.blocks4 <- sort(sample(136:139, n2))
    for(i in 1:n2) {
      dataset[dataset$block == out.blocks4[i] & dataset$rep == 2, 9] <- 
        dataset[dataset$block == out.blocks4[i] & dataset$rep == 2, 9] + 8 * sqrt(var_e)
    }
    
    as[[m]] = dataset               # Store the current dataset in 'as'
    a[[m]] = Z_replicated           # Store the Z matrix for this iteration
    bv[, m] = dataset[1:1430, 8]    # Store breeding values in 'bv' for this iteration
    res[, m] = dataset[1:1430, 9]   # Store the response effect from the dataset
    
  }, error = function(e) {
    cat("Error in iteration", m, ":", conditionMessage(e), "\n")
  })  # This closes tryCatch
  
} 

################################ END of block-contaminated data simulation #######################################