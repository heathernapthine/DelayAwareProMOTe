require(extraDistr)

# Input and Output Filenames
input <- "../../data/demodata.rds"
output <- "../../data/demooutput.rds"

# Read training data and source code
data_train_list <- readRDS(input)
source("../ProMOTe_VB.R")

# The data and number of clusters:
#   d - a NxM matrix containing data about presence of conditions
#   t - a NxM matrix containing data about onset times of conditions
#   rho - a length N vector containing study start dates
#   tau - a length N vector containing study end dates
#   iota - a length N vector indicating if individuals are alive/deceased at time tau
#   N - the number of individuals in the data
#   K - the number of clusters
#   M - the number of conditions in the data
#   sex - the sex of individuals in the data - set to NULL if not available
#   birth_conds - column indices for d indicating conditions which only occur at birth - set to NULL if not available
#   male_conds - column indices for d indicating conditions which only occur in males
#   female_conds - column indices for d indicating conditions which only occur in females
#   cond_list - a string vector containing the names of the conditions
d <- data_train_list$d
t <- data_train_list$t
rho <- data_train_list$rho
tau <- data_train_list$tau
iota <- data_train_list$iota
N <- data_train_list$N
M <- data_train_list$M
K <- 6
epsilon <- 0.1
sex <- data_train_list$sex
birth_conds <- data_train_list$birth_conds
male_conds <- data_train_list$male_conds
female_conds <- data_train_list$female_conds
cond_list <- data_train_list$disease_list

# Hyperparameters for the model
theta <- rep(1, K)
a <- matrix(1, M, K)
b <- matrix(1, M, K)
u <- matrix(50, M, K)
v <- matrix(0.3, M, K)
alpha <- matrix(5, M, K)
beta <- matrix(750, M, K)
hyperparameters <- list(theta, a, b, u, v, alpha, beta)

# Initial values for the Mean Field Variational Bayes(MFVB)
initial_Cstar <- matrix(runif(N*K, 0, 0.02), N, K)
initial_Dstar <- matrix(runif(N*M, 0, 1), N, M)
initial_pstar <- matrix(runif(N*M, 0, 10), N, M)
initial_qstar <- matrix(runif(N*M, 1,2), N, M)
initial_rstar <- matrix(runif(N*M, 0.01, 0.02), N, M)

# Fit the posterior via MFVB
posterior <- VB_gaussian_update(d, t, rho, tau, iota, hyperparameters, initial_Cstar, initial_Dstar, initial_pstar, initial_qstar, initial_rstar, N, M, K, epsilon, sex, birth_conds, male_conds, female_conds, cond_list)

# Save the posterior
saveRDS(posterior, output)
cat("\n")
