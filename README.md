# ProMOTe

## Overview

This project involves the application of Variational Bayes (VB) methods
to fit a probabilistic model to the presence and onset times of various
LTCs within a study population. The data is structured in an .rds file
containing a named list with multiple components, each representing
different aspects of the dataset.

## Dependencies

To run the code in this project, the following R libraries are required:

 - `extraDistr`: For working with additional distributions not available in
base R. 

## Data Structure

The `.rds` file contains a named list with the following components:

- `d`: An N x M matrix containing data about the presence of conditions for each individual. Each element represents whether a specific condition is present (1 for presence, 0 for absence).
- `t`: An N x M matrix containing data about the onset times of conditions for each individual. Each element represents the time at which the condition was observed or diagnosed.
- `rho`: A vector of length N containing the study start dates for each individual. Each element represents the start date of the study for a corresponding individual.
- `tau`: A vector of length N containing the study end dates for each individual. Each element represents the end date of the study for a corresponding individual.
- `iota`: A vector of length N indicating the status of each individual at the time of tau. Each element is a binary indicator, where 1 represents that the individual is alive and 0 represents that the individual is deceased.
- `N`: An integer representing the number of individuals in the dataset.
- `K`: An integer representing the number of clusters in the dataset. The clustering might be based on any relevant characteristic or analysis done on the dataset.
- `M`: An integer representing the number of conditions in the dataset.
- `sex`: A vector representing the sex of the individuals in the dataset. This is set to `NULL` if the information is not available.
- `birth_conds`: A vector containing the column indices in `d` that correspond to conditions which only occur at birth. This is set to `NULL` if there are no such conditions or the information is not available.
- `male_conds`: A vector containing the column indices in `d` that correspond to conditions which only occur in males. This is set to `NULL` if there are no such conditions or the information is not available.
- `female_conds`: A vector containing the column indices in `d` that correspond to conditions which only occur in females. This is set to `NULL` if there are no such conditions or the information is not available.
- `cond_list`: A string vector containing the names of the conditions included in the dataset. Each element in the vector corresponds to a condition listed in the matrix `d`.

## Functions

### 1. `VB_gaussian_update`

**Purpose:**
Perform Variational Bayes (VB) updates for the Gaussian latent class model with fixed K and censored data.

**Inputs:**

- `d`: A N x M matrix containing data about the presence of conditions.
- `t`: A N x M matrix containing data about the onset times of conditions.
- `rho`: A vector of length N containing study start dates.
- `tau`: A vector of length N containing study end dates.
- `iota`: A vector of length N indicating if individuals are alive/deceased at time tau.
- `hyperparameters`: A list of hyperparameters of the prior.
- `initial_Cstar`: A N x K matrix containing an initial value to initialize the latent variable `z`.
- `initial_Dstar`: A N x M matrix containing an initial value to initialize the latent variable `d`.
- `initial_pstar`: A N x M matrix containing an initial value to initialize the latent variable `d`.
- `initial_qstar`: A N x M matrix containing an initial value to initialize the latent variable `t`.
- `initial_rstar`: A N x M matrix containing an initial value to initialize the latent variable `t`.
- `N`: The number of individuals in the data.
- `K`: The number of clusters.
- `M`: The number of conditions in the data.
- `epsilon`: A number used to determine the stopping condition.
- `sex`: The sex of individuals in the data.
- `birth_conds`: Column indices for `d` indicating conditions that only occur at birth.
- `male_conds`: Column indices for `d` indicating conditions which only occur in males.
- `female_conds`: Column indices for `d` indicating conditions which only occur in females.
- `cond_list`: A string vector containing the names of the conditions.

**Outputs:**

- `theta_star`: A vector of parameters for the VB posterior of `gamma`.
- `a_star`: A M x K matrix of parameters for the VB posterior of `pi`.
- `b_star`: A M x K matrix of parameters for the VB posterior of `pi`.
- `u_star`, `v_star`, `alpha_star`, `beta_star`: M x K matrices of parameters for the VB posterior of `mu` and `sigma^2`.
- `C_star`: A N x K matrix of parameters for the VB posterior of `z`.
- `p_star`, `q_star`, `r_star`, `D_star`: N x M matrices of parameters for the VB posterior of `d` and `t`.
- `n_steps`: The number of iterations required to achieve the stopping condition.

### 2. `expected_lst_lefttrunc`

**Purpose:**
Calculate the expected value of a left-truncated location-scale t-distribution.

**Inputs:**

- `df`: Degrees of freedom of the t-distribution.
- `mu`: Location parameter of the t-distribution.
- `sigma`: Scale parameter of the t-distribution.
- `tau`: The left truncation point.

**Outputs:**

- `u`: The expected value.


### 3. `VB_gaussian_predictive_density`

**Purpose:**
Estimate the parameters of the posterior predictive distribution for a new individual.

**Inputs:**

- `hyperparameters`: A list of hyperparameters of the posterior.
- `M_obs`: Indices of the individual's fully observed conditions.
- `M_part`: Indices of the individual's partially observed conditions.
- `M_unobs`: Indices of the individual's unobserved conditions.
- `d_obs`: A vector containing the absence/presence of the fully observed conditions as 0s and 1s.
- `t_obs`: A vector containing the onset times of the fully observed conditions.
- `d_part`: A vector containing the absence/presence of the partially observed conditions as 0s and 1s.
- `rho`: A vector of length N containing observation start dates.
- `tau`: A vector of length N containing observation end dates.
- `M`: The number of conditions in the data.

**Outputs:**

- `phi`: A vector of cluster probabilities.
- `eta`: A matrix of condition probabilities given cluster.
- `varpi`: A vector of useful intermediate condition probabilities conditional on clusters.

### 4. `expected_LTHC_t_after_tau`

**Purpose:**
Estimate the expected time of occurrence of conditions in an individual after tau.

**Inputs:**

- `parameters`: A list of posterior predictive parameters (the output of `VB_gaussian_predictive_density`).
- `hyperparameters`: A list of hyperparameters of the posterior.
- `tau`: The individual's age at end of observation.
- `M`: The number of conditions.

**Outputs:**

- `Et`: Expected onset times.

### 5. `probability_LTHC_by_T`

**Purpose:**
Estimate the probability of conditions occurring in an individual by time T.

**Inputs:**

- `parameters`: A list of posterior predictive parameters (the output of `VB_gaussian_predictive_density`).
- `hyperparameters`: A list of hyperparameters of the posterior.
- `T`: The time by which the conditions should occur.
- `tau`: The individual's age at end of observation.
- `M`: The number of conditions.

**Outputs:**

- `prob`: A vector of probabilities.
