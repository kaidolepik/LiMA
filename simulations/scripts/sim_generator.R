suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(mvtnorm))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(matrixcalc))

##### GENERATE QTL NUMBERS OF MEDIATORS #####
generate_QTLs <- function(k, dist) {
    dist$n_eqtls[cut(runif(k), breaks = c(0, dist$cumulative_p), labels = FALSE)]
}

##### GENERATE MEDIATOR CORRELATION MATRIX SIGMA #####
generate_Sigma <- function(Sigma_type, k, in_Sigma_RDS = NULL) {
    # Could generate random Sigma as follows: 
    # Sigma <- cor(matrix(rnorm(k*20), nrow = 20, ncol = k)); Sigma <- 0.01*diag(k) + 0.99*Sigma
    # Instead we're either creating an identity Sigma or subsetting from the supplied RDS (e.g. Colaus gene expression Sigma)
    
    if (Sigma_type == "identity")
        Sigma <- diag(k)
    else {
        if (Sigma_type == "randomExtreme")
            Sigma <- cor(matrix(rnorm(k*10), nrow = 10, ncol = k))
        else if (Sigma_type == "Colaus") {
            Sigma <- readRDS(in_Sigma_RDS)
            ix <- sample(1:nrow(Sigma), k)
            Sigma <- Sigma[ix, ix]
            
            colnames(Sigma) <- NULL
            rownames(Sigma) <- NULL
        }
        
        lambda <- 1e-6
        Sigma <- lambda*diag(k) + (1-lambda)*Sigma
    }
    
    return(Sigma)
}

##### GENERATE 2x2 COVARIANCE MATRIX SIGMA FOR (gamma, delta); note (1-p_direct) and cov_mediated are related #####
generate_mediation_Sigma <- function(k, total_causal, p_direct, cor_mediated, var_explained_ZtoY) {
    E_var_explained_XtoZtoY <- ((1-p_direct) * total_causal)^2
    
    var_delta <- var_explained_ZtoY / k #sigmad2
    var_gamma <- ifelse(E_var_explained_XtoZtoY == 0, var_delta, E_var_explained_XtoZtoY / (k*var_delta + (k^2+k)*cor_mediated^2*var_delta)) #sigmag2
    cov_mediated <- cor_mediated * sqrt(var_delta) * sqrt(var_gamma) #sigmagd
    
    matrix(c(var_gamma, cov_mediated, cov_mediated, var_delta), nrow = 2, ncol = 2)
}

##### GENERATE CAUSAL EFFECTS (gamma, delta) MEDIATED BY Z FROM A BIVARIATE NORMAL (not scaled to total effect atm) #####
generate_mediated_effects <- function(k, total_causal, p_direct, cor_mediated, var_explained_ZtoY) {
    mediation_Sigma <- generate_mediation_Sigma(k, total_causal, p_direct, cor_mediated, var_explained_ZtoY)
    
    rmvnorm(k, mean = c(0, 0), sigma = mediation_Sigma)
}

##### GENERATE TRUE DATA ACCORDING TO THE DIRECTED ACYCLIC GRAPH #####
generate_true_data <- function(h2X, m, k, k_sig_perc, total_causal, p_direct, cor_mediated, var_explained_ZtoY, QTL_dist, seed = NULL) {
    set.seed(seed)
    
    k_sig <- k*k_sig_perc
    
    # Direct causal effects
    mediated_effects <- matrix(0, nrow = k, ncol = 2)
    if (k_sig > 0)
        mediated_effects[1:k_sig, ] <- generate_mediated_effects(k_sig, total_causal, p_direct, cor_mediated, var_explained_ZtoY)
    
    if (p_direct == 1) # Scaling for T1E by changing either gamma_i or delta_i to 0 for each i
        mediated_effects[cbind(1:k_sig, sample(1:2, k_sig, replace = TRUE))] <- 0
    # else { # This scales the total indirect effect appropriately to MP*total_causal
    #     const <- (1-p_direct)*total_causal / sum(mediated_effects[1:k_sig, 1] * mediated_effects[1:k_sig, 2])
    #     
    #     mediated_effects[1:k_sig, ] <- mediated_effects[1:k_sig, ] * sqrt(abs(const))
    #     mediated_effects[1:k_sig, 1] <- mediated_effects[1:k_sig, 1] * ifelse(sign(const) == sign(total_causal), 1, -1)
    # }
    
    gamma <- mediated_effects[, 1]
    delta <- mediated_effects[, 2]
    alpha <- ifelse(p_direct == 1, runif(1, -total_causal, 3*total_causal), sum(gamma*delta) / (1-p_direct) * p_direct)
    # If we're investigating T1E then we generate alpha uniformly around the expected total causal effect
    
    # Mediator heritabilities and QTL numbers
    h2Z <- generate_h2(k)
    QTLs <- generate_QTLs(k, QTL_dist)
    
    # Mediator instrument effects on the mediators and outcome, respectively
    B <- lapply(1:k, function(i) rnorm(QTLs[i], mean = 0, sd = sqrt(h2Z[i]/QTLs[i]))) %>%
        bdiag() %>%
        as.matrix() %>%
        (function(B) sweep(B, 2, sqrt(h2Z/colSums(B^2)), "*")) # Scale to h2Z
    b <- c(B %*% delta)
    
    # Exposure instrument effects on the exposure, mediators and outcome, respectively
    beta <- rnorm(m, mean = 0, sd = sqrt(h2X/m)) %>%
        (function(beta) beta * sqrt(h2X / sum(beta^2))) # Scale to h2X
    C <- beta %*% t(gamma)
    c <- beta*alpha + c(C %*% delta)
    
    list(gamma=gamma, delta=delta, alpha=alpha, B=B, b=b, C=C, c=c, beta=beta, QTLs=QTLs)
}

##### GENERATE EFFECT ESTIMATES BY ADDING RANDOM NOISE TO THE TRUE DATA #####
generate_estimate_noise <- function(B, b, C, c, beta, QTLs, nX, nZ, nY, k, m, sigmaC2, sigmac2, sigmab2, Sigma_chol) {
    # Add relevant noise to the effects
    beta_med_hat <- 0 + rnorm(sum(QTLs), mean = 0, sd = sqrt(1/nX))
    B_hat <- B + matrix(rnorm(sum(QTLs) * k), ncol = k) %*% (sqrt(1/nZ)*Sigma_chol)
    b_hat <- b + rnorm(sum(QTLs), mean = 0, sd = sqrt(1/nY + sigmab2))
    beta_hat <- beta + rnorm(m, mean = 0, sd = sqrt(1/nX)) 
    C_hat <- C + matrix(rnorm(m*k), ncol = k) %*% (sqrt(1/nZ + sigmaC2) * Sigma_chol)
    c_hat <- c + rnorm(m, mean = 0, sd = sqrt(1/nY + sigmac2))
    
    list(beta_med_hat=beta_med_hat, B_hat=B_hat, b_hat=b_hat, beta_hat=beta_hat, C_hat=C_hat, c_hat=c_hat)
}

generate_estimates <- function(true_data, m, nX, nZ, nY, k, p_value_mediators, sigmaC2, 
                               sigmac2, sigmab2, Sigma_chol, seed = NULL) {
    set.seed(seed)
    
    # Generate estimates and select mediators for the simulation
    estimates <- list()
    ix_mediators <- c()
    while (length(ix_mediators) == 0) { # Require at least one causal exposure-mediator relationship
        estimates <- generate_estimate_noise(true_data$B, true_data$b, true_data$C, true_data$c, true_data$beta, true_data$QTLs, 
                                             nX, nZ, nY, k, m, sigmaC2, sigmac2, sigmab2, Sigma_chol)
        
        ix_mediators <- 1:k %>%
            sapply(function(i) ivw_mr(estimates$beta_hat, estimates$C_hat[, i], rep(1/nZ, m))$P <= p_value_mediators) %>%
            which()
        
        seed <- seed + 1
        set.seed(seed)
    }
    
    # Mediator instruments
    ix_instruments <- tibble(ix_mediator = rep(1:k, true_data$QTLs), ix_instrument = 1:sum(true_data$QTLs)) %>%
        filter(ix_mediator %in% ix_mediators) %>%
        pull(ix_instrument)
    
    list(estimates = estimates,
         ix_instruments = ix_instruments, 
         ix_mediators = ix_mediators)
}

generate_simulation <- function(m, nX, nZ, nY, h2X, total_causal, p_direct, cor_mediated, var_explained_ZtoY, k, k_sig_perc, p_value_mediators, 
                                sigmaC2, sigmac2, sigmab2, QTL_dist, Sigma_chol, generate_seed = NULL, estimate_seed = NULL) {
    # Generate true data for the simulation
    true_data <- generate_true_data(h2X, m, k, k_sig_perc, total_causal, p_direct, cor_mediated, var_explained_ZtoY, QTL_dist, seed = generate_seed)
    
    # Generate estimates for the true data
    estimates <- generate_estimates(true_data, m, nX, nZ, nY, k, p_value_mediators, sigmaC2, sigmac2, sigmab2, Sigma_chol, seed = estimate_seed)

    ix_instruments <- estimates$ix_instruments
    ix_mediators <- estimates$ix_mediators
    estimates <- estimates$estimates
    
    list(true_total = sum(true_data$gamma * true_data$delta) + true_data$alpha,
         true_alpha = true_data$alpha,
         true_gamma_nonzero = true_data$gamma[1:(k*k_sig_perc)],
         true_delta_nonzero = true_data$delta[1:(k*k_sig_perc)],
         ix_mediators = ix_mediators,
         generate_seed = generate_seed,
         estimate_seed = estimate_seed,
         k = length(ix_mediators),
         l = length(ix_instruments),
         m = m,
         nX = nX, 
         nY = nY, 
         nZ = nZ,
         B_hat = estimates$B_hat[ix_instruments, ix_mediators, drop = FALSE],
         C_hat = estimates$C_hat[, ix_mediators, drop = FALSE],
         b_hat = estimates$b_hat[ix_instruments],
         c_hat = estimates$c_hat,
         beta_hat = estimates$beta_hat,
         beta_med_hat = estimates$beta_med_hat[ix_instruments],
         Sigma = estimate_Sigma(estimates$B_hat[, ix_mediators, drop = FALSE]), # B_hat works better than C_hat or rbind(B_hat, C_hat), probably C_hat effects can be too small/random as they're not mediator instruments
         var_B_hat = matrix(1/nZ, nrow = length(ix_instruments), ncol = length(ix_mediators)),
         var_C_hat = matrix(1/nZ, nrow = m, ncol = length(ix_mediators)),
         var_b_hat = rep(1/nY, length(ix_instruments)),
         var_c_hat = rep(1/nY, m),
         var_beta_hat = rep(1/nX, m),
         var_beta_med_hat = rep(1/nX, length(ix_instruments)))
}

parser <- ArgumentParser()
parser$add_argument("--m", type = "integer")
parser$add_argument("--nX", type = "integer")
parser$add_argument("--nZ", type = "integer")
parser$add_argument("--nY", type = "integer")
parser$add_argument("--h2X", type = "double")
parser$add_argument("--total_causal", type = "double")
parser$add_argument("--p_direct", type = "double")
parser$add_argument("--cor_mediated", type = "double")
parser$add_argument("--var_explained_ZtoY", type = "double")
parser$add_argument("--k", type = "integer")
parser$add_argument("--k_sig_perc", type = "double")
parser$add_argument("--p_value_mediators", type = "double")
parser$add_argument("--Sigma_type")
parser$add_argument("--sigmaC2", type = "double")
parser$add_argument("--sigmac2", type = "double")
parser$add_argument("--sigmab2", type = "double")
parser$add_argument("--seed", type = "integer")
parser$add_argument("--start", type = "integer")
parser$add_argument("--end", type = "integer")
parser$add_argument("--heritability_tools_R")
parser$add_argument("--ML_methods_R")
parser$add_argument("--in_mediator_QTLs_txt")
parser$add_argument("--in_mediator_Sigma_RDS")
parser$add_argument("--out_simulation_RDS")
args <- parser$parse_args()

source(args$heritability_tools_R)
source(args$ML_methods_R)

QTL_dist <- fread(args$in_mediator_QTLs_txt) # Distribution of QTLs (based on eQTLs from the eQTLGen Consortium)

# Generate the mediator correlation matrix
set.seed(args$seed)
Sigma_chol <- generate_Sigma(args$Sigma_type, args$k, args$in_mediator_Sigma_RDS) %>% # Generate mediator Sigma using Cholesky decomposition
    chol()

sim_bundle <- list()
for (i in args$start:args$end) {
    sim_bundle[[length(sim_bundle) + 1]] <- generate_simulation(args$m, args$nX, args$nZ, args$nY, args$h2X, args$total_causal, args$p_direct, args$cor_mediated, args$var_explained_ZtoY, args$k, 
                                                                args$k_sig_perc, args$p_value_mediators, args$sigmaC2, args$sigmac2, args$sigmab2, QTL_dist, Sigma_chol, generate_seed = i*10, estimate_seed = i*1e5)
    cat("Simulation", i, "generated:", format(Sys.time(), usetz = TRUE), "\n")
}

saveRDS(sim_bundle, args$out_simulation_RDS)
