suppressPackageStartupMessages(library(argparse))

simulate <- function(sim, method, pleiotropy = FALSE) {
    B_hat <- sim$B_hat
    C_hat <- sim$C_hat
    b_hat <- sim$b_hat
    c_hat <- sim$c_hat
    beta_hat <- sim$beta_hat
    beta_med_hat <- sim$beta_med_hat
    Sigma <- sim$Sigma
    var_B_hat <- sim$var_B_hat
    var_C_hat <- sim$var_C_hat
    var_b_hat <- sim$var_b_hat
    var_c_hat <- sim$var_c_hat
    var_beta_hat <- sim$var_beta_hat
    var_beta_med_hat <- sim$var_beta_med_hat
    nX <- sim$nX
    nY <- sim$nY
    nZ <- sim$nZ
    k_selected <- sim$k
    l <- sim$l
    m <- sim$m
    
    par <- mediation_analysis(B_hat, C_hat, b_hat, c_hat, beta_hat, beta_med_hat, Sigma, var_B_hat, var_C_hat, var_b_hat, 
                              var_c_hat, var_beta_hat, var_beta_med_hat, nX, nY, nZ, k_selected, l, m, pleiotropy, method)
    
    c(par, list(true_total = sim$true_total,
                true_alpha = sim$true_alpha,
                true_gamma_nonzero = sim$true_gamma_nonzero,
                true_delta_nonzero = sim$true_delta_nonzero,
                generate_seed = sim$generate_seed,
                estimate_seed = sim$estimate_seed))
}

parser <- ArgumentParser()
parser$add_argument("--method")
parser$add_argument("--ML_methods_R")
parser$add_argument("--in_simulation_RDS")
parser$add_argument("--out_simulation_RDS")
args <- parser$parse_args()

source(args$ML_methods_R)

sim_bundle <- readRDS(args$in_simulation_RDS)

sims <- list()
for (i in 1:length(sim_bundle)) {
    if (args$method == "naive") {
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "naive")
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "naive_zero")
    }
    if (args$method == "Burgess") {
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "Burgess")
    }
    else if (args$method == "original") {
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "original", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "original", pleiotropy = TRUE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "original_diagonal", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "original_diagonal", pleiotropy = TRUE)
    }
    else if (args$method == "integrated") {
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_free", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_free", pleiotropy = TRUE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_fixed_contained", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_fixed_contained", pleiotropy = TRUE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_fixed_naive", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_fixed_naive", pleiotropy = TRUE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_just_alpha", pleiotropy = FALSE)
        sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_just_alpha", pleiotropy = TRUE)
        #sims[[length(sims) + 1]] <- simulate(sim_bundle[[i]], method = "integrated_fixed_PCA")
    }
    
    cat("Simulation", i, "done:", format(Sys.time(), usetz = TRUE), "\n")
}

saveRDS(sims, args$out_simulation_RDS)
