suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))

analyze <- function(x, method, pleiotropy) {
    # Initialize the parameter vector of results
    par <- list(method=method, pleiotropy=pleiotropy, k_selected=x$k, l=x$l, m=x$m)
    
    # If the status was not zero then the summary statistics were not generated
    if (summary_statistics$status == 0) {
        par <- mediation_analysis(x$B_hat, x$C_hat, x$b_hat, x$c_hat, x$beta_hat, x$beta_med_hat, x$Sigma, 
                                  x$var_B_hat, x$var_C_hat, x$var_b_hat, x$var_c_hat, x$var_beta_hat, x$var_beta_med_hat,
                                  x$nX, x$nY, x$nZ, x$k, x$l, x$m, pleiotropy, method)
    }
    
    c(list(input_status = summary_statistics$status), par)
}

parser <- ArgumentParser()
parser$add_argument("--exposure_ID")
parser$add_argument("--exposure_type")
parser$add_argument("--outcome_ID")
parser$add_argument("--outcome_type")
parser$add_argument("--mediator_type")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--mediator_cor_threshold", type = "double")
parser$add_argument("--mediation_method")
parser$add_argument("--pleiotropy")
parser$add_argument("--ML_methods_R")
parser$add_argument("--in_summary_statistics_RDS")
parser$add_argument("--out_mediation_results_txt")
args <- parser$parse_args()
args$pleiotropy <- as.logical(args$pleiotropy)

source(args$ML_methods_R)
summary_statistics <- readRDS(args$in_summary_statistics_RDS)

# Perform the mediation analysis pipeline and output results
metadata <- list(exposure_ID = args$exposure_ID, exposure_type = args$exposure_type,
                 outcome_ID = args$outcome_ID, outcome_type = args$outcome_type, mediator_type = args$mediator_type,
                 instrument_pval_threshold = args$instrument_pval_threshold, mediator_cor_threshold = args$mediator_cor_threshold)

mediation_results <- c(metadata, analyze(summary_statistics, args$mediation_method, args$pleiotropy)) %>%
    unlist() %>%
    as_tibble_row()
colnames(mediation_results)[grepl("gamma", colnames(mediation_results))] <- paste0("gamma_", colnames(summary_statistics$B_hat))
colnames(mediation_results)[grepl("delta", colnames(mediation_results))] <- paste0("delta_", colnames(summary_statistics$B_hat))

if (mediation_results$input_status != 1) {
    mediation_results <- mediation_results %>%
        mutate_at(c("total", "alpha", "SE_total", "SE_alpha"), as.numeric) %>%
        mutate(MP = (total - alpha) / total,
               SE_MP = alpha^2/total^2 * (SE_alpha^2/alpha^2 + SE_total^2/total^2),
               P_MP = 2 * pnorm(-abs(MP / SE_MP)))
}

# Save the results
write_delim(mediation_results, args$out_mediation_results_txt, delim = "\t")
