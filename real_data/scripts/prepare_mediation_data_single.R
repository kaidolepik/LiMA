suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--exposure_ID")
parser$add_argument("--outcome_ID")
parser$add_argument("--in_exposure_gz")
parser$add_argument("--in_outcome_gz")
parser$add_argument("--in_mediators_as_exposure")
parser$add_argument("--in_mediators_as_outcome")
parser$add_argument("--in_mr_ivw_results_exposure")
parser$add_argument("--in_mr_ivw_results_outcome")
parser$add_argument("--mediator_pval_threshold_from_exposure")
parser$add_argument("--mediator_pval_threshold_to_outcome")
parser$add_argument("--out_mediation_data_RDS")
args <- parser$parse_args()

# Mediators that have a significant causal effect from the exposure
significant_mediators_exposure <- fread(args$in_mr_ivw_results_exposure) %>%
    drop_na(pval) %>%
    filter(exposure == args$exposure_ID)

pval_threshold_exposure <- ifelse(args$mediator_pval_threshold_from_exposure == "Bonferroni", 
    0.05 / nrow(significant_mediators_exposure), 
    as.numeric(args$mediator_pval_threshold_from_exposure))

significant_mediators_exposure <- filter(significant_mediators_exposure, pval <= pval_threshold_exposure) %>%
pull(outcome)

# Mediators that have a significant causal effect to the outcome
significant_mediators_outcome <- fread(args$in_mr_ivw_results_outcome) %>%
    drop_na(pval) %>%
    filter(outcome == args$outcome_ID)

pval_threshold_outcome <- ifelse(args$mediator_pval_threshold_to_outcome == "Bonferroni", 
    0.05 / nrow(significant_mediators_outcome), 
    as.numeric(args$mediator_pval_threshold_to_outcome))

significant_mediators_outcome <- filter(significant_mediators_outcome, pval <= pval_threshold_outcome) %>%
    pull(exposure)

significant_mediators <- intersect(significant_mediators_exposure, significant_mediators_outcome)

# Combine data together for mediators with significant causal effect from the exposure
dat_mediators_as_exposure <- fread(args$in_mediators_as_exposure, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos"))) %>%
    filter(id %in% significant_mediators)
dat_mediators_as_outcome <- fread(args$in_mediators_as_outcome, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos"))) %>%
    filter(id %in% significant_mediators)
dat_mediators <- rbind(dat_mediators_as_exposure, dat_mediators_as_outcome) %>%
    distinct()

# Read in harmonized exposure and outcome data of the mediator SNPs
dat_exposure <- fread(args$in_exposure_gz, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos"))) %>%
    merge(distinct(dat_mediators[, c("chr", "pos")]), by = c("chr", "pos"))
dat_outcome <- fread(args$in_outcome_gz, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos"))) %>%
    merge(distinct(dat_mediators[, c("chr", "pos")]), by = c("chr", "pos"))

mediation_data <- list(dat_exposure = dat_exposure, dat_outcome = dat_outcome, dat_mediators = dat_mediators)
saveRDS(mediation_data, args$out_mediation_data_RDS)
