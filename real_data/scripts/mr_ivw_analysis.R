suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--exposure_ID")
parser$add_argument("--in_exposure_gz")
parser$add_argument("--outcome_ID")
parser$add_argument("--in_outcome_gz")
parser$add_argument("--reference_dir")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--out_mr_ivw_results_txt")
args <- parser$parse_args()

source(args$instrument_selection_tools_R)

dat_exposure <- fread(args$in_exposure_gz, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos")))
dat_outcome <- fread(args$in_outcome_gz, colClasses = list(character = c("SNP", "effect_allele", "other_allele"), integer = c("chr", "pos")))

dat_mr <- prepare_mr_data(dat_exposure, dat_outcome, args$reference_dir, args$instrument_pval_threshold, clump = !(args$outcome_ID %in% c("INTERVAL_plasma_proteins")))

mr_ivw_results <- mr_ivw_analysis(dat_mr, args$exposure_ID, args$outcome_ID)
write_delim(mr_ivw_results, args$out_mr_ivw_results_txt)
