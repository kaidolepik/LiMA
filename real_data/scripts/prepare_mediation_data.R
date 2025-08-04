suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(argparse))

significant_mediator_SNPs <- function(mediators_dir, instrument_pval_threshold) {
    dat_SNPs <- NULL
    
    for (mediator_gz in list.files(mediators_dir, full.names = TRUE)) {
        dat_mediator <- fread(mediator_gz, select = c("SNP", "chr", "pos", "pval")) %>%
            filter(pval <= instrument_pval_threshold) %>%
            select(SNP, chr, pos)
        
        if (is.null(dat_SNPs))
            dat_SNPs <- dat_mediator
        else
            dat_SNPs <- merge(dat_SNPs, dat_mediator, by = c("SNP", "chr", "pos"), all = TRUE)
    }
    
    return(dat_SNPs)
}

mediator_data <- function(mediators_dir, significant_SNPs) {
    dat_mediators <- list.files(mediators_dir, full.names = TRUE) %>%
        lapply(function(mediator_gz) {
            fread(mediator_gz) %>%
                merge(significant_SNPs, by = c("SNP", "chr", "pos"))
        }) %>%
        rbindlist()
    
    return(dat_mediators)
}

parser <- ArgumentParser()
parser$add_argument("--in_exposure_gz")
parser$add_argument("--in_outcome_gz")
parser$add_argument("--mediator_dir")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--out_mediation_data_RDS")
args <- parser$parse_args()

# Read in the exposure harmonized data and significant exposure SNPs
dat_exposure <- fread(args$in_exposure_gz)
exposure_SNPs <- dat_exposure %>%
    filter(pval <= args$instrument_pval_threshold) %>%
    select(SNP, chr, pos)

# Read in significant mediator SNPs and merge with exposure SNPs
significant_SNPs <- significant_mediator_SNPs(args$mediator_dir, args$instrument_pval_threshold) %>%
    merge(exposure_SNPs, by = c("SNP", "chr", "pos"), all = TRUE)
stopifnot(nrow(significant_SNPs) == nrow(distinct(significant_SNPs)))

# Leave only significant SNPs of exposure, outcome and mediators data
dat_exposure <- merge(dat_exposure, significant_SNPs, by = c("SNP", "chr", "pos"))
dat_outcome <- merge(fread(args$in_outcome_gz), significant_SNPs, by = c("SNP", "chr", "pos"))
dat_mediators <- mediator_data(args$mediator_dir, significant_SNPs)

# Save the data for memory efficiency
mediation_data <- list(dat_exposure = dat_exposure, dat_outcome = dat_outcome, dat_mediators = dat_mediators)
saveRDS(mediation_data, args$out_mediation_data_RDS)
