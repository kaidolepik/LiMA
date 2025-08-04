suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))
suppressPackageStartupMessages(library(argparse))

prune_mediators <- function(dat_mediators, mr_exposure_to_mediators, mediator_cor_threshold) {
    if (nrow(dat_mediators) == 0)
        return(character())
    
    mediator_cor <- dat_mediators %>%
        pivot_wider(id_cols = SNP, names_from = id, values_from = beta) %>%
        column_to_rownames("SNP") %>%
        cor(use = "pairwise.complete.obs")
    
    mediators <- mr_exposure_to_mediators
    
    i = 1
    while (i <= nrow(mediators)) {
        mediator <- mediators$mediator_ID[i]
        pruned_mediators <- colnames(mediator_cor)[abs(mediator_cor[mediator, ]) >= mediator_cor_threshold & colnames(mediator_cor) != mediator]
        mediators <- filter(mediators, mediator_ID == mediator | !(mediator_ID %in% pruned_mediators))
        
        i = i + 1
    }
    
    return(mediators$mediator_ID)
}

prepare_mediator_data <- function(dat_mediators, dat_exposure, dat_outcome, mr_exposure_to_mediators, reference_dir, instrument_pval_threshold, mediator_cor_threshold, clump) {
    # Get instrument strength as R2, make sure instrument is measured in the exposure, and filter by the Steiger test on the exposure
    mediator_instruments <- potential_instruments(dat_mediators, dat_exposure, instrument_pval_threshold) %>%
        select(SNP, chr, pos, id = id.exposure, rsq = rsq.exposure)
    
    # Make sure instrument is measured in the outcome and filter by the Steiger test on the outcome
    mediator_instruments <- potential_instruments(dat_mediators, dat_outcome, instrument_pval_threshold) %>%
        select(SNP, chr, pos, id = id.exposure) %>%
        merge(mediator_instruments, by = c("SNP", "chr", "pos", "id"))
    
    # Make sure instruments are measured in all mediators that have at least a single instrument of their own
    mediator_instruments <- dat_mediators %>%
        merge(distinct(mediator_instruments, id), by = "id") %>% # leave only those mediators that have instruments
        merge(distinct(mediator_instruments, SNP, chr, pos), by = c("SNP", "chr", "pos")) %>%
        group_by(SNP, chr, pos) %>%
        summarize(N = n(), .groups = "drop") %>%
        filter(N == n_distinct(mediator_instruments$id)) %>%
        select(-N) %>%
        merge(mediator_instruments)
    
    # Update the mediator data by filtering out mediators with no instruments
    dat_mediators <- filter(dat_mediators, id %in% unique(mediator_instruments$id))
    mr_exposure_to_mediators <- filter(mr_exposure_to_mediators, mediator_ID %in% unique(mediator_instruments$id))
    
    # Prune correlated mediators, keep those that have stronger causal effects from the exposure
    mediators <- prune_mediators(dat_mediators, mr_exposure_to_mediators, mediator_cor_threshold) # Mediators with at least one significant pre-clumping instrument (could consider post-clumping)
    
    # Update mediator data by filtering out pruned mediators
    dat_mediators <- filter(dat_mediators, id %in% mediators)
    mediator_instruments <- filter(mediator_instruments, id %in% mediators)
    
    # Clump the mediator instruments
    mediator_instruments_clumped <- select(mediator_instruments, -id)
    if (clump && nrow(mediator_instruments) > 0) {
        mediator_instruments_clumped <- mediator_instruments_clumped %>%
            group_by(SNP, chr, pos) %>%
            summarize(rsq = max(rsq)) %>%
            arrange(-rsq) %>%
            mutate(pval.exposure = 1:n() / (n() + 1)) %>%
            plink_clump(reference_dir) %>%
            select(SNP, chr, pos, rsq)
    }
    
    # Add respective mediator IDs for each mediator instrument
    mediator_instrument_IDs <- inner_join(mediator_instruments, mediator_instruments_clumped, 
                                          by = c("SNP" = "SNP", "chr" = "chr", "pos" = "pos")) %>%
        select(SNP, chr, pos, id)
    
    list(dat_mediators = dat_mediators, mediator_instruments = mediator_instruments_clumped, mediator_instrument_IDs = mediator_instrument_IDs)
}

prepare_exposure_instruments <- function(dat_mediators, dat_exposure, dat_outcome, reference_dir, instrument_pval_threshold, clump) {
    n_mediators <- n_distinct(dat_mediators$id)
    
    # Get instrument strength as R2, make sure instrument is measured in the outcome, and filter by the Steiger test on the outcome
    exposure_instruments <- potential_instruments(dat_exposure, dat_outcome, instrument_pval_threshold) %>%
        select(SNP, chr, pos, rsq = rsq.exposure)
    
    # Make sure instrument is measured in the mediators and filter by the Steiger test on the mediators
    exposure_instruments <- potential_instruments(dat_exposure, dat_mediators, instrument_pval_threshold) %>%
        select(SNP, chr, pos) %>%
        merge(exposure_instruments, by = c("SNP", "chr", "pos")) %>%
        group_by(SNP, chr, pos) %>%
        filter(n() == n_mediators) %>%
        summarize(rsq = mean(rsq), .groups = "drop") # mean(rsq) is the exposure instrument rsq, same for any single instrument
    
    if (clump && nrow(exposure_instruments) > 0) {
        exposure_instruments <- exposure_instruments %>%
            arrange(-rsq) %>%
            mutate(pval.exposure = 1:n() / (n() + 1)) %>%
            plink_clump(reference_dir) %>%
            select(SNP, chr, pos, rsq)
    }
    
    return(exposure_instruments)
}

clump_instruments <- function(mediator_instruments, exposure_instruments, reference_dir) {
    if (nrow(mediator_instruments) != 0 && nrow(exposure_instruments) != 0) {
        # Concatenate exposure and mediator instruments together for clumping
        instruments <- mediator_instruments %>%
            group_by(SNP, chr, pos) %>%
            summarize(rsq = 1 + max(rsq), .groups = "drop") %>% # remove duplicate mediator instruments and prioritize mediator instruments over exposure instruments
            rbind(exposure_instruments) %>% # concatenate exposure instruments
            group_by(SNP, chr, pos) %>%
            filter(n() == 1) %>% # exposure and mediators can't have the same instrument
            summarize(rsq = max(rsq), .groups = "drop") %>%
            arrange(desc(rsq)) %>%
            mutate(pval.exposure = 1:n() / (n() + 1)) # Add a dummy P-value for clumping
        
        if (nrow(instruments) > 0)
            instruments <- plink_clump(instruments, reference_dir) %>%
                select(SNP, chr, pos)
        
        exposure_instruments <- exposure_instruments %>%
            merge(instruments, by = c("SNP", "chr", "pos")) %>%
            select(SNP, chr, pos)
        
        mediator_instruments <- mediator_instruments %>%
            distinct(SNP, chr, pos) %>%
            merge(instruments, by = c("SNP", "chr", "pos")) %>%
            select(SNP, chr, pos)
    }
    
    list(exposure_instruments = exposure_instruments, mediator_instruments = mediator_instruments)
}

prepare_summary_statistics <- function(dat_mediators, dat_exposure, dat_outcome, exposure_instruments, mediator_instruments) {
    dat_exposure <- data.frame(dat_exposure, row.names = dat_exposure$SNP, stringsAsFactors = FALSE)
    dat_outcome <- data.frame(dat_outcome, row.names = dat_outcome$SNP, stringsAsFactors = FALSE)
    mediators <- unique(dat_mediators$id)
    k <- length(mediators)
    l <- nrow(mediator_instruments)
    m <- nrow(exposure_instruments)
    
    if (l == 0 || m == 0)
        return(list(status = 1, k=k, l=l, m=m))
        
    med_hat <- dat_mediators %>%
        pivot_wider(id_cols = c(SNP, chr, pos), names_from = id, values_from = beta) %>%
        column_to_rownames("SNP")
    Sigma <- estimate_Sigma(med_hat[mediators]) # Estimated on all the available SNPs, Ideally, it would need to be estimated on a bigger set of SNPs
    
    beta_hat <- pull(dat_exposure[exposure_instruments$SNP, ], beta, SNP)
    beta_med_hat <- pull(dat_exposure[mediator_instruments$SNP, ], beta, SNP)
    B_hat <- as.matrix(med_hat[mediator_instruments$SNP, mediators, drop = FALSE])
    C_hat <- as.matrix(med_hat[exposure_instruments$SNP, mediators, drop = FALSE])
    b_hat <- pull(dat_outcome[mediator_instruments$SNP, ], beta, SNP)
    c_hat <- pull(dat_outcome[exposure_instruments$SNP, ], beta, SNP)
    
    var_mediators <- dat_mediators %>%
        pivot_wider(id_cols = c(SNP, chr, pos), names_from = id, values_from = se) %>%
        column_to_rownames("SNP")
    var_B_hat <- as.matrix(var_mediators[mediator_instruments$SNP, mediators, drop = FALSE])^2
    var_C_hat <- as.matrix(var_mediators[exposure_instruments$SNP, mediators, drop = FALSE])^2
    var_b_hat <- pull(dat_outcome[mediator_instruments$SNP, ], se, SNP)^2
    var_c_hat <- pull(dat_outcome[exposure_instruments$SNP, ], se, SNP)^2
    var_beta_hat <- pull(dat_exposure[exposure_instruments$SNP, ], se, SNP)^2
    var_beta_med_hat <- pull(dat_exposure[mediator_instruments$SNP, ], se, SNP)^2
    
    # Danger: if a matrix didn't have column names then colnames(matrix) would return NULL, and all(NULL == c("name1", "name2")) == TRUE
    if (!(all(names(beta_hat) == names(c_hat)) && all(names(beta_hat) == rownames(C_hat)) && all(names(beta_med_hat) == names(b_hat)) && all(names(beta_med_hat) == rownames(B_hat)) && all(colnames(C_hat) == colnames(B_hat)) &&
          all(rownames(var_C_hat) == rownames(C_hat)) && all(colnames(var_C_hat) == colnames(C_hat)) && all(names(var_c_hat) == names(c_hat)) && all(names(var_b_hat) == names(b_hat))))
        return(list(status = 2, k=k, l=l, m=m))
    
    instruments <- rbind(exposure_instruments, mediator_instruments) %>%
        distinct()
    
    nX <- pull(dat_exposure[instruments$SNP, ], samplesize) %>%
        mean() %>%
        round()
    nY <- pull(dat_outcome[instruments$SNP, ], samplesize) %>%
        mean() %>%
        round()
    nZ <- dat_mediators %>%
        merge(instruments) %>%
        pull(samplesize) %>%
        mean() %>%
        round()
    
    list(status = 0, k=k, l=l, m=m, nX=nX, nY=nY, nZ=nZ, B_hat=B_hat, C_hat=C_hat, b_hat=b_hat, c_hat=c_hat, beta_hat=beta_hat, beta_med_hat=beta_med_hat, 
         Sigma=Sigma, var_B_hat=var_B_hat, var_C_hat=var_C_hat, var_b_hat=var_b_hat, var_c_hat=var_c_hat, var_beta_hat=var_beta_hat, var_beta_med_hat=var_beta_med_hat)
}

parser <- ArgumentParser()
parser$add_argument("--exposure_ID")
parser$add_argument("--outcome_ID")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--mediator_cor_threshold", type = "double")
parser$add_argument("--clump")
parser$add_argument("--reference_dir")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--ML_methods_R")
parser$add_argument("--in_mediation_data_RDS")
parser$add_argument("--in_mr_results_txt")
parser$add_argument("--out_summary_statistics_RDS")
args <- parser$parse_args()
args$clump <- as.logical(args$clump)

source(args$instrument_selection_tools_R)
source(args$ML_methods_R)

# Read in harmonized data about the exposure, mediators and outcome
mediation_data <- readRDS(args$in_mediation_data_RDS)
dat_exposure <- mediation_data$dat_exposure
dat_outcome <- mediation_data$dat_outcome
dat_mediators <- mediation_data$dat_mediators

# Read in MR IVW results from the exposure to the mediators
mr_exposure_to_mediators <- fread(args$in_mr_results_txt) %>%
    filter(exposure == args$exposure_ID,
           outcome %in% unique(dat_mediators$id)) %>%
    arrange(pval) %>%
    select(mediator_ID = outcome, pval)

# Select mediator instruments and keep only instrumented mediators in the mediator data
mediators <- prepare_mediator_data(dat_mediators, dat_exposure, dat_outcome, mr_exposure_to_mediators, args$reference_dir, args$instrument_pval_threshold, args$mediator_cor_threshold, args$clump)
dat_mediators <- mediators$dat_mediators                # mediators have instruments and are pruned by cor_threshold
mediator_instruments <- mediators$mediator_instruments  # instruments are measured in all mediators, pass Steiger filters and are clumped

# Select exposure instruments: measured in all traits, pass Steiger filters, clumped
exposure_instruments <- prepare_exposure_instruments(dat_mediators, dat_exposure, dat_outcome, args$reference_dir, args$instrument_pval_threshold, args$clump)

# Clump exposure and mediator instruments together
instruments <- clump_instruments(mediator_instruments, exposure_instruments, args$reference_dir)
mediator_instruments <- instruments$mediator_instruments
exposure_instruments <- instruments$exposure_instruments

# Prepare summary statistics
summary_statistics <- prepare_summary_statistics(dat_mediators, dat_exposure, dat_outcome, exposure_instruments, mediator_instruments)
summary_statistics[["mediator_instrument_IDs"]] <- mediator_instruments %>%
    inner_join(mediators$mediator_instrument_IDs, by = c("SNP" = "SNP", "chr" = "chr", "pos" = "pos")) %>% 
    select(SNP, id) # Make sure we know which mediators have which instruments
    
# Save the data
saveRDS(summary_statistics, args$out_summary_statistics_RDS)
