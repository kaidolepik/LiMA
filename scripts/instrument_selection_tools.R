suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(TwoSampleMR))

std_beta <- function(Z, N) {
    Z / sqrt(N + Z^2)
}

std_se <- function(Z, N) {
    1 / sqrt(N + Z^2)
}

is_biallelic <- function(chr, pos) {
    dat <- tibble(chr, pos)
    
    !(duplicated(dat, fromLast = FALSE) | duplicated(dat, fromLast = TRUE))
}

is_nonambiguous_SNP <- function(A1, A2) {
    (A1 %in% c("G", "C") & A2 %in% c("A", "T")) | (A1 %in% c("A", "T") & A2 %in% c("G", "C"))
}

alleles_match <- function(A1, A2, A1_bim, A2_bim) {
    (A1 == A1_bim & A2 == A2_bim) | (A1 == A2_bim & A2 == A1_bim)
}

harmonize_chr_data <- function(bim_filename, dat_trait) {
    dat_harmonized <- fread(bim_filename, col.names = c("chr", "SNP", "dummy", "pos", "A1_bim", "A2_bim")) %>%
        merge(dat_trait) %>%
        group_by(id) %>% # In case there are multiple traits in dat_trait
        filter(is_biallelic(chr, pos),
               is_nonambiguous_SNP(A1_bim, A2_bim),
               alleles_match(effect_allele, other_allele, A1_bim, A2_bim)) %>%
        mutate(flip = (effect_allele != A1_bim),
               beta = ifelse(flip, -beta, beta), 
               z = ifelse(flip, -z, z), 
               eaf = ifelse(flip, 1 - eaf, eaf), 
               effect_allele = ifelse(flip, A1_bim, effect_allele), 
               other_allele = ifelse(flip, A2_bim, other_allele)) %>%
        select(-c(dummy, A1_bim, A2_bim, flip)) %>%
        ungroup()
    
    return(dat_harmonized)
}

harmonize_data <- function(dat_trait, reference_dir) {
    reference_chr_bims <- list.files(path = reference_dir, pattern = ".*\\.bim", full.names = TRUE, recursive = TRUE)
    
    dat_harmonized_trait <- reference_chr_bims %>%
        lapply(harmonize_chr_data, dat_trait) %>%
        rbindlist() %>%
        unique()
    
    return(dat_harmonized_trait)
}

plink_clump_chr <- function(dat_mr, in_reference_chr_prefix, p1 = 1, p2 = 1, r2 = 0.001, kb = 10000) {
    cat(paste0("---Clumping: ", format(Sys.time(), "%X"), "\n"))
    
    tmp_dir <- basename(in_reference_chr_prefix)
    if (!dir.exists(tmp_dir))
        dir.create(tmp_dir)
    
    file_extract <- paste0(tmp_dir, "/extract.txt")
    write_delim(select(dat_mr, SNP), file_extract, col_names = FALSE)
    
    file_potential_instruments <- paste0(tmp_dir, "/potential_instruments.txt")
    write_delim(select(dat_mr, SNP, P = pval.exposure), file_potential_instruments)
    
    file_out <- paste0(tmp_dir, "/clump")
    command <- paste0("plink --bfile ", in_reference_chr_prefix, " ",
                      "--extract ", file_extract, " ",
                      "--clump ", file_potential_instruments, " ",
                      "--clump-p1 ", p1, " ",
                      "--clump-p2 ", p2, " ",
                      "--clump-r2 ", r2, " ",
                      "--clump-kb ", kb, " ",
                      "--out ", file_out)
    system(command, ignore.stdout = TRUE)
    
    dat_mr <- fread(paste0(file_out, ".clumped")) %>%
        select(SNP) %>%
        inner_join(dat_mr, by = c("SNP"))
    
    system(paste0("rm -r ", tmp_dir), ignore.stdout = TRUE)
    
    return(dat_mr)
}

plink_clump <- function(dat_mr, reference_dir) {
    reference_chr_prefixes <- list.files(path = reference_dir, pattern = ".*\\.bim", full.names = TRUE, recursive = TRUE) %>%
        str_replace("\\.bim$", "") %>%
        set_names(str_match(., ".*UK10K_chr(\\d+)$")[, 2])
    
    dat_mr <- dat_mr %>%
        group_by(chr) %>%
        group_modify(~ plink_clump_chr(., in_reference_chr_prefix = reference_chr_prefixes[as.character(.y)])) %>%
        ungroup()
    
    return(dat_mr)
}

potential_instruments <- function(dat_exposure, dat_outcome, instrument_pval_threshold, steiger = TRUE, steiger_pval_threshold = 0.05) {
    dat_mr <- dat_exposure %>%
        filter(pval <= instrument_pval_threshold) %>%
        merge(dat_outcome, by = c("SNP", "chr", "pos", "effect_allele", "other_allele"), suffixes = c(".exposure", ".outcome")) %>%
        mutate(exposure = "exposure", outcome = "outcome")
    
    if (nrow(dat_mr) == 0) 
        dat_mr <- mutate(dat_mr, units.outcome = NA, units.exposure = NA, rsq.exposure = NA, effective_n.exposure = NA, 
                         rsq.outcome = NA, effective_n.outcome = NA, steiger_dir = NA, steiger_pval = NA)
    else if (steiger)
        dat_mr <- TwoSampleMR::steiger_filtering(dat_mr) %>%
            filter(!(!steiger_dir & steiger_pval <= steiger_pval_threshold))
    
    return(dat_mr)
}

prepare_mr_data <- function(dat_exposure, dat_outcome, reference_dir, instrument_pval_threshold, clump = TRUE) {
    dat_mr <- potential_instruments(dat_exposure, dat_outcome, instrument_pval_threshold)
    
    if (clump && nrow(dat_mr) > 0)
        dat_mr <- plink_clump(dat_mr, reference_dir)
    
    dat_mr <- mutate(dat_mr, mr_keep = TRUE)
    
    return(dat_mr)
}

mr_ivw_analysis <- function(dat_mr, exposure_ID, outcome_ID) {
    if (nrow(dat_mr) == 0)
        mr_ivw_results <- tibble(exposure = exposure_ID, outcome = outcome_ID, nsnp = 0, b = NA, se = NA, pval = NA)
    else {
        mr_ivw_results <- mr(dat_mr, method_list = c("mr_wald_ratio", "mr_ivw")) %>%
            select(exposure = id.exposure, outcome = id.outcome, nsnp, b, se, pval)
    }
    
    return(mr_ivw_results)
}
