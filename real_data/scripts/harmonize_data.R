suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(magrittr))

read_besd_fileset <- function(mediator_ID, phenotype, filename) {
    cat(paste0("-Read in besd-fileset ", phenotype, " from ", filename, ": ", format(Sys.time(), "%X"), "\n"))
    
    dat <- fread(filename, select = c("Chr", "BP", "A1", "A2", "Freq", "b", "SE", "p", "Probe"),
                 col.names = c("chr", "pos", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "id"),
                 colClasses = list(character = c("A1", "A2", "Probe"), integer = c("Chr", "BP"), numeric = c("Freq", "b", "SE", "p"))) %>%
        filter(id == mediator_ID | mediator_ID == phenotype) %>%
        mutate(z = beta / se,
               samplesize = round(1 / se**2),
               Phenotype = phenotype) %>%
        select(chr, pos, effect_allele, other_allele, eaf, beta, se, z, pval, samplesize, Phenotype, id)
    
    return(dat)
}

read_Lotta_et_al_2021 <- function(metabolite_ID, filename) {
    cat(paste0("-Read in Lotta et al. 2021 metabolite ", metabolite_ID, " from ", filename, ": ", format(Sys.time(), "%X"), "\n"))

    dat <- fread(filename, col.names = c("id", "chr", "pos", "effect_allele", "other_allele", "eaf", "samplesize", "z")) %>%
        mutate(effect_allele = toupper(effect_allele),
               other_allele = toupper(other_allele),
               beta = std_beta(z, samplesize),
               se = std_se(z, samplesize),
               pval = 2*(1 - pnorm(abs(z))),
               Phenotype = "Lotta_et_al_2021") %>%
        select(chr, pos, effect_allele, other_allele, eaf, beta, se, z, pval, samplesize, Phenotype, id)
    
    return(dat)
}

read_UKBB <- function(UKBB_ID, filename, type = c("quantitative", "binary")) {
    cat(paste0("-Read in UKBB trait ", UKBB_ID, " from ", filename, ": ", format(Sys.time(), "%X"), "\n"))
    
    drop_cols <- NULL
    if (match.arg(type) == "binary")
        drop_cols <- c(4)
    
    dat <- fread(cmd = paste0("gzip -cd ", filename), drop = drop_cols, col.names = c("variant", "minor_allele", "minor_AF", "low_confidence_variant", "samplesize", "AC", "ytx", "beta", "se", "z", "pval"))
    chr_pos_other_effect <- str_split_fixed(dat$variant, ":", n = 4) %>%
        set_colnames(c("chr", "pos", "other_allele", "effect_allele")) %>%
        as_tibble()
    dat <- bind_cols(chr_pos_other_effect, dat) %>%
        filter(!low_confidence_variant,
               chr %in% paste0(1:22)) %>%
        mutate(chr = as.integer(chr),
               pos = as.integer(pos),
               eaf = ifelse(effect_allele == minor_allele, minor_AF, 1 - minor_AF),
               beta = std_beta(z, samplesize),
               se = std_se(z, samplesize),
               Phenotype = paste0("UKBB_", type),
               id = UKBB_ID) %>%
        select(chr, pos, effect_allele, other_allele, eaf, beta, se, z, pval, samplesize, Phenotype, id)
    
    return(dat)
}

read_trait <- function(trait_ID, filename, type = c("UKBB_quantitative", "UKBB_binary", "Lotta_et_al_2021", "INTERVAL_plasma_proteins")) {
    type <- match.arg(type)
    
    if (type == "UKBB_quantitative")
        dat <- read_UKBB(trait_ID, filename, "quantitative")
    else if (type == "UKBB_binary")
        dat <- read_UKBB(trait_ID, filename, "binary")
    else if (type == "Lotta_et_al_2021")
        dat <- read_Lotta_et_al_2021(trait_ID, filename)
    else if (type %in% c("INTERVAL_plasma_proteins"))
        dat <- read_besd_fileset(trait_ID, type, filename)
    
    return(dat)
}

parser <- ArgumentParser()
parser$add_argument("--trait_ID")
parser$add_argument("--trait_type")
parser$add_argument("--in_trait_gz")
parser$add_argument("--reference_dir")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--out_trait_gz")
args <- parser$parse_args()

source(args$instrument_selection_tools_R)

dat_harmonized_trait <- read_trait(args$trait_ID, args$in_trait_gz, args$trait_type) %>%
    harmonize_data(args$reference_dir)

write_delim(dat_harmonized_trait, gzfile(args$out_trait_gz), delim = "\t")
