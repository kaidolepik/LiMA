suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

# Match SNPs with common_SNPs (by Name) and return with the minimum P-value
match_SNPs <- function(common_SNPs, SNPs) {
    common_SNPs <- common_SNPs %>%
        mutate(Pval.y = SNPs$Pval[match(Name, SNPs$Name)]) %>%
        drop_na(Pval.y) %>%
        mutate(Pval = pmin(Pval, Pval.y)) %>%
        select(-Pval.y)
    
    return(common_SNPs)
}

# Return common SNPs from harmonized gz-files (e.g. harmonized UKBB trait files)
common_trait_SNPs <- function(files_dir) {
    common_SNPs <- data.table()
    
    files <- list.files(files_dir, full.names = TRUE)
    for (i in 1:length(files)) {
        SNPs <- fread(files[i], select = c("chr", "pos", "SNP", "effect_allele", "other_allele", "pval"), col.names = c("chr", "pos", "SNP", "A1", "A2", "Pval")) %>%
            mutate(Name = paste0("chr", chr, ":", pos, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A2 > A1, A2, A1)))
        
        if (i == 1)
            common_SNPs <- SNPs
        else
            common_SNPs <- match_SNPs(common_SNPs, SNPs)
    }

    common_SNPs <- select(common_SNPs, Name, SNP)

    return(common_SNPs)
}

# Return common SNPs from mediator gz-files (requires Name to match the format below in the code)
# We don't make sure that alleles are ordered lexicographically because
# 1) we can expect all mediator files to follow the same standard
# 2) we need to create a besd-fileset so we can't change the names
common_mediator_SNPs_gz <- function(files_dir) {
    common_SNPs <- data.table()
    
    files <- list.files(files_dir, full.names = TRUE)
    for (i in 1:length(files)) {
        SNPs <- fread(files[i], select = c("Name", "Pval")) %>%
            filter(grepl("^chr[0-9][0-9]?:[0-9]+:[ACGT]:[ACGT]$", Name)) %>%
            mutate(chr_pos = substr(Name, 1, nchar(Name) - 4)) %>%
            filter(duplicated(chr_pos, fromLast = FALSE) + duplicated(chr_pos, fromLast = TRUE) == 0) %>%
            select(Name, Pval)
        
        if (i == 1)
            common_SNPs <- SNPs
        else
            common_SNPs <- match_SNPs(common_SNPs, SNPs)
    }
    
    return(common_SNPs)
}

# Return common SNPs from the bed-format mediator files and save with the minimum P-value for each SNP
common_mediator_SNPs_bed <- function(files_dir) {
    common_SNPs <- data.table()
    
    files <- list.files(files_dir, full.names = TRUE)
    for (i in 1:length(files)) {
        SNPs <- fread(files[i], header = FALSE, col.names = c("Name", "Pval")) # colClasses = list(integer = c(2, 3)) created numeric classes
        
        if (i == 1)
            common_SNPs <- SNPs
        else
            common_SNPs <- match_SNPs(common_SNPs, SNPs)
    }
    
    common_SNPs <- common_SNPs %>%
        mutate(chr_pos = substr(Name, 1, nchar(Name) - 4),
               A1_A2 = substr(Name, nchar(Name) - 2, nchar(Name))) %>%
        filter(duplicated(chr_pos, fromLast = FALSE) + duplicated(chr_pos, fromLast = TRUE) == 0) %>% # Leave only unique chr-pos combinations
        mutate(chr = sub(":.*", "", chr_pos),
               pos = substr(chr_pos, nchar(chr) + 2, nchar(chr_pos)),
               A1 = substr(A1_A2, 1, 1),
               A2 = substr(A1_A2, 3, 3)) %>%
        select(chr, pos, Pval, A1, A2) %>%
        arrange(as.integer(sub("^chr", "", chr)), pos) %>%
        mutate(chr_assembly = chr, # Add additional columns reflecting genome assembly
               pos_assembly = pos, # This standardizes the file regardless of whether liftOver follows or not
               ID = 1:n(),
               .after = Pval)
    
    return(common_SNPs)
}

# Return common mediator and trait SNPs merged together
common_SNPs_merged <- function(files_dir) {
    # Common mediator SNPs that legally passed liftOver
    common_mediator_SNPs <- paste0(files_dir, "/hg19.bed") %>%
        fread(header = FALSE, select = c(1, 2, 7, 8, 3, 4, 5), col.names = c("chr", "pos", "A1", "A2", "Pval", "chr_assembly", "pos_assembly")) %>%
        filter(chr == chr_assembly) %>%
        mutate(chr = as.integer(sub("^chr", "", chr)))
    
    # Common trait SNPs from harmonized files
    common_trait_SNPs <- paste0(files_dir, "/common_trait_SNPs.txt") %>%
        fread(header = FALSE, select = 1, col.names = "Name")
    
    # Merge common mediator and trait SNPs
    common_SNPs <- common_mediator_SNPs %>%
        mutate(chr_pos = paste0("chr", chr, ":", pos),
               Name = paste0(chr_pos, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A2 > A1, A2, A1))) %>%
        merge(common_trait_SNPs, by = "Name") %>%
        filter(duplicated(chr_pos, fromLast = FALSE) + duplicated(chr_pos, fromLast = TRUE) == 0) %>%
        select(chr, pos, Pval, chr_assembly, pos_assembly) # There's no need to include the Name anymore because we know we're dealing with biallelic unique SNPs
    
    return(common_SNPs)
}

# Return significant mediator SNPs + clumped trait SNPs so that all necessary SNPs would be available in the final besd-fileset
common_significant_SNPs <- function(files_dir, reference_dir, instrument_pval_threshold) {
    common_SNPs <- fread(paste0(files_dir, "/common_SNPs_merged.txt"), header = FALSE, col.names = c("chr", "pos", "Pval", "chr_assembly", "pos_assembly"))

    # Common mediator SNPs passing the P-value threshold    
    significant_SNPs <- common_SNPs %>%
        filter(Pval <= instrument_pval_threshold) %>%
        select(chr, pos, chr_assembly, pos_assembly)
    
    common_SNPs <- common_SNPs %>%
        select(chr, pos, chr_assembly, pos_assembly)
    
    # Add in significant clumped SNPs of traits
    harmonized_trait_files <- list.files(files_dir, pattern = ".gz", full.names = TRUE)
    for (i in 1:length(harmonized_trait_files)) {
        significant_SNPs <- fread(harmonized_trait_files[i], header = TRUE, select = c("chr", "pos", "SNP", "pval"), col.names = c("chr", "pos", "SNP", "pval.exposure")) %>%
            filter(pval.exposure <= instrument_pval_threshold) %>%
            merge(common_SNPs, by = c("chr", "pos")) %>%
            plink_clump(reference_dir) %>%
            select(chr, pos, chr_assembly, pos_assembly) %>%
            as.data.table() %>%
            funion(significant_SNPs)
    }

    significant_SNPs <- significant_SNPs %>%
        arrange(chr, pos)
    
    return(significant_SNPs)
}

# Return common SNPs from harmonized gz-files (e.g. harmonized UKBB trait files)
common_analysis_SNPs <- function(files_dir, prefix) {
    # Common trait SNPs from harmonized files
    common_trait_SNPs <- paste0(files_dir, "/common_trait_SNPs.txt") %>%
        fread(header = FALSE, col.names = c("Name", "ID_trait"))
    
    # Common mediator SNPs from the esi-file
    common_esi_SNPs <- paste0(files_dir, "/", prefix, ".esi") %>%
        fread(select = c(2, 1, 4, 5, 6), col.names = c("ID_esi", "chr", "pos", "A1", "A2")) %>%
        mutate(chr_pos = paste0("chr", chr, ":", pos),
               Name = paste0(chr_pos, ":", ifelse(A1 < A2, A1, A2), ":", ifelse(A2 > A1, A2, A1))) %>%
        filter(duplicated(chr_pos, fromLast = FALSE) + duplicated(chr_pos, fromLast = TRUE) == 0)

    common_SNPs <- common_esi_SNPs %>%
        merge(common_trait_SNPs, by = c("Name")) %>%
        select(Name, ID_esi, ID_trait)
    
    return(common_SNPs)
}

parser <- ArgumentParser()
parser$add_argument("--in_files_dir")
parser$add_argument("--type")
parser$add_argument("--reference_dir")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--out_SNPs")
args <- parser$parse_args()

if (args$type == "harmonized_trait") {
    common_SNPs <- common_trait_SNPs(args$in_files_dir)
} else if (args$type == "mediator_gz") {
    common_SNPs <- common_mediator_SNPs_gz(args$in_files_dir)
} else if (args$type == "mediator_bed") {
    common_SNPs <- common_mediator_SNPs_bed(args$in_files_dir)
} else if (args$type == "merge") {
    common_SNPs <- common_SNPs_merged(args$in_files_dir)
} else if (args$type == "significant") {
    source(args$instrument_selection_tools_R)
    common_SNPs <- common_significant_SNPs(args$in_files_dir, args$reference_dir, args$instrument_pval_threshold)
} else if (args$type %in% c("INTERVAL_plasma_proteins")) {
    common_SNPs <- common_analysis_SNPs(args$in_files_dir, args$type)
}

write_delim(common_SNPs, args$out_SNPs, col_names = FALSE)
