suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--in_mediator_gz")
parser$add_argument("--in_common_SNPs")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--out_esd")
args <- parser$parse_args()

source(args$instrument_selection_tools_R)

common_SNPs <- fread(args$in_common_SNPs, header = FALSE, col.names = c("chr_hg19", "pos_hg19", "chr_assembly", "pos_assembly"))
MIN_NONZERO <- format(.Machine$double.xmin, digits = 6, trim = TRUE)

dat_esd <- fread(args$in_mediator_gz) %>%
    filter(grepl(paste0("^chr[0-9][0-9]?:[0-9]+:[ACGT]:[ACGT]$"), Name)) %>%
    mutate(N_Name = nchar(Name),
           chr = sub(":.*", "", Name),
           pos = as.integer(substr(Name, nchar(chr) + 2, N_Name - 4)),
           A1 = substr(Name, N_Name - 2, N_Name - 2),
           A2 = substr(Name, N_Name, N_Name),
           Freq = 0.5,
           z = Beta / SE,
           Beta = format(std_beta(z, N), digits = 6, trim = TRUE),
           SE = format(std_se(z, N), digits = 6, trim = TRUE)) %>%
    merge(common_SNPs, by.x = c("chr", "pos"), by.y = c("chr_assembly", "pos_assembly")) %>%
    mutate(Name = paste0(chr_hg19, ":", pos_hg19, ":", A1, ":", A2),
           Beta = ifelse(Pval == 1 | as.numeric(Beta) == 0, MIN_NONZERO, Beta),
           Pval = ifelse(Beta == MIN_NONZERO, 1, Pval),
           Pval = ifelse(Pval == 0, MIN_NONZERO, Pval)) %>%
    select(Chr = chr_hg19, SNP = Name, Bp = pos_hg19, A1, A2, Freq, Beta, se = SE, p = Pval)

fwrite(dat_esd, args$out_esd, sep = " ", col.names = TRUE, row.names = FALSE, quote = FALSE)
