suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--in_harmonized_gz")
parser$add_argument("--reference_dir")
parser$add_argument("--instrument_selection_tools_R")
parser$add_argument("--instrument_pval_threshold", type = "double")
parser$add_argument("--out_clumped_txt")
args <- parser$parse_args()

source(args$instrument_selection_tools_R)

dat_clumped <- fread(args$in_harmonized_gz) %>%
    filter(pval <= args$instrument_pval_threshold) %>%
    rename(pval.exposure = pval) %>%
    plink_clump(args$reference_dir) %>%
    rename(pval = pval.exposure) %>%
    relocate(SNP, .after = pos)

write_delim(dat_clumped, args$out_clumped_txt, delim = "\t")
