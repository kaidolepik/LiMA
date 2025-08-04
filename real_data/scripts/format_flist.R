suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()
parser$add_argument("--in_esd_dir")
parser$add_argument("--in_annotation_file")
parser$add_argument("--out_flist")
args <- parser$parse_args()

probes <- list.files(args$in_esd_dir, pattern = "*.esd", full.names = FALSE) %>%
    as_tibble_col("PathOfEsd") %>%
    mutate(SeqId = sub(".esd", "", basename(PathOfEsd)), .before = 1)

dat_flist <- fread(args$in_annotation_file) %>%
    mutate(GeneticDistance = 0,
           Orientation = "+",
           TSS = ifelse(TSS != 0, TSS, 1)) %>%
    merge(probes, by = "SeqId") %>%
    select(Chr, ProbeID = SeqId, GeneticDistance, ProbeBp = TSS, Gene, Orientation, PathOfEsd)

write.table(dat_flist, args$out_flist, col.names = TRUE, row.names = FALSE, quote = FALSE)
