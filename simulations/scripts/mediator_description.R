suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(gridExtra))

##########################################################################################################################
# The simulation should handle many mediators, even the entire set of a certain omics layer (e.g. all the protein levels).
# Here, we (rather naively) determine the distributions of mediator heritabilities and the number of QTLs to sample from.
##########################################################################################################################

# As proxy to heritabilities, we use the R2 of each eGene's strongest eQTL from the eQTLGen whole-blood gene expression data.
# The distribution is similar to the variance explained by conditionally significant protein QTLs in 
# Sun et. al 2018 (Fig 1e in Genomic atlas of the human plasma protein), measured in adjusted R2. 
# It roughly resembles a Weibull distribution with shape 0.5 and scale 0.05, truncated to [0, 1].

# As proxy to the mediator QTL number distribution, we use the number of significant variants remaining
# after pruning eQTLGen's cis-eQTLs a certain distance from lead variants (for each gene separately).

# Overlay actual and generated heritabilities in ggplot for a sanity check
visually_verify_h2_dist <- function(dat, mediator_heritabilities_pdf, left_truncate = 0, right_truncate = 1) {
    actual_h2 <- dat %>%
        arrange(-abs(Zscore)) %>%
        #filter(row_number() == 1, .by = Gene) %>%
        group_by(Gene) %>%
        filter(row_number() == 1) %>%
        ungroup() %>%
        mutate(R = sign(Zscore) * sqrt(Zscore**2 / (NrSamples + Zscore**2)),
               R2 = R**2)
    
    generated_h2 <- data.frame(x = generate_h2(10000000))
    
    gg <- ggplot() +
        geom_histogram(aes(x = R2, y = stat(density)), data = actual_h2, binwidth = 0.001) +
        geom_density(aes(x = x), data = generated_h2, colour = "red") +
        xlab("Mediator R2") +
        theme_bw() +
        theme(panel.grid.minor = element_blank())
    
    ggsave(mediator_heritabilities_pdf, gg)
}

# In order to simulate the B matrix, determine the distribution of the number of QTLs mediators have
determine_QTL_dist <- function(dat, mediator_QTLs_pdf) {
    dist <- dat %>%
        group_by(Gene) %>%
        summarize(n_eqtls = n()) %>%
        count(n_eqtls) %>%
        mutate(p = prop.table(n),
               cumulative_p = cumsum(p)) %>%
        select(n_eqtls, p, cumulative_p)
    
    pdf(mediator_QTLs_pdf, width = 2.5, height = 2.5)
    grid.table(mutate(dist, p = round(p, 3), cumulative_p = round(cumulative_p, 3)), rows = NULL)
    dev.off()
        
    return(dist)
}

# A quick shortcut to clumping: without considering actual LD, prune everything a certain distance from lead variants
prune <- function(dat, bp_threshold = 500000) {
    dat_prune <- dat %>%
        mutate(reference_pos = SNPPos + bp_threshold) %>%
        arrange(-abs(Zscore))
    
    dat_selected <- data.frame()
    while (nrow(dat_prune) != 0) {
        dat_prune <- dat_prune %>%
            filter(abs(reference_pos - SNPPos) >= bp_threshold) %>%
            #mutate(reference_pos = first(SNPPos), .by = Gene)
            group_by(Gene) %>%
            mutate(reference_pos = first(SNPPos)) %>%
            ungroup()
        
        dat_selected <- dat_prune %>%
            filter(reference_pos == SNPPos) %>%
            bind_rows(dat_selected)
    }
    
    select(dat_selected, Gene, SNP, Zscore, NrSamples)
}

parser <- ArgumentParser()
parser$add_argument("--heritability_tools_R")
parser$add_argument("--in_mediator_summary_data")
parser$add_argument("--in_mediator_individual_data")
parser$add_argument("--out_mediator_heritabilities_pdf")
parser$add_argument("--out_mediator_QTLs_pdf")
parser$add_argument("--out_mediator_QTLs_txt")
parser$add_argument("--out_mediator_Sigma_RDS")
args <- parser$parse_args()

source(args$heritability_tools_R)

### Distributions for mediator heritabilities and QTL counts are currently only based on eQTLGen's cis-eQTL data
dat_summary <- fread(args$in_mediator_summary_data) %>%
    filter(BonferroniP <= 0.05) %>%
    prune(bp_threshold = 250000)

visually_verify_h2_dist(dat_summary, args$out_mediator_heritabilities_pdf) # Verify the <generate_h2> function defined in simulation_tools.R
QTL_dist <- determine_QTL_dist(dat_summary, args$out_mediator_QTLs_pdf) # Determine the distribution for mediator QTL counts

### Correlations between mediators are currently based on Colaus gene expression rpkm data
Sigma <- fread(args$in_mediator_individual_data) %>%
    column_to_rownames("ID") %>%
    filter(rowSums(. >= 0.1) >= 0.5 * (ncol(.) - 1)) %>%
    t() %>%
    cor()

write_delim(QTL_dist, args$out_mediator_QTLs_txt)
saveRDS(Sigma, args$out_mediator_Sigma_RDS)
