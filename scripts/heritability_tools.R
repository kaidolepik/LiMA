####################################################################################################
# Auxiliary functions to generate mediator heritabilities (in mediator_description.R and simulate.R)
####################################################################################################

# If we used the R2 of each eGene's strongest eQTL from the eQTLGen whole-blood gene expression as proxies to heritabilities,
# the distribution of mediator heritabilities roughly resembles a Weibull distribution with shape 0.5 and scale 0.05, truncated to [0, 1].

##### GENERATE MEDIATOR HERITABILITIES (modify accordingly) #####
generate_untruncated_h2 <- function(k) {
    rweibull(k, shape = 0.5, scale = 0.05)
}

##### GENERATE APPROPRIATELY TRUNCATED HERITABILITIES #####
generate_h2 <- function(k, left_truncate = 0, right_truncate = 1) {
    h2 <- generate_untruncated_h2(k)
    
    out_of_range <- !between(h2, left_truncate, right_truncate)
    while (sum(out_of_range) > 0) {
        h2[out_of_range] <- generate_untruncated_h2(sum(out_of_range))
        out_of_range <- !between(h2, left_truncate, right_truncate)
    }
    
    return(h2)
}
# For the moment, we only consider h2 through B; additional heritability explained through X would be small anyway
