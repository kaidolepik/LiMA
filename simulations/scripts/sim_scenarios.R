suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(tools))

gather_sim_scenarios <- function(filename, cols, N_sims) {
    sim_name <- file_path_sans_ext(basename(filename))
    
    # Set the initial seed
    initial_seed <- charToRaw(sim_name) %>%
        readBin(what = "integer")
    
    sim_scenarios <- getSheetNames(filename) %>% 
        lapply(function(sheet) {
            read.xlsx(filename, sheet = sheet, cols = cols) %>% 
                   mutate(sheet = sheet, .before = 1)
            }) %>%
        rbindlist() %>%
        distinct() %>%
        mutate(seed = initial_seed + row_number(),
               N_mediators = as.integer(k*k_sig_perc + k*(1-k_sig_perc)*p_value_mediators),
               piece_size = case_when(k <= 500 & N_mediators <= 30 ~ N_sims,
                                      k <= 1000 & N_mediators <= 50 ~ 100,
                                      k <= 3000 & N_mediators <= 100 ~ 10,
                                      TRUE ~ 1),
               N_pieces = N_sims / piece_size, .before = 1) %>%
        uncount(N_pieces, .remove = FALSE) %>% # Duplicate rows by N_pieces
        group_by_all() %>%
        mutate(piece = 1:n()) %>% # Row/piece number
        ungroup() %>%
        mutate(start = piece_size*(piece-1) + 1,
               end = piece_size*piece,
               sim_name = sim_name, .after = 1) %>%
        select(-N_mediators, -piece_size, -N_pieces, -piece)
    
    sim_scenarios$out_simulation_RDS <- paste0(sapply(1:nrow(sim_scenarios), function(i) paste0(unlist(sim_scenarios[i, ]), collapse = "_")), ".RDS")
    
    return(sim_scenarios)
}

parser <- ArgumentParser()
parser$add_argument("--N_sims", type = "double")
parser$add_argument("--in_sims_xlsx")
parser$add_argument("--out_scenarios_txt")
parser$add_argument("--out_sims_txt")
args <- parser$parse_args()

# Read all the simulation scenarios into a data frame
sim_scenarios <- gather_sim_scenarios(args$in_sims_xlsx, cols = 1:16, N_sims = args$N_sims)

# Add the method which is used to solve the sim scenario
sims <- crossing(method = c("naive", "original", "integrated"), sim_scenarios) %>%
    mutate(out_simulation_RDS = paste0(method, "_", out_simulation_RDS))

write_delim(sim_scenarios, args$out_scenarios_txt, delim = "\t")
write_delim(sims, args$out_sims_txt, delim = "\t")
    