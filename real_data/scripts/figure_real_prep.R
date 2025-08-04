library(data.table)
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(gt)
library(scales)
library(patchwork)
library(circlize)
library(latex2exp)
library(ggrepel)
library(ggtext)
library(GGally)
library(MetBrewer)
library(RColorBrewer)
library(metafor)

# Calculate MP, filter and attach annotations (exposure and outcome name to the data)
prepare_data_with_annotations <- function(filename, dat_annotations) {
    dat_annotated <- fread(filename, colClasses = list(character = c("exposure_ID", "outcome_ID")), fill = TRUE) %>%
        drop_na(total, alpha) %>%
        mutate(MP_lower = MP - 1.96*SE_MP,
               MP_upper = MP + 1.96*SE_MP) %>%
        merge(select(dat_annotations, exposure_ID = Field_ID, exposure = Abbreviation, exposure_annotation = Annotation, exposure_color = color), by = "exposure_ID") %>%
        merge(select(dat_annotations, outcome_ID = Field_ID, outcome = Abbreviation, outcome_annotation = Annotation, outcome_color = color), by = "outcome_ID")
    
    return(dat_annotated)
}

forest_data <- function(mediators = "Lotta_et_al_2021", mediator_pval_threshold_from_exposure = 0.05, mediator_pval_threshold_to_outcome = 0.05, 
                        mediator_cor_threshold = 0.1, n_sig = 1, with_selected_mediators = TRUE) {
    mediators_basename <- paste0(mediators, "_5E-8_", mediator_pval_threshold_from_exposure, "_", mediator_pval_threshold_to_outcome, "_", mediator_cor_threshold, "_", "true.txt")
    
    mr_exposure_to_mediators <- fread(paste0("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_", mediators, "_5E-8.txt"), select = c("exposure", "outcome", "pval"), colClasses = c(exposure = "character", outcome = "character"), col.names = c("exposure_ID", "mediator_ID", "pval"))
    if (mediators == "Lotta_et_al_2021")
        map <- Lotta_map
    else if (mediators == "INTERVAL_plasma_proteins") {
        map <- data.table(mediator_ID = unique(mr_exposure_to_mediators$mediator_ID)) %>%
            mutate(mediator = sub("\\..*", "", mediator_ID))
    }
    
    if (with_selected_mediators) {
        selected_mediators <- fread(paste0("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mediation_by_method/", mediators_basename), colClasses = c(exposure_ID = "character", outcome_ID = "character")) %>% 
            filter(method == "naive") %>%
            select(exposure_ID, outcome_ID, starts_with("gamma")) %>%
            pivot_longer(cols = starts_with("gamma"),
                         names_to = "mediator_ID",
                         names_prefix = "gamma_",
                         values_to = "gamma") %>%
            drop_na(gamma) %>%
            select(-gamma) %>%
            merge(map, by = "mediator_ID", all.x = TRUE) %>%
            merge(mr_exposure_to_mediators, by = c("exposure_ID", "mediator_ID"), all.x = TRUE)
    }
    else
        selected_mediators = data.frame()
    
    dat_forest <- prepare_data_with_annotations(paste0("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mediation_by_method/", mediators_basename), dat_annotations) %>%
        filter(method %in% c("naive", "integrated_fixed_contained")) %>%
        mutate(SE_total = ifelse(is.na(SE_total), SE_total_naive, SE_total),
               SE_alpha = ifelse(is.na(SE_alpha), SE_alpha_naive, SE_alpha),
               SE_MP = alpha^2/total^2 * (SE_alpha^2/alpha^2 + SE_total^2/total^2),
               P_MP = 2 * pnorm(-abs(MP / SE_MP)),
               MP_lower = MP - 1.96*SE_MP,
               MP_upper = MP + 1.96*SE_MP,
               is_significant = MP > 0 & P_MP <= 0.05 / nrow(distinct(., exposure, outcome)),
               is_significant_nominal = MP > 0 & P_MP <= 0.05,
               line_size = ifelse(is_significant_nominal, 1.2, 0.5),
               point_size = ifelse(is_significant_nominal, 2.5, 1),
               alpha_color = ifelse(is_significant_nominal, 1, 0.1)) %>%
        filter(sum(is_significant & between(MP, 0, 1) & k_selected >= 1) >= n_sig, .by = c(exposure, outcome)) %>%
        arrange(exposure_annotation, exposure, outcome_annotation, outcome)
    
    not_converged <- dat_forest %>%
        filter(!(convergence == 0 | method %in% c("naive", "naive_zero", "Burgess") &
                     sigmag2_convergence %in% c(0, NA) &
                     sigmad2_convergence %in% c(0, NA)) | 
                   (sign(alpha) != sign(total) & method == "integrated_fixed_contained")) %>%
        select(exposure_ID, outcome_ID)
    
    dat_forest <- anti_join(dat_forest, not_converged, by = c("exposure_ID", "outcome_ID"))
    
    list(dat_forest = dat_forest, selected_mediators = selected_mediators)
}

gtab_traits <- function(dat_annotations) {
    gtab <- dat_annotations %>%
        select(Annotation, `UKBB trait` = Description, Abbreviation, `UKBB Field ID` = Field_ID) %>%
        arrange(Annotation, Abbreviation) %>%
        mutate(Annotation = ifelse(duplicated(Annotation), "", toupper(as.character(Annotation)))) %>%
        gt() %>%
        data_color(columns = Abbreviation,
                   target_columns = everything(),
                   colors = col_factor(dat_annotations$color, levels = dat_annotations$Abbreviation)) %>%
        cols_width(Annotation ~ px(220),
                   `UKBB trait` ~ px(300),
                   Abbreviation ~ px(150),
                   `UKBB Field ID` ~ px(150)) %>%
        tab_style(locations = cells_body(columns = everything()), style = cell_borders(sides = c("top"), weight = px(0))) %>%
        tab_style(locations = cells_body(columns = Annotation), style = list(cell_text(color = "white", size = px(16)), cell_borders(sides = c("right"), weight = px(7), style = "double"))) %>%
        tab_style(locations = cells_body(columns = c(`UKBB trait`, Abbreviation, `UKBB Field ID`)), style = cell_text(color = "white")) %>%
        tab_style(locations = cells_body(columns = c(Abbreviation, `UKBB Field ID`)), style = list(cell_text(align = "center"), cell_borders(sides = c("left"), weight = px(2)))) %>%
        tab_style(locations = cells_column_labels(), style = list(cell_text(align = "center", weight = "bold"), cell_borders(sides = c("bottom"), weight = px(3)))) %>%
        tab_style(locations = cells_body(rows = nrow(dat_annotations)), style = cell_borders(sides = c("bottom"), weight = px(3))) %>%
        tab_options(table.border.top.style = "hidden",
                    table.font.names = FONT_FAMILY,
                    table.font.size = 14,
                    column_labels.font.size = 20,
                    data_row.padding = 6)
    
    return(gtab)
}

gg_MR <- function(dat_MR) {
    dat_hlines <- distinct(dat_MR, exposure_annotation, exposure) %>%
        group_by(exposure_annotation) %>%
        mutate(i = row_number() - 0.5) %>%
        rbind(summarize(., exposure = "", i = n() + 0.5))
    dat_vlines <- distinct(dat_MR, outcome_annotation, outcome) %>%
        group_by(outcome_annotation) %>%
        mutate(i = row_number() - 0.5) %>%
        rbind(summarize(., outcome = "", i = n() + 0.5))
    
    gg <- ggplot() +
        geom_point(data = dat_MR, aes(y = fct_rev(exposure), x = outcome), size = 0, color = "white") +
        geom_point(data = filter(dat_MR, is_sig), aes(y = fct_rev(exposure), x = outcome, size = -log10(pval)), alpha = 0.8) +
        geom_point(data = filter(dat_MR, is_sig & !is_sig_rev), aes(y = fct_rev(exposure), x = outcome, size = -log10(pval)), color = "firebrick3", show.legend = FALSE) +
        geom_hline(data = dat_hlines, aes(yintercept = i), linewidth = 0.1, alpha = 0.1) +
        geom_vline(data = dat_vlines, aes(xintercept = i), linewidth = 0.1, alpha = 0.1) +
        facet_grid(rows = vars(exposure_annotation), cols = vars(outcome_annotation), switch = "y", scales = "free", space = "free") +
        scale_x_discrete(position = "top", expand = expansion(add = c(0.505, 0.505))) +
        scale_y_discrete(expand = expansion(add = c(0.505, 0.505))) +
        scale_size(name = TeX("$-log_{10}(P):$"), range = c(1, 4), breaks = c(-log10(5e-8), -log10(2e-308)), labels = c("5e-8", "2e-308"), limits = c(-log10(5e-8), -log10(2e-308))) +
        theme_classic() +
        theme(legend.position = c(-0.095, 1.055),
              legend.title = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
              legend.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
              panel.spacing = unit(0.2, "cm"),
              strip.background = element_rect(linewidth = 0.1),
              strip.placement = "outside",
              strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(10, 10, 10, 10), color = "white"),
              strip.text.y.left = element_text(angle = 0),
              axis.title = element_blank(),
              axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY, angle = 45, hjust = 0),
              axis.text.y = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
        guides(size = guide_legend(title.position = "left"))
    
    gtab <- ggplot_gtable(ggplot_build(gg))
    strips <- which(grepl("strip-", gtab$layout$name))
    colors <- c(unique(dat_MR$exposure_color), unique(dat_MR$outcome_color))
    
    for (k in 1:length(strips)) {
        i <- strips[k]
        j <- which(grepl("rect", gtab$grobs[[i]]$grobs[[1]]$childrenOrder))
        gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colors[k]
        gtab$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[k]
    }
    
    return(gtab)
}

gg_Fig3 <- function(dat_forest, selected_mediators, dat_annotations) {
    gg_forest_Fig3 <- function(dat_forest, dat_fills) {
        gg <- ggplot(dat_forest, aes(x = MP, y = fct_rev(factor(outcome)), color = method, alpha = alpha_color)) +
            geom_rect(data = dat_fills, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
            geom_vline(xintercept = 0, linewidth = 0.5, alpha = 0.5) +
            geom_linerange(aes(xmin = MP_lower, xmax = MP_upper), linewidth = 0.6, position = position_dodge(width = 0.8)) +
            geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
            facet_grid(exposure ~ ., scales = "free", space = "free", switch = "y") +
            coord_cartesian(xlim = c(-0.25, 1)) +
            labs(tag = "Exposure  Outcome") + #Exposure   Outcome
            scale_x_continuous(name = TeX("\\widehat{MP}"), expand = c(0, 0), breaks = c(0, 0.2, 0.5, 1), labels = c("0", "0.2", "0.5", "1")) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_colour_manual(values = colors) +
            scale_fill_identity() +
            scale_alpha_identity() +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag.position = c(0.13, 0.025), #c(0.12, 0.015)
                  plot.tag = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  #plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold", margin = margin(0, 8, 8, 0)),
                  panel.spacing = unit(0.4, "cm"),
                  panel.border = element_rect(linewidth = 0.5, fill = "transparent", color = "black"),
                  strip.background.y = element_rect(linewidth = 0.1),
                  strip.placement = "outside",
                  strip.text.y.left = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE, color = "white", angle = 0, margin = margin(0, 10, 0, 10)),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
                  axis.text.y = element_markdown(family = FONT_FAMILY, size = AXIS_TEXT_SIZE, fill = "black", color = "white", margin = margin(0, 10, 0, 10), padding = unit(c(4, 2, 2, 2), "pt")),
                  axis.text.x = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE))
        
        return(gg)
    }
    
    gg_k_Fig3 <- function(dat_forest, dat_fills) {
        dat_k <- filter(dat_forest, method == "naive")
        
        gg <- ggplot(dat_k, aes(x = k_selected, y = fct_rev(factor(outcome)))) +
            geom_rect(data = dat_fills, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white") +
            geom_col(width = 0.6) +
            facet_grid(exposure ~ ., scales = "free", space = "free", switch = "y") +
            coord_cartesian(xlim = c(0, NA)) +
            scale_x_continuous(name = "Nr of mediators", expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            theme_classic() +
            theme(panel.spacing = unit(0.4, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_blank(),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
                  axis.text.y = element_blank(),
                  axis.text.x = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE),
                  axis.ticks.y = element_blank())
        
        return(gg)
    }
    
    gg_labels_Fig3 <- function(dat_forest, dat_fills, selected_mediators, top = 5) {
        dat_labels <- distinct(dat_forest, exposure_ID, exposure, outcome_ID, outcome) %>%
            merge(selected_mediators, by = c("exposure_ID", "outcome_ID")) %>%
            arrange(pval) %>%
            group_by(exposure, outcome) %>%
            mutate(n = n()) %>%
            slice_head(n = top) %>%
            summarize(mediators = paste0(mediator, collapse = ", "), n = unique(n), .groups = "drop") %>%
            mutate(mediators = ifelse(n > top, paste0(mediators, ", ..."), mediators))
        
        gg <- ggplot(dat_labels, aes(x = 0, y = fct_rev(factor(outcome)), label = mediators)) +
            geom_rect(data = dat_fills, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "white") +
            geom_text(hjust = 0, family = FONT_FAMILY, size = (AXIS_TEXT_SIZE-2) * 5/14) +
            facet_grid(exposure ~ ., scales = "free", space = "free", switch = "y") +
            coord_cartesian(xlim = c(0, 1)) +
            scale_x_continuous(name = "Mediators in the model", expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            theme_classic() +
            theme(panel.spacing = unit(0.4, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_blank(),
                  axis.line.y = element_blank(),
                  axis.title.y = element_blank(),
                  axis.title.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
                  axis.text = element_blank(),
                  axis.ticks = element_blank())
        
        return(gg)
    }
    
    gg_legend_Fig3 <- function(dat_annotations) {
        dat_dummy_annotation <- dat_annotations %>%
            filter(Field_ID %in% dat_forest$exposure_ID | Field_ID %in% dat_forest$outcome_ID) %>%
            distinct(Annotation, color)
        dat_dummy_method <- tibble(x = c(1, 2, 1, 2), y = c(1, 2, 1, 2), method = c("MR+MVMR", "MR+MVMR", "I-LiMA", "I-LiMA"))
        
        gg <- ggplot() +
            geom_bar(data = dat_dummy_annotation, aes(x = 1, y = 1, fill = Annotation), position = "dodge", stat = "identity") +
            geom_line(data = dat_dummy_method, aes(x = x, y = y, color = method), linewidth = 1) +
            geom_point(data = dat_dummy_method, aes(x = x, y = y, color = method), size = 3) +
            scale_fill_manual(values = dat_dummy_annotation$color) +
            scale_color_manual(values = c("#D45626", "#9C9C9C")) +
            guides(fill = guide_legend(title = NULL, label.position = "top", order = 1, keyheight = unit(0.8, "cm"), label.vjust = -8.5, nrow = 1,
                                       label.theme = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE, color = "white", margin = margin(5, 10, 5, 10))),
                   color = guide_legend(title = NULL, order = 2, keywidth = unit(2, "cm"),
                                        label.theme = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE))) +
            theme(legend.direction = "horizontal",
                  legend.background = element_rect(fill = "transparent"),
                  legend.box.just = "center",
                  legend.box.margin = margin(0, 0, 20, 0),
                  legend.key = element_rect(fill = "transparent", linewidth = 0))
        
        gg_legend <- ggpubr::get_legend(gg) %>%
            ggpubr::as_ggplot()
        
        return(gg_legend)
    }
    
    dat_fills <- distinct(dat_forest, exposure, outcome) %>%
        mutate(fill = rep(c("#f5f5f5", "white"), length.out = n())) %>%
        mutate(id = n() - row_number() + 1, .by = exposure) %>%
        mutate(ymin = id - 0.5,
               ymax = id + 0.5,
               xmin = -Inf,
               xmax = Inf) %>%
        merge(distinct(dat_forest, exposure, outcome, mediator_type), by = c("exposure", "outcome")) %>%
        arrange(exposure, desc(id))
    
    gg_forest <- gg_forest_Fig3(dat_forest, dat_fills)
    gg_k <- gg_k_Fig3(dat_forest, dat_fills)
    gg_labels <- gg_labels_Fig3(dat_forest, dat_fills, selected_mediators)
    
    gg_patch <- gg_forest + gg_k + gg_labels + plot_layout(widths = c(3, 1, 2)) #c(4, 1, 4)
    
    gtab <- patchworkGrob(gg_patch)
    strip <- which(grepl("strip-l", gtab$layout$name))[1]
    colors <- distinct(select(dat_forest, exposure, exposure_color))$exposure_color
    for (i in 1:length(colors)) {
        j <- which(grepl("rect", gtab$grobs[[strip]]$grobs[[i]]$grobs[[1]]$childrenOrder))
        gtab$grobs[[strip]]$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- colors[i]
        gtab$grobs[[strip]]$grobs[[i]]$grobs[[1]]$children[[j]]$gp$col <- colors[i]
    }
    
    axis <- which(grepl("axis-l", gtab$layout$name))[1]
    colors <- distinct(select(dat_forest, exposure, outcome, outcome_color))$outcome_color
    colored <- 0
    for (i in 1:length(gtab$grobs[[axis]]$grobs)) {
        n_labels <- length(gtab$grobs[[axis]]$grobs[[i]]$children[[2]]$grobs[[1]]$children)
        for (j in 1:n_labels) {
            gtab$grobs[[axis]]$grobs[[i]]$children[[2]]$grobs[[1]]$children[[j]]$children[[1]]$gp$fill <- colors[colored + n_labels - j + 1]
            gtab$grobs[[axis]]$grobs[[i]]$children[[2]]$grobs[[1]]$children[[j]]$children[[1]]$width <- unit(48.1449731445313, "points")
            gtab$grobs[[axis]]$grobs[[i]]$children[[2]]$grobs[[1]]$children[[j]]$children[[1]]$x <- unit(-54.3764770507813, "points")
        }
        colored = colored + n_labels
    }
    
    gg_legend <- gg_legend_Fig3(dat_annotations)
    
    gg <- gg_legend / ggplotify::as.ggplot(gtab) + plot_layout(heights = c(1, 5)) #c(1, 8)
    
    return(gg)
}

gg_scatter_Fig4 <- function(dat_x, dat_y, tag, axis_title_x, axis_title_y, color_x = "dodgerblue3", color_y = "firebrick3") {
    dat <- merge(dat_x, dat_y, by = c("exposure", "outcome", "method"), suffixes = c("_x", "_y")) %>%
        filter(method == "naive") %>%
        mutate(alpha = ifelse(is_significant_x | is_significant_y, 0.6, 0.05),
               color = ifelse(is_significant_x & !is_significant_y, color_x,
                              ifelse(!is_significant_x & is_significant_y, color_y, "black")))
    
    gg <- ggplot(dat, aes(x = MP_x, y = MP_y, color = color, alpha = alpha)) +
        annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, alpha = 0.1, fill = color_y) +
        annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, alpha = 0.1, fill = color_x) +
        annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.5, fill = "#f5f5f5") +
        geom_linerange(y = 0, xmin = -Inf, xmax = 0, linewidth = 0.1, color = "#f5f5f5") +
        geom_linerange(x = 0, ymin = 0, ymax = Inf, linewidth = 0.1, color = "#f5f5f5") +
        geom_linerange(y = 0, xmin = 0, xmax = Inf, linewidth = 0.1, color = "#f5f5f5") +
        geom_linerange(x = 0, ymin = -Inf, ymax = 0, linewidth = 0.1, color = "#f5f5f5") +
        geom_abline(intercept = 0, slope = 1, linewidth = 0.1, alpha = 0.5) +
        #geom_smooth(inherit.aes = FALSE, aes(x = MP_x, y = MP_y), method = lm, se = FALSE, color = "#970E53", linewidth = 1, fullrange = TRUE) +
        geom_linerange(aes(xmin = MP_lower_x, xmax = MP_upper_x), linewidth = 0.6) +
        geom_linerange(aes(ymin = MP_lower_y, ymax = MP_upper_y), linewidth = 0.6) +
        geom_point(size = 2.5) +
        #annotate("text", x = 0.53, y = 0.53, hjust = 1, label = paste0("n = ", sum(dat$is_significant_y & dat$is_significant_x)), color = "black", size = 5/14 * AXIS_TITLE_SIZE, fontface = "bold") +
        #annotate("text", x = -0.2, y = 0.53, hjust = 0, label = paste0("n = ", sum(dat$is_significant_y & !dat$is_significant_x)), color = color_y, size = 5/14 * AXIS_TITLE_SIZE, fontface = "bold") +
        #annotate("text", x = 0.53, y = -0.2, hjust = 1, label = paste0("n = ", sum(!dat$is_significant_y & dat$is_significant_x)), color = color_x, size = 5/14 * AXIS_TITLE_SIZE, fontface = "bold") +
        coord_cartesian(xlim = c(-0.2, 0.53), ylim = c(-0.2, 0.53)) +
        scale_x_continuous(axis_title_x, breaks = c(0, 0.3), labels = c("0", "0.3")) +
        scale_y_continuous(axis_title_y, breaks = c(0, 0.3), labels = c("0", "0.3")) +
        scale_alpha_identity() +
        scale_color_identity() +
        labs(tag = tag) +
        theme_classic() +
        theme(plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
              axis.title = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              axis.text = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE))
    
    return(gg)
}

gg_histogram_Fig4 <- function(dat_x, dat_y, legend_title, color_x = "dodgerblue3", color_y = "firebrick3", legend_position = c(0.35, 0.77)) {
    dat_histogram <- rbind(dat_x, dat_y) %>%
        filter(method == "naive") %>%
        mutate(k = pmin(k_selected, 10))
    
    k_max <- max(dat_histogram$k_selected)
    breaks <- seq(1, max(dat_histogram$k))
    labels <- as.character(breaks)
    if (k_max == 10)
        labels[k_max] <- " 10"
    if (k_max > 10)
        labels[10] <- "  10+"
    
    gg <- ggplot(dat_histogram, aes(x = k, y = after_stat(density), fill = type, color = type)) +
        geom_histogram(position = "stack", binwidth = 0.5, alpha = 0.5) +
        scale_x_continuous("Nr of mediators", breaks = breaks, labels = labels) +
        scale_y_continuous("Density", expand = c(0, 0)) +
        scale_fill_manual(legend_title, values = c(color_y, color_x)) +
        scale_color_manual(legend_title, values = c(color_y, color_x)) +
        theme_classic() +
        theme(legend.position = legend_position,
              legend.justification = c("left", "bottom"),
              legend.title = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              legend.text = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE),
              legend.background = element_rect(fill = "transparent"),
              legend.key.width = unit(0.5, "cm"),
              legend.key.height = unit(0.5, "cm"),
              axis.ticks.y = element_blank(),
              axis.text.y = element_blank(),
              axis.text.x = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE),
              axis.title.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              axis.title.y = element_blank(),
              axis.line.y = element_blank())
    
    return(gg)
}

gg_Fig4 <- function(dat_forest) {
    gg_methods_Fig4 <- function(dat_forest, tag) {
        dat_naive <- filter(dat_forest, method == "naive") %>%
            select(exposure, outcome, MP, MP_lower, MP_upper, is_significant_nominal)
        dat_LiMA <- filter(dat_forest, method == "integrated_fixed_contained") %>%
            select(exposure, outcome, MP, MP_lower, MP_upper, is_significant_nominal)
        
        dat_methods <- merge(dat_naive, dat_LiMA, by = c("exposure", "outcome"), suffixes = c("_naive", "_integrated")) %>%
            mutate(alpha = ifelse(is_significant_nominal_naive & is_significant_nominal_integrated, 1, 0.05))
        
        gg <- ggplot(dat_methods, aes(x = MP_naive, y = MP_integrated, alpha = alpha)) +
            annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = "#f5f5f5", alpha = 0.5) +
            geom_linerange(x = 0, ymin = 0, ymax = Inf, linewidth = 0.1, color = "#f5f5f5") +
            geom_linerange(y = 0, xmin = 0, xmax = Inf, linewidth = 0.1, color = "#f5f5f5") +
            geom_abline(intercept = 0, slope = 1, linewidth = 0.1, alpha = 0.5) +
            geom_linerange(aes(xmin = MP_lower_naive, xmax = MP_upper_naive), linewidth = 0.8) +
            geom_linerange(aes(ymin = MP_lower_integrated, ymax = MP_upper_integrated), linewidth = 0.8) +
            geom_point(size = 3) +
            coord_cartesian(xlim = c(-0.1, 0.5), ylim = c(-0.1, 0.5)) +
            scale_x_continuous(TeX("$\\widehat{MP}$ of MR+MVMR"), expand = c(0, 0), breaks = c(0, 0.4), labels = c("0", "0.4")) +
            scale_y_continuous(TeX("$\\widehat{MP}$ of I-LiMA"), expand = c(0, 0), breaks = c(0, 0.4), labels = c("0", "0.4")) +
            scale_alpha_identity() +
            labs(tag = tag) +
            theme_classic() +
            theme(plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  axis.title = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
                  axis.text = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE))
        
        return(gg)
    }
    
    gg_methods <- gg_methods_Fig4(dat_forest, "a")
    
    dat_P_strict <- forest_data(mediator_pval_threshold_from_exposure = "Bonferroni", mediator_pval_threshold_to_outcome = "Bonferroni", n_sig = 0, with_selected_mediators = FALSE)$dat_forest %>%
        select(exposure, outcome, method, MP, MP_lower, MP_upper, is_significant = is_significant_nominal, k_selected) %>%
        mutate(type = "Strict")
    dat_P_relaxed <- forest_data(n_sig = 0, with_selected_mediators = FALSE)$dat_forest %>%
        select(exposure, outcome, method, MP, MP_lower, MP_upper, is_significant = is_significant_nominal, k_selected) %>%
        mutate(type = "Relaxed")
    gg_scatter_filtering <- gg_scatter_Fig4(dat_P_strict, dat_P_relaxed, "b", TeX("$\\widehat{MP}$ from P-strict filtering"), TeX("$\\widehat{MP}$ from P-relaxed filtering"), "#0f7ba2", "#dd5129")
    gg_histogram_filtering <- gg_histogram_Fig4(dat_P_strict, dat_P_relaxed, "P-filtering:", "#0f7ba2", "#dd5129")
    
    dat_LD_strict <- forest_data(n_sig = 0, with_selected_mediators = FALSE)$dat_forest %>%
        select(exposure, outcome, method, MP, SE_MP, MP_lower, MP_upper, is_significant = is_significant_nominal, k_selected) %>%
        mutate(type = "Strict")
    dat_LD_relaxed <- forest_data(mediator_cor_threshold = 0.5, n_sig = 0, with_selected_mediators = FALSE)$dat_forest %>%
        select(exposure, outcome, method, MP, SE_MP, MP_lower, MP_upper, is_significant = is_significant_nominal, k_selected) %>%
        mutate(type = "Relaxed")
    gg_scatter_LD <- gg_scatter_Fig4(dat_LD_strict, dat_LD_relaxed, "c", TeX("$\\widehat{MP}$ from LD-strict filtering"), TeX("$\\widehat{MP}$ from LD-relaxed filtering"), "#0f7ba2", "#dd5129")
    gg_histogram_LD <- gg_histogram_Fig4(dat_LD_strict, dat_LD_relaxed, "LD-filtering:", "#0f7ba2", "#dd5129")
    
    gg <- (gg_methods | ((gg_scatter_filtering + gg_histogram_filtering + plot_layout(widths = c(2, 1))) /
                             (gg_scatter_LD + gg_histogram_LD + plot_layout(widths = c(2, 1))))) + plot_layout(widths = c(7, 5))
    
    return(gg)
}

dat_annotations <- read.xlsx("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/data/raw/conf/UKBB_Neale_traits_ALL.xlsx") %>%
    filter(Include) %>%
    arrange(Annotation, Abbreviation) %>%
    mutate(Abbreviation = factor(Abbreviation, levels = Abbreviation),
           Annotation = factor(Annotation, levels = unique(Annotation))) %>%
    mutate(color = c(met.brewer(name = "Juarez", n = n_distinct(.$Annotation), type = "continuous"))[cur_group_id()], .by = Annotation)

Lotta_map <- fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/data/raw/metabolites/Lotta_et_al_2021/metabolites_label.txt", header = FALSE) %>%
    #mutate(mediator = ifelse(V4 %in% c("Lyso-phosphatidylcholines", "Phosphatidylcholines", "Sphingomyelins"), V2, V3)) %>%
    select(mediator_ID = V1, mediator = V3)

##### Supplementary figure about UKBB complex traits
FONT_FAMILY <- "Helvetica"

gtab <- gtab_traits(dat_annotations)
print(gtab)
#gtsave(gtab, filename = "UKBB_traits.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper", zoom = 6)


##### Supplementary figure about UKBB pairwise MR
FONT_FAMILY <- "Helvetica"
AXIS_TITLE_SIZE <- 12
AXIS_TEXT_SIZE <- 10

dat_MR <- rbind(fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_UKBB_binary_5E-8.txt", col.names = c("exposure_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")),
                fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_UKBB_quantitative_5E-8.txt", col.names = c("exposure_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character"))) %>%
    merge(select(dat_annotations, exposure_ID = Field_ID, exposure = Abbreviation, exposure_annotation = Annotation, exposure_color = color), by = "exposure_ID") %>%
    merge(select(dat_annotations, outcome_ID = Field_ID, outcome = Abbreviation, outcome_annotation = Annotation, outcome_color = color), by = "outcome_ID") %>%
    merge(select(., exposure, outcome, pval_rev = pval), by.x = c("exposure", "outcome"), by.y = c("outcome", "exposure"), all.x = TRUE) %>%
    mutate(pval = ifelse(pval == 0, 2e-308, pval),
           is_sig = pval <= 5e-8,
           is_sig_rev = ifelse(is.na(pval_rev), FALSE, pval_rev <= 5e-8))

gtab <- gg_MR(dat_MR)
grid::grid.draw(gtab)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/MR.pdf", gtab, width = 14, height = 10)


##### Main application Fig3, forest plot with MP estimates and selected mediators
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12
colors = c("#D45626", "#9C9C9C")

forest_data_Fig3 <- forest_data()
gg_Fig3(forest_data_Fig3$dat_forest, forest_data_Fig3$selected_mediators, dat_annotations)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/Fig3.pdf", width = 15.5, height = 10.2)

forest_data_Fig3_INTERVAL <- forest_data(mediators = "INTERVAL_plasma_proteins")
gg_Fig3(forest_data_Fig3_INTERVAL$dat_forest, forest_data_Fig3_INTERVAL$selected_mediators, dat_annotations)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/Fig3_INTERVAL.pdf", width = 12.5, height = 6.4)


##### Main application Fig4, general comparisons between methods and filtering strategies
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 16
AXIS_TEXT_SIZE <- 14

gg_Fig4(forest_data_Fig3$dat_forest)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/Fig4.pdf", width = 16, height = 9)


##### Supplementary .xlsx tables
prepare_MP_sup_data <- function(forest_data_sup, mediator_pval_threshold_from_exposure, mediator_pval_threshold_to_outcome) {
    selected_mediators <- forest_data_sup$selected_mediators %>%
        arrange(pval) %>%
        group_by(exposure_ID, outcome_ID) %>%
        summarize(mediators = paste0(mediator, collapse = ", "), .groups = "drop")
    
    dat_sup <- forest_data_sup$dat_forest %>%
        merge(selected_mediators, by = c("exposure_ID", "outcome_ID")) %>%
        select(exposure_annotation, exposure, exposure_ID, mediator_type, mediators, outcome_annotation, outcome, outcome_ID,
               instrument_pval_threshold, mediator_cor_threshold, k = k_selected, l, m, nX, nZ, nY, method,
               total, alpha, SE_total, SE_alpha, MP, SE_MP, P_MP) %>%
        mutate(method = case_when(method == "naive" ~ "MR+MVMR",
                                  method == "integrated_fixed_contained" ~ "I-LiMA",
                                  TRUE ~ method)) %>%
        mutate(mediator_pval_threshold_from_exposure = mediator_pval_threshold_from_exposure,
               mediator_pval_threshold_to_outcome = mediator_pval_threshold_to_outcome, .before = mediator_cor_threshold) %>%
        arrange(exposure_annotation, exposure, outcome_annotation, outcome, desc(method))
    
    return(dat_sup)
}

prepare_MP_sup_data_with_indirect_effects <- function(forest_data_sup, mediator_pval_threshold_from_exposure, mediator_pval_threshold_to_outcome) {
    dat_sup <- forest_data_sup$dat_forest %>%
        select(exposure_annotation, exposure, exposure_ID, mediator_type, outcome_annotation, outcome, outcome_ID,
               instrument_pval_threshold, mediator_cor_threshold, k = k_selected, l, m, nX, nM=nZ, nY, method,
               total, alpha, SE_total, SE_alpha, MP, SE_MP, P_MP, starts_with("gamma"), starts_with("delta"))
    
    selected_mediators <- dat_sup %>%
        pivot_longer(cols = c(starts_with("gamma"), starts_with("delta")), names_to = c(".value", "mediator_ID"), names_pattern = "([^_]*)_(.*)") %>%
        merge(forest_data_sup$selected_mediators, by = c("exposure_ID", "outcome_ID", "mediator_ID")) %>%
        filter(method == "naive") %>%
        arrange(pval) %>%
        group_by(exposure_ID, outcome_ID, method) %>%
        summarize(mediators = paste0(mediator, collapse = ", "), 
                  gammas = paste0(format(gamma, digits = 5), collapse = ", "),
                  deltas = paste0(format(delta, digits = 5), collapse = ", "), .groups = "drop")
    
    dat_sup <- select(dat_sup, -starts_with("gamma_"), -starts_with("delta_")) %>%
        left_join(selected_mediators) %>%
        mutate(method = case_when(method == "naive" ~ "MR+MVMR",
                                  method == "integrated_fixed_contained" ~ "I-LiMA",
                                  TRUE ~ method)) %>%
        mutate(mediator_pval_threshold_from_exposure = mediator_pval_threshold_from_exposure,
               mediator_pval_threshold_to_outcome = mediator_pval_threshold_to_outcome, .before = mediator_cor_threshold) %>%
        arrange(exposure_annotation, exposure, outcome_annotation, outcome, desc(method))
    
    return(dat_sup)
}

Lotta_default <- prepare_MP_sup_data_with_indirect_effects(forest_data("Lotta_et_al_2021", 0.05, 0.05, 0.1, n_sig = 0), 0.05, 0.05)
Lotta_P_strict <- prepare_MP_sup_data(forest_data("Lotta_et_al_2021", "Bonferroni", "Bonferroni", 0.1, n_sig = 0), "Bonferroni", "Bonferroni")
Lotta_LD_relaxed <- prepare_MP_sup_data(forest_data("Lotta_et_al_2021", 0.05, 0.05, 0.5, n_sig = 0), 0.05, 0.05)
INTERVAL_default <- prepare_MP_sup_data_with_indirect_effects(forest_data("INTERVAL_plasma_proteins", 0.05, 0.05, 0.1, n_sig = 0), 0.05, 0.05)

MP_result_tables <- list("Lotta_default" = Lotta_default, "Lotta_P_strict" = Lotta_P_strict, "Lotta_LD_relaxed" = Lotta_LD_relaxed, "INTERVAL_default" = INTERVAL_default)
#write.xlsx(MP_result_tables, file = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/result_tables/MP_results.xlsx")


UKBB_to_UKBB <- rbind(fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_UKBB_binary_5E-8.txt", col.names = c("exposure_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")),
                      fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_UKBB_quantitative_5E-8.txt", col.names = c("exposure_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character"))) %>%
    merge(select(dat_annotations, exposure_ID = Field_ID, exposure = Abbreviation, exposure_annotation = Annotation), by = "exposure_ID") %>%
    merge(select(dat_annotations, outcome_ID = Field_ID, outcome = Abbreviation, outcome_annotation = Annotation), by = "outcome_ID") %>%
    select(exposure_annotation, exposure, exposure_ID, outcome_annotation, outcome, outcome_ID, nsnp, b, se, pval) %>%
    arrange(exposure_annotation, exposure, outcome_annotation, outcome)

UKBB_to_Lotta <- fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_Lotta_et_al_2021_5E-8.txt", col.names = c("exposure_ID", "mediator_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")) %>%
    merge(select(dat_annotations, exposure_ID = Field_ID, exposure = Abbreviation, exposure_annotation = Annotation), by = "exposure_ID") %>%
    merge(Lotta_map, by = "mediator_ID") %>%
    mutate(outcome_annotation = "Lotta_et_al_2021") %>%
    select(exposure_annotation, exposure, exposure_ID, outcome_annotation, outcome = mediator, outcome_ID = mediator_ID, nsnp, b, se, pval)

UKBB_to_INTERVAL <- fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/UKBB_quantitative_INTERVAL_plasma_proteins_5E-8.txt", col.names = c("exposure_ID", "mediator_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")) %>%
    merge(select(dat_annotations, exposure_ID = Field_ID, exposure = Abbreviation, exposure_annotation = Annotation), by = "exposure_ID") %>%
    mutate(mediator = sub("\\..*", "", mediator_ID),
           outcome_annotation = "INTERVAL_plasma_proteins") %>%
    select(exposure_annotation, exposure, exposure_ID, outcome_annotation, outcome = mediator, outcome_ID = mediator_ID, nsnp, b, se, pval)

Lotta_to_UKBB <- rbind(fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/Lotta_et_al_2021_UKBB_binary_5E-8.txt", col.names = c("mediator_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")),
                       fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/Lotta_et_al_2021_UKBB_quantitative_5E-8.txt", col.names = c("mediator_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character"))) %>%
    merge(Lotta_map, by = "mediator_ID") %>%
    merge(select(dat_annotations, outcome_ID = Field_ID, outcome = Abbreviation, outcome_annotation = Annotation), by = "outcome_ID") %>%
    mutate(exposure_annotation = "Lotta_et_al_2021") %>%
    select(exposure_annotation, exposure = mediator, exposure_ID = mediator_ID, outcome_annotation, outcome, outcome_ID, nsnp, b, se, pval)
    
INTERVAL_to_UKBB <- rbind(fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/INTERVAL_plasma_proteins_UKBB_binary_5E-8.txt", col.names = c("mediator_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character")),
                          fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/results/urblauna/mr_ivw/INTERVAL_plasma_proteins_UKBB_quantitative_5E-8.txt", col.names = c("mediator_ID", "outcome_ID", "nsnp", "b", "se", "pval"), colClasses = c(exposure = "character"))) %>%
    merge(select(dat_annotations, outcome_ID = Field_ID, outcome = Abbreviation, outcome_annotation = Annotation), by = "outcome_ID") %>%
    mutate(mediator = sub("\\..*", "", mediator_ID),
           exposure_annotation = "INTERVAL_plasma_proteins") %>%
    select(exposure_annotation, exposure = mediator, exposure_ID = mediator_ID, outcome_annotation, outcome, outcome_ID, nsnp, b, se, pval)

MR_IVW_tables <- list("UKBB_to_UKBB" = UKBB_to_UKBB, "UKBB_to_Lotta" = UKBB_to_Lotta, "UKBB_to_INTERVAL" = UKBB_to_INTERVAL, "Lotta_to_UKBB" = Lotta_to_UKBB, "INTERVAL_to_UKBB" = INTERVAL_to_UKBB)
#write.xlsx(MR_IVW_tables, file = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/real_data/figures/paper/result_tables/MR_IVW_results.xlsx")
