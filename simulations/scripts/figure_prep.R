library(tidyverse)
library(data.table)
library(magrittr)
library(patchwork)
library(MetBrewer)
library(latex2exp)
library(ggh4x)
library(ggpubr)
library(ggtext)
library(gt)
library(scales)
library(RColorBrewer)

gg_line_supplementary <- function(dat_gg, x_var, label_var, x_title, x_breaks, x_labels, y_limits, y_breaks, x_annotate_MP, y_annotate_MP, x_min, tag, colors) {
    y_intercept <- 1 - unique(dat_gg$p_direct)
    
    dat_rect <- tibble(xmin = x_min, xmax = Inf, ymin = -Inf, ymax = Inf)
    dat_text <- tibble(x = x_annotate_MP, y = y_annotate_MP, label = "True MP")
    
    gg <- ggplot(dat_gg, aes(x = get(x_var), y = slope, colour = get(label_var))) +
        geom_rect(data = dat_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.03, color = NA, inherit.aes = FALSE) +
        geom_hline(yintercept = y_intercept, linetype = "solid", linewidth = 0.5, alpha = 0.5) +
        geom_point(size = 3) +
        geom_line(linewidth = 1.4) +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 1.4, width = 0.04*max(dat_gg[x_var])) +
        geom_text(data = dat_text, aes(x = x, y = y, label = label), alpha = 0.8, size = 5, colour = "black", family = FONT_FAMILY, inherit.aes = FALSE) +
        coord_cartesian(ylim = y_limits) +
        scale_y_continuous(TeX("$\\widehat{MP}$"), breaks = y_breaks) +
        scale_x_continuous(x_title, breaks = x_breaks, labels = x_labels) +
        scale_colour_manual(values = colors) +
        labs(tag = tag) +
        theme_bw() +
        theme(legend.position = "none",
              axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
              axis.title = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
              panel.grid.minor = element_blank(),
              panel.grid.major = element_line(linewidth = 0.4),
              plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE))
    
    print(gg)
    
    return(gg)
}

gg_legend_Figs <- function(dat, colors, margin_r = 1.5) {
    gg_dummy <- ggplot(dat, aes(x = 1, y = 1, fill = method_label_succinct)) +
        geom_bar(position = "dodge", stat = "identity") +
        scale_fill_manual(values = colors) +
        theme(legend.title = element_blank(),
              legend.direction = "horizontal",
              legend.text = element_text(size = 20, family = FONT_FAMILY, margin = margin(r = margin_r, unit = "cm")),
              legend.background = element_rect(fill = "transparent"),
              legend.key = element_rect(fill = "transparent", linewidth = 0),
              legend.key.width = unit(0.7, "cm"),
              legend.key.height = unit(0.6, "cm"))
    
    gg_legend <- get_legend(gg_dummy) %>%
        as_ggplot()
    
    return(gg_legend)
}

gg_Fig1 <- function(dat_Fig1, dat_time, colors) {
    gg_sample_sizes_Fig1 <- function(dat_Fig1, tag, colors) {
        dat_nZ <- filter(dat_Fig1, purpose == "nZ", nZ %in% c(1000, 3000, 10000, 30000))
        dat_nX <- filter(dat_Fig1, purpose == "nX", nX %in% c(5000, 10000, 30000, 50000))
        dat_nY <- filter(dat_Fig1, purpose == "nY", nY %in% c(5000, 10000, 30000, 50000))
        
        mediator_panel_coef <- 1.67
        
        dat_gg <- rbind(dat_nX, dat_nZ, dat_nY) %>%
            mutate(purpose = factor(purpose, levels = c("nZ", "nX", "nY")), # To order the panels
                   x = case_when(purpose == "nX" ~ nX, 
                                 purpose == "nZ" ~ nZ,
                                 purpose == "nY" ~ nY)) %>%
            mutate(x = x + ifelse(purpose == "nZ", 0.34, 1) * ifelse(method == "original", 0, ifelse(method == "naive", -650, 650)),
                   width = ifelse(purpose == "nZ", 1/mediator_panel_coef, 1) * 0.04 * max(x), .by = purpose)
        dat_text <- tibble(purpose = factor("nY"), x = 10000, y = 0.162, label = "True MP")
        
        gg <- ggplot(dat_gg, aes(x = x, y = slope, colour = method_label_succinct)) +
            facet_wrap(vars(purpose), strip.position = "bottom", scales = "free_x", labeller = labeller(purpose = c(nZ = "Mediator sample size", nX = "Exposure sample size", nY = "Outcome sample size"))) +
            geom_hline(yintercept = 1 - unique(dat_gg$p_direct), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1.4) +
            geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, width = width), linewidth = 1.4) +
            geom_point(size = 3) +
            geom_text(data = dat_text, aes(x = x, y = y, label = label), alpha = 0.8, size = 5, colour = "black", family = FONT_FAMILY, inherit.aes = FALSE) +
            coord_cartesian(xlim = c(0, NA), ylim = c(0, NA)) +
            scale_y_continuous(TeX("$\\widehat{MP}$"), breaks = c(0.05, 0.15, 0.25), expand = expansion(mult = c(0, 0.05))) +
            scale_colour_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY)) +
            force_panelsizes(cols = c(mediator_panel_coef, 1, 1)) +
            facetted_pos_scales(x = list(
                purpose == "nZ" ~ scale_x_continuous(breaks = c(1000, 3000, 10000, 30000), labels = c("1000    ", "    3000", "10000", "30000"), expand = expansion(mult = c(0, 0.02))),
                purpose == "nX" ~ scale_x_continuous(breaks = c(5000, 10000, 30000, 50000), labels = c("5000     ", "     10000", "30000", "50000"), expand = expansion(mult = c(0, 0.04))),
                purpose == "nY" ~ scale_x_continuous(breaks = c(5000, 10000, 30000, 50000), labels = c("5000     ", "     10000", "30000", "50000"), expand = expansion(mult = c(0, 0.04)))
            ))
        
        return(gg)
    }
    
    gg_bar_Fig1 <- function(dat_Fig1, tag, colors) {
        dat_p_value_mediators <- filter(dat_Fig1, purpose == "p-value-mediators") %>%
            mutate(panel = factor("widehat(MP)"), # Plotmath of TeX("$\\widehat{MP}$") is widehat(MP)
                   y = slope, 
                   y_intercept = 1 - p_direct)
        dat_k_selected <- mutate(dat_p_value_mediators, panel = factor("Mediators"), y = k_selected, y_intercept = k * k_sig_perc)
        
        dat_hline <- rbind(dat_k_selected, dat_p_value_mediators) %>%
            summarize(y_intercept = unique(y_intercept), .by = panel)
        dat_hline_best <- data.frame(panel = factor("widehat(MP)"), y_intercept = max(dat_p_value_mediators$slope))
        dat_texts <- rbind(tibble(panel = factor("widehat(MP)"), x = 0.75, y = 0.16, label = "True MP"),
                           tibble(panel = factor("Mediators"), x = 4.5, y = 48, label = "Number of non-zero mediators"))
        
        gg <- ggplot() +
            facet_grid(panel ~ ., scales = "free_y", switch = "y", labeller = label_parsed) + 
            geom_hline(data = dat_hline, aes(yintercept = y_intercept), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
            geom_hline(data = dat_hline_best, aes(yintercept = y_intercept), linetype = "dashed", linewidth = 0.5, alpha = 0.5, colour = colors[3]) +
            geom_bar(data = dat_p_value_mediators, aes(x = factor(p_value_mediators), y = y, fill = method_label_succinct), position = "dodge", stat = "identity", width = 0.7) +
            geom_errorbar(data = dat_p_value_mediators, aes(x = factor(p_value_mediators), group = method_label_succinct, y = y, ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.7), linewidth = 0.7, width = 0.4) +
            geom_point(data = dat_k_selected, aes(x = factor(p_value_mediators), y = y), size = 3, colour = "black") +
            geom_line(data = dat_k_selected, aes(x = as.numeric(factor(p_value_mediators)), y = y), linewidth = 1.2, colour = "black") +
            geom_text(data = dat_texts, aes(x = x, y = y, label = label), alpha = 0.8, size = 5, colour = "black", family = FONT_FAMILY) +
            scale_x_discrete("P-value threshold for selecting mediators in the model") +
            scale_fill_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  panel.spacing = unit(0.6, "cm"),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY)) +
            force_panelsizes(rows = c(1, 2)) +
            facetted_pos_scales(y = list(
                panel == "Mediators" ~ scale_y_continuous(limits = c(0, max(dat_k_selected$y) + 10), breaks = c(25, 100)),
                panel == "widehat(MP)" ~ scale_y_continuous(breaks = c(0.05, 0.15), expand = expansion(mult = c(0, 0.05)))
            ))
        
        return(gg)
    }
    
    gg_other_params_Fig1 <- function(dat_Fig1, tag, colors) {
        dat_m <- filter(dat_Fig1, purpose == "m", m >= 10)
        dat_h2X <- filter(dat_Fig1, purpose == "h2X")
        dat_total_causal <- filter(dat_Fig1, purpose == "total-causal", total_causal >= 0.05)
        
        dat_gg <- rbind(dat_m, dat_h2X, dat_total_causal) %>%
            mutate(x = case_when(purpose == "h2X" ~ h2X + ifelse(method == "original", 0, ifelse(method == "naive", -0.013, 0.013)),
                                 purpose == "m" ~ m + ifelse(method == "original", 0, ifelse(method == "naive", -14.5, 14.5)),
                                 purpose == "total-causal" ~ total_causal + ifelse(method == "original", 0, ifelse(method == "naive", -0.0047, 0.0047)))) %>%
            mutate(width = 0.05*max(x), .by = purpose)
        dat_rects <- rbind(tibble(purpose = "h2X", xmin = -Inf, xmax = 0.1, ymin = -Inf, ymax = Inf),
                           tibble(purpose = "h2X", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf),
                           tibble(purpose = "total-causal", xmin = 0.2, xmax = Inf, ymin = -Inf, ymax = Inf))
        dat_text <- tibble(purpose = "total-causal", x = 0.08, y = 0.156, label = "True MP")
        
        gg <- ggplot(dat_gg, aes(x = x, y = slope, colour = method_label_succinct)) +
            facet_wrap(vars(purpose), strip.position = "bottom", scales = "free_x", labeller = labeller(purpose = c(m = "Exposure polygenicity", h2X = "Exposure heritability", `total-causal` = "Expected total causal effect"))) +
            geom_rect(data = dat_rects, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.03, color = NA, inherit.aes = FALSE) +
            geom_hline(yintercept = 1 - unique(dat_gg$p_direct), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1.4) +
            geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, width = width), linewidth = 1.4) +
            geom_point(size = 3) +
            geom_text(data = dat_text, aes(x = x, y = y, label = label), alpha = 0.8, size = 5, colour = "black", family = FONT_FAMILY, inherit.aes = FALSE) +
            coord_cartesian(ylim = c(0, NA)) +
            scale_y_continuous(TeX("$\\widehat{MP}$"), breaks = c(0.05, 0.15), expand = expansion(mult = c(0, 0.05))) +
            scale_colour_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(x = list(
                purpose == "m" ~ scale_x_continuous(breaks = c(10, 100, 200, 400, 800), labels = c("10  ", "100", "   200", "400", "800")),
                purpose == "h2X" ~ scale_x_continuous(breaks = c(0.01, 0.1, 0.2, 0.35, 0.5, 0.7), labels = c("0.01   ", "0.1", "0.2", "0.35", "0.5", "0.7")),
                purpose == "total-causal" ~ scale_x_continuous(breaks = c(0.01, 0.05, 0.1, 0.15, 0.2, 0.3), labels = c("0.01    ", "0.05", "0.1", "0.15", "0.2", "0.3"))
            ))
        
        return(gg)
    }
    
    gg_MP_Fig1 <- function(dat_Fig1, tag, colors) {
        dat_gg <- filter(dat_Fig1, purpose == "p-direct", p_direct %in% c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95)) %>% 
            mutate(true_MP = 1 - p_direct + ifelse(method == "original", 0, ifelse(method == "naive", -0.012, 0.012)),
                   panel = "widehat(MP)")
        dat_rect <- tibble(xmin = 0.3, xmax = Inf, ymin = -Inf, ymax = Inf)
        
        gg <- ggplot(dat_gg, aes(x = true_MP, y = slope, colour = method_label_succinct)) +
            facet_grid(panel ~ ., scales = "free_y", switch = "y", labeller = label_parsed) +
            geom_rect(data = dat_rect, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", alpha = 0.03, color = NA, inherit.aes = FALSE) +
            geom_abline(slope = 1, intercept = 0, linetype = "solid", linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1.4) +
            geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), linewidth = 1.4, width = 0.035*max(dat_gg$true_MP)) +
            geom_point(size = 3) +
            scale_x_continuous("True MP", limits = c(0, 1), expand = c(0, 0), breaks = c(0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95), labels = c("0.05    ", " 0.15", "0.3", "0.5", "0.7", "0.85 ", "    0.95")) +
            scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = c(0.15, 0.5, 0.85), labels = c("0.15", "0.5", "0.85")) +
            scale_colour_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY))
    }
    
    gg_box_Fig1 <- function(dat_time, tag, colors) {
        x_breaks <- unique(dat_time$k_selected)
        
        dat_gg <- mutate(dat_time, panel = factor("Time (sec)"))
        
        gg <- ggplot(dat_gg, aes(x = k_selected, fill = method_label_succinct, group = interaction(k_selected, method_label_succinct))) + 
            facet_grid(panel ~ ., scales = "free_y", switch = "y") + 
            geom_boxplot(aes(ymin = time_0, lower = time_25, middle = time_50, upper = time_75, ymax = time_100), stat = "identity", width = 10) +
            scale_x_continuous("Mediators in the model", breaks = x_breaks) +
            scale_y_continuous("Time (sec)", trans = "log10", expand = expansion(mult = c(0, 0.05)), breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
            scale_fill_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY))
        
        return(gg)
    }
    
    gg_sample_sizes <- gg_sample_sizes_Fig1(dat_Fig1, "a", colors)
    gg_p_value_mediators <- gg_bar_Fig1(dat_Fig1, "b", colors)
    gg_other_params <- gg_other_params_Fig1(dat_Fig1, "c", colors)
    gg_MP <- gg_MP_Fig1(dat_Fig1, "d", colors)
    gg_time <- gg_box_Fig1(dat_time, "e", colors)
    gg_legend <- gg_legend_Figs(dat_Fig1, colors)
    
    gg <- (gg_sample_sizes / 
               gg_legend /
               (gg_p_value_mediators + gg_MP + plot_layout(widths = c(5, 2))) /
               (gg_other_params + gg_time + plot_layout(widths = c(5, 2)))) +
        plot_layout(heights = c(5, 0.9, 5.1, 5))
    
    print(gg)
    
    return(gg)
}

gg_Fig2 <- function(dat_Fig2, dat_metrics_Fig2, dat_metrics_filtered_Fig2, colors) {
    gg_powcovT1E_Fig2 <- function(dat_metrics_Fig2, tag, colors) {
        dat_T1E <- filter(dat_metrics_Fig2, purpose == "T1E", metric == "Power") %>%
            mutate(metric = "Type I error")
        dat_powcov <- filter(dat_metrics_Fig2, purpose == "Truth", metric %in% c("Power", "Coverage"))
        
        dat_gg <- rbind(dat_T1E, dat_powcov)
        dat_hlines <- data.frame(metric = unique(dat_gg$metric), yintercept = c(0.05, 0.8, 0.8))
        
        gg <- ggplot(dat_gg, aes(x = as.factor(k), y = estimate, fill = method_label_succinct)) +
            geom_bar(position = "dodge", stat = "identity", width = 0.7) +
            geom_errorbar(aes(group = method_label_succinct, ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.7), linewidth = 0.7, width = 0.4) +
            geom_hline(data = dat_hlines, aes(yintercept = yintercept), linewidth = 0.5, alpha = 0.5) +
            facet_wrap(vars(metric), nrow = 1, strip.position = "left", scales = "free_y") +
            coord_cartesian(ylim = c(0, 1)) +
            scale_fill_manual(values = colors) +
            xlab("Mediators in the model") +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_markdown(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(y = list(
                metric == "Coverage" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Power" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Type I error" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.05, 1), labels = c("0", paste0("<sup style='font-size: ", AXIS_TEXT_SIZE, "pt'>0.05</sup>"), "1"))
            ))
        
        return(gg)
    }
    
    gg_oracle_Fig2 <- function(dat_Fig2, tag, colors, purpose) {
        dat_gg <- dat_Fig2 %>%
            filter(purpose == !!purpose, default) %>%
            mutate(SE_MP = ifelse(is.nan(SE_MP), SE_MP_naive, SE_MP),
                   MP_lower = MP - 1.96*SE_MP,
                   MP_upper = MP + 1.96*SE_MP,
                   contains_zero = between(0, MP_lower, MP_upper))
        
        x_breaks <- c(-0.1, 0, 0.1)
        x_labels <- c("    -0.1", "0", "0.1   ")
        x_lim <- c(-0.1, 0.1)
        if (purpose == "Truth") {
            x_breaks <- c(0, 1 - unique(dat_gg$p_direct))
            x_labels <- c("0", paste0("     ", as.character(x_breaks[2])))
            x_lim <- c(-x_breaks[2], 3*x_breaks[2])
        }
        
        gg <- ggplot() + 
            geom_vline(xintercept = 0, linewidth = 0.2, alpha = 0.5, linetype = "longdash") +
            geom_vline(xintercept = 1 - unique(dat_gg$p_direct), linewidth = 0.2, alpha = 0.5) +
            geom_linerange(data = filter(dat_gg, contains_zero), aes(y = true_total, xmin = MP_lower, xmax = MP_upper, colour = method_label_succinct), linewidth = 0.6, alpha = 0.2) +
            geom_linerange(data = filter(dat_gg, !contains_zero), aes(y = true_total, xmin = MP_lower, xmax = MP_upper, colour = method_label_succinct), linewidth = 0.8, alpha = 1) +
            geom_linerange(linewidth = 0.3, alpha = 1) +
            facet_wrap(vars(method_label_succinct), nrow = 1) +
            coord_cartesian(xlim = x_lim, ylim = c(-0.15, 0.3)) +
            scale_x_continuous(TeX("95% CI of $\\widehat{MP}$"), expand = c(0, 0), breaks = x_breaks, labels = x_labels) +
            scale_y_continuous(TeX(paste0("True total causal effect $(\\theta)$ of ", ifelse(purpose == "T1E", "null", "oracle"), " model $(k=", unique(dat_gg$k), ")$")), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "0.15")) +
            scale_color_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major = element_blank(),
                  panel.spacing = unit(0.6, "cm"),
                  strip.placement = "outside",
                  strip.text = element_blank(),
                  axis.title = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY))
        
        return(gg)
    }
    
    gg_predictions_Fig2 <- function(dat_Fig2, tag, colors) {
        dat_gg <- dat_Fig2 %>%
            filter(default) %>%
            mutate(Coverage = 1 * between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP),
                   Power = 1 * !between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP))
        
        methods <- unique(dat_gg$method)
        metrics <- c("Predicted coverage", "Predicted power", "Predicted type I error")
        fits <- list()
        for (method in methods) {
            for (metric in metrics) {
                dat_metric_method <- dat_gg %>%
                    filter(method == !!method, purpose == ifelse(metric == "Predicted type I error", "T1E", "Truth")) %>%
                    mutate(metric = !!metric,
                           obs = ifelse(metric == "Predicted coverage", Coverage, Power)) %>%
                    select(metric, true_total, method_label_succinct, obs)
                
                model <- glm(obs ~ I(abs(true_total)), data = dat_metric_method, family = "binomial")
                
                new_totals <- seq(min(dat_metric_method$true_total), max(dat_metric_method$true_total), 0.001)
                fits_model <- predict(model, newdata = list(true_total = new_totals), se.fit = TRUE, type = "response") %>%
                    bind_cols() %>%
                    mutate(true_total = new_totals,
                           CI_lower = fit - 1.96*se.fit,
                           CI_upper = fit + 1.96*se.fit,
                           method_label_succinct = unique(dat_metric_method$method_label_succinct),
                           metric = !!metric)
                
                fits[[length(fits) + 1]] <- fits_model
            }
        }
        dat_fits <- bind_rows(fits)
        dat_hlines <- data.frame(metric = metrics, yintercept = c(0.8, 0.8, 0.05))
        
        gg <- ggplot(dat_fits, aes(x = true_total, y = fit, color = method_label_succinct, fill = method_label_succinct)) +
            geom_hline(data = dat_hlines, aes(yintercept = yintercept), linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1) +
            geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, linewidth = 0) +
            facet_wrap(vars(metric), scales = "free", strip.position = "left") +
            coord_cartesian(ylim = c(0, 1), xlim = c(-0.15, 0.3)) +
            scale_x_continuous(TeX("True total causal effect $(\\theta)$"), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "0.15")) +
            scale_fill_manual(values = colors) +
            scale_color_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_markdown(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(y = list(
                metric == "Predicted coverage" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Predicted power" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Predicted type I error" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.05, 1), labels = c("0", paste0("<sup style='font-size: ", AXIS_TEXT_SIZE, "pt'>0.05</sup>"), "1"))
            ))
        
        return(gg)
    }
    
    gg_biasvar_Fig2 <- function(dat_metrics_Fig2, dat_metrics_filtered_Fig2, tag, colors) {
        dat_gg_bias <- filter(dat_metrics_Fig2, purpose == "Truth", metric %in% c("Bias")) %>%
            mutate(metric = "|Bias|")
        dat_gg_variance <- filter(dat_metrics_filtered_Fig2, purpose == "Truth", metric %in% c("Variance"))
        
        gg <- ggplot() +
            geom_bar(data = dat_gg_bias, aes(x = as.factor(k), y = estimate, fill = method_label_succinct), position = "dodge", stat = "identity", width = 0.7) +
            geom_bar(data = dat_gg_variance, aes(x = as.factor(k), y = estimate, fill = method_label_succinct), position = "dodge", stat = "identity", width = 0.7) +
            geom_errorbar(data = dat_gg_bias, aes(group = method_label_succinct, x = as.factor(k), ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.7), linewidth = 0.7, width = 0.4) +
            geom_errorbar(data = dat_gg_variance, aes(group = method_label_succinct, x = as.factor(k), ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.7), linewidth = 0.7, width = 0.4) +
            facet_wrap(vars(metric), ncol = 1, strip.position = "left", scales = "free_y") +
            xlab("Mediators in the model") +
            coord_cartesian(ylim = c(0, 0.108)) +
            scale_fill_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(y = list(
                metric == "|Bias|" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.1), labels = c("0", "0.1")),
                metric == "Variance" ~ scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.1), labels = c("0", "0.1"))
            ))
        
        return(gg)
    }
    
    gg_lm_Fig2 <- function(dat_Fig2, tag, colors) {
        dat_gg <- dat_Fig2 %>%
            filter(purpose == "Truth", default == TRUE) %>%
            mutate(method_label_succinct = factor(method_label_succinct, labels = sapply(1:n_distinct(method), function(i) paste0(paste0(rep(" ", i-1), collapse = ""), "widehat(alpha)", paste0(rep(" ", i-1), collapse = "")))))
        
        gg <- ggplot(dat_gg, aes(x = total, y = alpha, colour = method_label_succinct)) +
            geom_abline(intercept = 0, slope = 0.85, linewidth = 0.5, alpha = 0.5) +
            geom_linerange(aes(ymin = alpha - SE_alpha, ymax = alpha + SE_alpha), linewidth = 0.6, alpha = 0.4) +
            geom_linerange(aes(xmin = total - SE_total, xmax = total + SE_total), linewidth = 0.6, alpha = 0.4) +
            geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, linewidth = 0.8, fullrange = TRUE) +
            facet_wrap(vars(method_label_succinct), ncol = 1, strip.position = "left", labeller = label_parsed) +
            scale_color_manual(values = colors) +
            scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.15*0.85), labels = c("0", TeX("                     0.15$\\cdot(1-MP)$"))) +
            scale_x_continuous(TeX("$\\widehat{\\theta}$"), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "     0.15")) +
            coord_cartesian(xlim = c(-0.15, 0.45), ylim = c(-0.15, 0.45)) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_text(size = c(AXIS_TEXT_SIZE, AXIS_TEXT_SIZE - 2), family = FONT_FAMILY, angle = 90, hjust = 0.5))
        
        return(gg)
    }
    
    gg_violin_Fig2 <- function(dat_Fig2, dat_metrics_Fig2, colors) {
        dat_gg <- filter(dat_Fig2, purpose == "Truth", default)
        
        gg <- ggplot(dat_gg, aes(y = 1, x = MP, fill = method_label_succinct)) + 
            geom_violin(linewidth = 0.8) +
            geom_vline(xintercept = 0.15, linewidth = 0.5, alpha = 0.5) +
            geom_jitter(alpha = 0.05) +
            facet_wrap(vars(method_label_succinct), ncol = 1, strip.position = "left") +
            coord_cartesian(xlim = c(-0.15, 0.45)) +
            scale_x_continuous(TeX("$\\widehat{MP}$"), limits = c(-10, 10), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "     0.15")) +
            scale_fill_manual(values = colors) +
            theme_classic() +
            theme(legend.position = "none",
                  panel.spacing = unit(0.6, "cm"),
                  strip.placement = "outside",
                  strip.text = element_blank(),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())
        
        dat_violin_slope <- dat_metrics_Fig2 %>%
            filter(purpose == "Truth", metric == "Bias", k == unique(dat_gg$k)) %>%
            arrange(method_label_succinct) %>%
            mutate(group = 1:n()) %>%
            inner_join(ggplot_build(gg)$data[[1]], by = c("group" = "group")) %>%
            filter(abs(slope - x) == min(abs(slope - x)), .by = method_label_succinct)
        
        gg <- gg + 
            geom_linerange(data = dat_violin_slope, aes(x = slope, ymin = y - width/2, ymax = y + width/2), linewidth = 0.8, colour = "grey20")
        
        return(gg)
    }
    
    gg_biasvar <- gg_biasvar_Fig2(dat_metrics_Fig2, dat_metrics_filtered_Fig2, "a", colors)    
    gg_lm <- gg_lm_Fig2(dat_Fig2, "b", colors)
    gg_violin <- gg_violin_Fig2(dat_Fig2, dat_metrics_Fig2, colors)
    gg_powcovT1E <- gg_powcovT1E_Fig2(dat_metrics_Fig2, "e", colors)
    gg_oracle <- gg_oracle_Fig2(dat_Fig2, "c", colors, "Truth")
    gg_null <- gg_oracle_Fig2(dat_Fig2, "", colors, "T1E")
    gg_predictions <- gg_predictions_Fig2(dat_Fig2, "d", colors)
    gg_legend <- gg_legend_Figs(dat_Fig2, colors, 0.7)
    
    gg <- (((gg_biasvar / 
                 gg_legend /
                 (gg_lm + gg_violin)) + plot_layout(heights = c(5, 1, 5))) |
               (((gg_oracle + gg_null) /
                     gg_predictions /
                     gg_powcovT1E) + plot_layout(heights = c(5, 3, 3)))) +
        plot_layout(widths = c(5, 11))
    
    return(gg)
}

gg_sup <- function(dat_sup, highlight_scenario, xlim_upper, colors, zeroed = FALSE) {
    highlight_colors <- distinct(dat_sup, purpose_plotmath, scenario) %>%
        filter(scenario == highlight_scenario)
    other_colors <- distinct(dat_sup, purpose_plotmath, scenario) %>%
        filter(scenario != highlight_scenario)
    vlines <- distinct(dat_sup, purpose_plotmath, scenario) %>%
        filter(purpose_plotmath != "MP") %>%
        mutate(x = 0.15)
    ablines <- filter(dat_sup, purpose_plotmath == "MP") %>%
        mutate(a = 0, b = 20)
    ablines2 <- mutate(dat_sup, a = 1.5, b = 0)
    ablines3 <- mutate(dat_sup, a = 2.5, b = 0)
    ablines4 <- mutate(dat_sup, a = 3.5, b = 0)
    
    gg <- ggplot(dat_sup, aes(x = slope, y = x, color = method_label_succinct)) +
        geom_rect(data = highlight_colors, inherit.aes = FALSE, fill = "#f5f5f5", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 1) +
        geom_rect(data = other_colors, inherit.aes = FALSE, fill = "#fcfcfc", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 1) +
        geom_vline(data = vlines, aes(xintercept = x), linewidth = 0.5, alpha = 0.5) +
        geom_abline(data = ablines, aes(intercept = a, slope = b), linewidth = 0.5, alpha = 0.5) +
        geom_abline(data = ablines2, aes(intercept = a, slope = b), linewidth = 0.1, alpha = 0.15) +
        geom_abline(data = ablines3, aes(intercept = a, slope = b), linewidth = 0.1, alpha = 0.15) +
        geom_abline(data = ablines4, aes(intercept = a, slope = b), linewidth = 0.1, alpha = 0.15) +
        geom_linerange(aes(xmin = CI_lower, xmax = CI_upper), linewidth = 1, position = position_dodge(width = 0.8)) +
        geom_point(size = 1.5, position = position_dodge(width = 0.8)) +
        scale_x_continuous(TeX("$\\widehat{MP}$"), limits = c(NA, NA), breaks = c(0, 0.15, 0.3), labels = c("0", "0.15", "0.3  ")) +
        coord_cartesian(xlim = c(0, xlim_upper)) +
        scale_colour_manual(values = colors) +
        theme_classic() +
        theme(legend.position = "top",
              legend.title = element_blank(),
              legend.text = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE, margin = margin(r = 1.5, unit = "cm")),
              legend.key.width = unit(1.2, "cm"),
              legend.key.height = unit(0.6, "cm"),
              strip.placement = "outside",
              strip.background = element_blank(),
              strip.text.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              strip.text.y = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              strip.text.y.left = element_text(angle = 0),
              axis.title.y = element_blank(),
              axis.title.x = element_text(family = FONT_FAMILY, size = AXIS_TITLE_SIZE),
              axis.text.y = ggtext::element_markdown(family = FONT_FAMILY, size = AXIS_TEXT_SIZE),
              axis.text.x = element_text(family = FONT_FAMILY, size = AXIS_TEXT_SIZE),
              panel.grid.major.y = element_line(linewidth = 0.2),
              panel.spacing = unit(0.3, "cm"))
    
    if (zeroed)
        gg <- gg + facet_wrap(vars(purpose_plotmath), nrow = 6, ncol = 2, dir = "v", strip.position = "left", scales = "free_y", labeller = labeller(purpose_plotmath = label_parsed, scenario = waiver()))
    else
        gg <- gg + facet_grid(purpose_plotmath ~ scenario, switch = "y", scales = "free_y", labeller = labeller(purpose_plotmath = label_parsed, scenario = waiver()))
    
    print(gg)
    
    return(gg)
}

gg_Zhu <- function(dat_regression, dat_filtered, methods, scenario, pleiotropy, colors, y_text) {
    gg_sample_sizes <- function(dat_Zhu_bias, tag, colors, y_text) {
        dat_nZ <- filter(dat_Zhu_bias, purpose == "nZ", nZ %in% c(3000, 10000, 30000))
        dat_nX <- filter(dat_Zhu_bias, purpose == "nX", nX %in% c(10000, 30000, 50000))
        dat_nY <- filter(dat_Zhu_bias, purpose == "nY", nY %in% c(10000, 30000, 50000))
        
        mediator_panel_coef <- 1.67
        
        dat_gg <- rbind(dat_nX, dat_nZ, dat_nY) %>%
            mutate(purpose = factor(purpose, levels = c("nZ", "nX", "nY")), # To order the panels
                   x = case_when(purpose == "nX" ~ nX, 
                                 purpose == "nZ" ~ nZ,
                                 purpose == "nY" ~ nY)) %>%
            mutate(x = x + ifelse(purpose == "nZ", 0.34, 1) * ifelse(!(method %in% c("integrated_fixed_contained", "naive")), 0, ifelse(method == "naive", -650, 650)),
                   width = ifelse(purpose == "nZ", 1/mediator_panel_coef, 1) * 0.07 * max(x), .by = purpose)
        dat_text <- tibble(purpose = factor("nZ"), x = 22000, y = y_text, label = "True MP")
        
        gg <- ggplot(dat_gg, aes(x = x, y = slope, colour = method_label_succinct)) +
            facet_wrap(vars(purpose), strip.position = "bottom", scales = "free_x", labeller = labeller(purpose = c(nZ = "Mediator sample size", nX = "Exposure sample size", nY = "Outcome sample size"))) +
            geom_hline(yintercept = 1 - unique(dat_gg$p_direct), linetype = "solid", linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1.4) +
            geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, width = width), linewidth = 1.4) +
            geom_point(size = 3) +
            geom_text(data = dat_text, aes(x = x, y = y, label = label), alpha = 0.8, size = 5, colour = "black", family = FONT_FAMILY, inherit.aes = FALSE) +
            coord_cartesian(xlim = c(0, NA), ylim = c(NA, NA)) +
            scale_y_continuous(TeX("$\\widehat{MP}$"), breaks = c(-0.5, 0, 0.15, 0.5), labels = c("-0.5", "0", "0.15", "0.5"), expand = expansion(mult = c(0.04, 0.04))) +
            scale_colour_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY)) +
            force_panelsizes(cols = c(mediator_panel_coef, 1, 1)) +
            facetted_pos_scales(x = list(
                purpose == "nZ" ~ scale_x_continuous(breaks = c(3000, 10000, 30000), labels = c("3000", "10000", "30000"), expand = expansion(mult = c(0, 0.04))),
                purpose == "nX" ~ scale_x_continuous(breaks = c(10000, 30000, 50000), labels = c("10000", "30000", "50000"), expand = expansion(mult = c(0, 0.04))),
                purpose == "nY" ~ scale_x_continuous(breaks = c(10000, 30000, 50000), labels = c("10000", "30000", "50000"), expand = expansion(mult = c(0, 0.04)))
            ))
        
        return(gg)
    }
    
    gg_time <- function(dat_Zhu_time, tag, colors) {
        x_breaks <- unique(dat_Zhu_time$k_selected)
        
        dat_gg <- mutate(dat_Zhu_time, panel = factor("Time (sec)"))
        
        gg <- ggplot(dat_gg, aes(x = k_selected, fill = method_label_succinct, group = interaction(k_selected, method_label_succinct))) + 
            facet_grid(panel ~ ., scales = "free_y", switch = "y") + 
            geom_boxplot(aes(ymin = time_0, lower = time_25, middle = time_50, upper = time_75, ymax = time_100), stat = "identity", width = 10) +
            scale_x_continuous("Mediators in the model", breaks = x_breaks) +
            scale_y_continuous("Time (sec)", trans = "log10", limits = c(NA, 10), expand = expansion(mult = c(0, 0)), breaks = c(0.01, 0.1, 1, 10), labels = c("0.01", "0.1", "1", "10")) +
            scale_fill_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  axis.text = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  plot.tag = element_text(family = FONT_FAMILY, face = "bold", size = PLOT_TAG_SIZE),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY))
        
        return(gg)
    }
    
    gg_metrics <- function(dat_Zhu_metrics, tag, colors) {
        dat_T1E <- filter(dat_Zhu_metrics, purpose == "T1E", metric == "Power") %>%
            mutate(metric = "Type I error") %>%
            filter(k %in% c(1, 5, 10))
        dat_bias <- filter(dat_Zhu_metrics, purpose == "Truth", metric %in% c("Bias")) %>%
            mutate(metric = ifelse(metric == "Bias", "|Bias|", metric)) %>%
            filter(k %in% c(5, 10, 25))
        dat_powcov <- filter(dat_Zhu_metrics, purpose == "Truth", metric %in% c("Power", "Coverage")) %>%
            filter(k %in% c(5, 10, 25))
        
        dat_gg <- rbind(dat_T1E, dat_bias, dat_powcov)
        dat_hlines <- data.frame(metric = c("Type I error", "Coverage", "Power"), yintercept = c(0.05, 0.8, 0.8))
        
        gg <- ggplot(dat_gg, aes(x = as.factor(k), y = estimate, fill = method_label_succinct)) +
            geom_bar(position = "dodge", stat = "identity", width = 0.7) +
            geom_errorbar(aes(group = method_label_succinct, ymin = CI_lower, ymax = CI_upper), position = position_dodge(width = 0.7), linewidth = 0.7, width = 0.4) +
            geom_hline(data = dat_hlines, aes(yintercept = yintercept), linewidth = 0.5, alpha = 0.5) +
            facet_wrap(vars(metric), nrow = 1, strip.position = "left", scales = "free") +
            scale_fill_manual(values = colors) +
            xlab("Mediators in the model") +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_markdown(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(y = list(
                metric == "|Bias|" ~ scale_y_continuous(limits = c(0, max(0.1, max(abs(dat_bias$CI_upper)))), expand = c(0, 0), breaks = c(0, 0.1, 0.2, 0.3), labels = c("0", "0.1", "0.2", "0.3")),
                metric == "Coverage" ~ scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Power" ~ scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Type I error" ~ scale_y_continuous(limits = c(0, 1), expand = expansion(add = c(0, 0)), breaks = c(0, 0.05, 1), labels = c("0", paste0("<sup style='font-size: ", AXIS_TEXT_SIZE, "pt'>0.05</sup>"), "1"))
            ))
        
        return(gg)
    }
    
    gg_predictions <- function(dat_Zhu, tag, colors) {
        dat_gg <- dat_Zhu %>%
            filter(default) %>%
            mutate(Coverage = 1 * between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP),
                   Power = 1 * !between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP))
        
        methods <- unique(dat_gg$method)
        metrics <- c("Predicted coverage", "Predicted power", "Predicted type I error")
        fits <- list()
        for (method in methods) {
            for (metric in metrics) {
                dat_metric_method <- dat_gg %>%
                    filter(method == !!method, purpose == ifelse(metric == "Predicted type I error", "T1E", "Truth")) %>%
                    mutate(metric = !!metric,
                           obs = ifelse(metric == "Predicted coverage", Coverage, Power)) %>%
                    select(metric, true_total, method_label_succinct, obs)
                
                model <- glm(obs ~ I(abs(true_total)), data = dat_metric_method, family = "binomial")
                
                new_totals <- seq(min(dat_metric_method$true_total), max(dat_metric_method$true_total), 0.001)
                fits_model <- predict(model, newdata = list(true_total = new_totals), se.fit = TRUE, type = "response") %>%
                    bind_cols() %>%
                    mutate(true_total = new_totals,
                           CI_lower = fit - 1.96*se.fit,
                           CI_upper = fit + 1.96*se.fit,
                           method_label_succinct = unique(dat_metric_method$method_label_succinct),
                           metric = !!metric)
                
                fits[[length(fits) + 1]] <- fits_model
            }
        }
        dat_fits <- bind_rows(fits)
        dat_hlines <- data.frame(metric = metrics, yintercept = c(0.8, 0.8, 0.05))
        
        gg <- ggplot(dat_fits, aes(x = true_total, y = fit, color = method_label_succinct, fill = method_label_succinct)) +
            geom_hline(data = dat_hlines, aes(yintercept = yintercept), linewidth = 0.5, alpha = 0.5) +
            geom_line(linewidth = 1) +
            geom_ribbon(aes(ymin = CI_lower, ymax = CI_upper), alpha = 0.5, linewidth = 0) +
            facet_wrap(vars(metric), scales = "free", strip.position = "left") +
            coord_cartesian(ylim = c(0, 1), xlim = c(-0.15, 0.3)) +
            scale_x_continuous(TeX("True total causal effect $(\\theta)$"), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "0.15")) +
            scale_fill_manual(values = colors) +
            scale_color_manual(values = colors) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.minor = element_blank(),
                  panel.grid.major.y = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_markdown(size = AXIS_TEXT_SIZE, family = FONT_FAMILY)) +
            facetted_pos_scales(y = list(
                metric == "Predicted coverage" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Predicted power" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.8, 1), labels = c("0", "0.8", "1")),
                metric == "Predicted type I error" ~ scale_y_continuous(expand = expansion(add = c(0, 0)), breaks = c(0, 0.05, 1), labels = c("0", paste0("<sup style='font-size: ", AXIS_TEXT_SIZE, "pt'>0.05</sup>"), "1"))
            ))
        
        return(gg)
    }
    
    gg_lm <- function(dat_Zhu, tag, colors) {
        dat_gg <- dat_Zhu %>%
            filter(purpose == "Truth", default == TRUE) %>%
            mutate(method_label_succinct = factor(method_label_succinct, labels = sapply(1:n_distinct(method), function(i) paste0(paste0(rep(" ", i-1), collapse = ""), "widehat(alpha)", paste0(rep(" ", i-1), collapse = "")))))
        
        gg <- ggplot(dat_gg, aes(x = total, y = alpha, colour = method_label_succinct)) +
            geom_abline(intercept = 0, slope = 0.85, linewidth = 0.5, alpha = 0.5) +
            geom_linerange(aes(ymin = alpha - SE_alpha, ymax = alpha + SE_alpha), linewidth = 0.6, alpha = 0.4) +
            geom_linerange(aes(xmin = total - SE_total, xmax = total + SE_total), linewidth = 0.6, alpha = 0.4) +
            geom_smooth(method = lm, formula = y ~ 0 + x, se = FALSE, linewidth = 0.8, fullrange = TRUE) +
            facet_wrap(vars(method_label_succinct), ncol = 1, strip.position = "left", labeller = label_parsed) +
            scale_color_manual(values = colors) +
            scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.15*0.85), labels = c("0", TeX("                   0.15$\\cdot(1-MP)$"))) +
            scale_x_continuous(TeX("$\\widehat{\\theta}$"), expand = c(0, 0), breaks = c(0, 0.15), labels = c("0", "0.15")) +
            coord_cartesian(xlim = c(-0.15, 0.45), ylim = c(-0.15, 0.45)) +
            labs(tag = tag) +
            theme_classic() +
            theme(legend.position = "none",
                  plot.tag = element_text(size = PLOT_TAG_SIZE, family = FONT_FAMILY, face = "bold"),
                  panel.grid.major = element_line(linewidth = 0.2),
                  panel.spacing = unit(0.6, "cm"),
                  strip.background = element_blank(),
                  strip.placement = "outside",
                  strip.text = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY, margin = margin(0, 0, 0, 0)),
                  axis.title.x = element_text(size = AXIS_TITLE_SIZE, family = FONT_FAMILY),
                  axis.title.y = element_blank(),
                  axis.text.x = element_text(size = AXIS_TEXT_SIZE, family = FONT_FAMILY),
                  axis.text.y = element_text(size = c(AXIS_TEXT_SIZE, AXIS_TEXT_SIZE - 2), family = FONT_FAMILY, angle = 90, hjust = 0.5))
        
        return(gg)
    }
    
    dat_Zhu_bias <- filter(dat_regression, method %in% methods, scenario == !!scenario, pleiotropy == !!pleiotropy)
    dat_Zhu_time <- dat_filtered %>%
        filter(method %in% methods, scenario == !!scenario, pleiotropy == !!pleiotropy,
               m < 800, k_selected %in% c(1, 25, 50)) %>%
        summarize(k_selected = unique(k_selected),
                  N = n(),
                  time_0 = min(time),
                  time_25 = quantile(time, 0.25),
                  time_50 = quantile(time, 0.5),
                  time_75 = quantile(time, 0.75),
                  time_100 = max(time),
                  .by = c(method_label_succinct, k_selected))
    
    dat_Zhu <- dat_filtered %>%
        filter(method %in% methods, 
               scenario == !!scenario, 
               pleiotropy == !!pleiotropy,
               purpose %in% c("T1E", "Truth"), 
               k %in% c(1, 5, 10, 25))
    dat_Zhu_metrics <- calculate_metrics_Fig2(dat_Zhu)
    
    gg1 <- gg_sample_sizes(dat_Zhu_bias, "a", colors, y_text = y_text)
    gg2 <- gg_metrics(dat_Zhu_metrics, "b", colors)
    gg3 <- gg_predictions(dat_Zhu, "c", colors)
    gg4 <- gg_lm(dat_Zhu, "d", colors)
    gg5 <- gg_time(dat_Zhu_time, "e", colors)
    gg_legend <- gg_legend_Figs(dat_Zhu, colors)
    
    gg <- (((gg1 / gg_legend / gg2 / gg3) + plot_layout(heights = c(5, 1, 5, 5))) |
        (gg4 / gg5) + plot_layout(heights = c(5, 2))) + 
        plot_layout(widths = c(4, 1))
    
    return(gg)
}

calculate_metrics_sup <- function(dat) {
    dat_metrics <- dat %>%
        summarize(slope = 1 - (t(total) %*% alpha / t(total) %*% total)[1],
                  bias_alpha = mean(alpha - true_alpha),
                  bias_total = mean(total - true_total),
                  bias_MP = slope - unique(1 - p_direct),
                  bias_alpha_prc = mean((alpha - true_alpha) / true_alpha),
                  bias_total_prc = mean((total - true_total) / true_total),
                  bias_MP_prc = bias_MP / unique(1 - p_direct),
                  variance_alpha = sum((alpha - true_alpha)^2) / (n() - 1),
                  variance_total = sum((total - true_total)^2) / (n() - 1),
                  variance_MP = sum((MP - slope)^2) / (n() - 1),
                  coverage_alpha = mean(between(true_alpha, alpha - 1.96*SE_alpha, alpha + 1.96*SE_alpha)),
                  coverage_total = mean(between(true_total, total - 1.96*SE_total, total + 1.96*SE_total)),
                  coverage_MP = mean(between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
                  power_alpha = mean(!between(0, alpha - 1.96*SE_alpha, alpha + 1.96*SE_alpha)),
                  power_total = mean(!between(0, total - 1.96*SE_total, total + 1.96*SE_total)),
                  power_MP = mean(!between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
                  not_converged = mean(convergence != 0 | optim_value > 9.99e+307),
                  sigmag2_not_converged = mean(k_selected > 1 & sigmag2_convergence != 0),
                  sigmad2_not_converged = mean(k_selected > 1 & sigmad2_convergence != 0),
                  sigmag2_on_the_bound = mean(k_selected > 1 & abs(sigmag2 - 1e-8) < 1e-8),
                  sigmad2_on_the_bound = mean(k_selected > 1 & abs(sigmad2 - 1e-8) < 1e-8),
                  # SE_Bias = summary(lm(alpha ~ 0 + total))$coef[1, 2],
                  # bias_MP_CI_lower = ((slope - 1.96*SE_Bias) - unique(1 - p_direct)) / unique(1 - p_direct),
                  # bias_MP_CI_upper = ((slope + 1.96*SE_Bias) - unique(1 - p_direct)) / unique(1 - p_direct),
                  # SE_Coverage = sqrt(var(between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP)) / n()),
                  # coverage_MP_CI_lower = coverage_MP - 1.96*SE_Coverage,
                  # coverage_MP_CI_upper = coverage_MP + 1.96*SE_Coverage,
                  # SE_Power = sqrt(var(!between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP)) / n()),
                  # power_MP_CI_lower = power_MP - 1.96*SE_Power,
                  # power_MP_CI_upper = power_MP + 1.96*SE_Power,
                  .groups = "drop") %>%
        pivot_wider(id_cols = !slope, names_from = method, values_from = c(bias_alpha, bias_total, bias_MP, 
                                                                           bias_alpha_prc, bias_total_prc, bias_MP_prc,
                                                                           variance_alpha, variance_total, variance_MP,
                                                                           coverage_alpha, coverage_total, coverage_MP, 
                                                                           power_alpha, power_total, power_MP, 
                                                                           not_converged, sigmag2_not_converged, sigmad2_not_converged, 
                                                                           sigmag2_on_the_bound, sigmad2_on_the_bound), names_glue = "{method}_{.value}")
    
    return(dat_metrics)
}

calculate_metrics_Fig2 <- function(dat_Fig2, total_filter = 0) {
    dat <- dat_Fig2 %>%
        filter(abs(total) > total_filter) %>%
        group_by(purpose, method_label_succinct, k) %>%
        summarize(slope = 1 - (t(total) %*% alpha / t(total) %*% total)[1],
                  estimate_Bias = abs(slope - mean(1 - p_direct)),
                  estimate_Variance = mean((MP - slope)^2),
                  estimate_Coverage = mean(between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
                  estimate_Power = mean(!between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
                  SE_Bias = summary(lm(alpha ~ 0 + total))$coef[1, 2],
                  SE_Variance = sqrt(var((MP - slope)^2) / n()), # For Gaussian, it would be SE_Variance = sqrt(2 * var(MP)^2 / (n()-1))
                  SE_Coverage = sqrt(var(between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP)) / n()),
                  SE_Power = sqrt(var(!between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP)) / n()),
                  .groups = "drop") %>%
        pivot_longer(cols = contains(c("Bias", "Variance", "Coverage", "Power")), 
                     names_to = c(".value", "metric"), 
                     names_pattern = "(.*)_(.*)") %>%
        mutate(CI_lower = estimate - 1.96*SE,
               CI_upper = estimate + 1.96*SE)
    
    return(dat)
}

stylize_gtab <- function(gtab, nrow) {
    gtab_styled <- gtab %>%
        fmt_percent(columns = -c(purpose, x), decimals = 1, drop_trailing_zeros = TRUE) %>%
        data_color(columns = starts_with(c("naive", "original", "integrated")), colors = col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000), domain = c(0, 1), na.color = "#f1f1f1")) %>%
        cols_width(starts_with(c("naive", "original", "integrated")) ~ px(63),
                   x ~ px(90),
                   purpose ~ px(40)) %>%
        cols_align(align = "center") %>%
        text_transform(locations = cells_body(columns = purpose),
                       fn = function(x) {
                           purrr::map(x, ~ gt::html(ifelse(.x == "default", "<span style='font-size:14pt; line-height:40px; overflow:visible'>default</span>",
                                                           paste0("<span style='font-size:14pt; line-height:10px; overflow:hidden;'>",
                                                                  case_when(.x == "nX" ~ "n<sub style='font-size:65%;'>X</sub>",
                                                                            .x == "nZ" ~ "n<sub style='font-size:65%;'>Z</sub>",
                                                                            .x == "nY" ~ "n<sub style='font-size:65%;'>Y</sub>",
                                                                            .x == "m" ~ "m",
                                                                            .x == "h2X" ~ "h<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>X</sub>",
                                                                            .x == "k" ~ "k",
                                                                            .x == "k-sig-perc" ~ "p<sub style='font-size:65%;'>k</sub>",
                                                                            .x == "p-value-mediators" ~ "P",
                                                                            .x == "cor-mediated" ~ "&rho;<sub style='font-size:65%;'>&gamma;,&delta;</sub>",
                                                                            .x == "var-explained-ZtoY" ~ "&sigma;<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>Y,Z</sub>",
                                                                            .x == "total-causal" ~ "&theta;",
                                                                            .x == "p-direct" ~ "MP",
                                                                            TRUE ~ .x),
                                                                  "</span>"))))
                       }) %>%
        text_transform(locations = cells_body(columns = x),
                       fn = function(x) {
                           purrr::map_chr(x, ~ ifelse(.x == "default", "<span style='font-size:16pt; line-height:63px; overflow:visible'>default</span>",
                                                      paste0("<span style='font-size:12pt;'>", xml2::xml_text(xml2::read_html(paste0("<x>", .x, "</x>"))), "</span>")))
                       }) %>%
        tab_style(locations = cells_body(columns = purpose), style = cell_text(v_align = "bottom", align = "left")) %>%
        tab_style(locations = cells_body(rows = c(TRUE, 1:nrow %% 4 == 0)), style = cell_borders(sides = c("bottom"), weight = px(2))) %>%
        tab_style(locations = cells_body(rows = x == "default"), style = list(cell_text(weight = "bold", size = px(14)),
                                                                              cell_borders(sides = "bottom", weight = px(7), style = "double"))) %>%
        tab_style(locations = cells_column_labels(), style = cell_borders(sides = c("top", "bottom"), weight = px(2))) %>%
        tab_options(table.border.top.style = "hidden",
                    table.font.names = "Helvetica",
                    table.font.size = 10,
                    table_body.border.bottom.color = "black",
                    column_labels.border.bottom.color = "black",
                    column_labels.font.size = 14,
                    data_row.padding = 2)
    
    return(gtab_styled)
}

parameter_mat <- function(sim_files) {
    parameters <- sim_files %>%
        sapply(function(x) {
            sub(".RDS", "", x) %>%
                strsplit("_")
        }) %>%
        bind_rows() %>%
        transpose() %>%
        set_rownames(sim_files) %>%
        set_colnames(c("method", "seed", "start", "end", "scenario", "purpose", "nX", "nZ", "nY", "m", "h2X", "k", "k_sig_perc", "p_value_mediators",
                       "cor_mediated", "var_explained_ZtoY", "total_causal", "p_direct", "Sigma_type", "sigmaC2", "sigmac2", "sigmab2"))
    parameters <- parameters %<>% 
        mutate_if(!(colnames(parameters) %in% c("method", "scenario", "purpose", "Sigma_type")), as.numeric)
    
    return(parameters)
}

sim_data <- function(params, sim_filename) {
    sim_params <- tibble(params) %>% # "tibble" eliminates the unwanted rowname
        select(scenario, purpose, m, nX, nZ, nY, h2X, total_causal, p_direct, cor_mediated, var_explained_ZtoY, k, k_sig_perc, p_value_mediators, Sigma_type, true_sigmaC2 = sigmaC2, true_sigmac2 = sigmac2, true_sigmab2 = sigmab2)
    
    sim_results <- readRDS(sim_filename) %>%
        lapply(function(sim) data.table(method = sim$method, generate_seed = sim$generate_seed, estimate_seed = sim$estimate_seed, true_total = sim$true_total, true_alpha = sim$true_alpha, total = sim$total, alpha = sim$alpha, 
                                        SE_total = sim$SE_total, SE_alpha = sim$SE_alpha, SE_total_naive = sim$SE_total_naive, SE_alpha_naive = sim$SE_alpha_naive, sigmag2 = sim$sigmag2, sigmad2 = sim$sigmad2, 
                                        sigmaC2 = sim$sigmaC2, sigmac2 = sim$sigmac2, sigmab2 = sim$sigmab2, k_selected = sim$k_selected, pleiotropy = sim$pleiotropy, convergence = sim$convergence, optim_value = sim$optim_value, 
                                        time = sim$time, sigmag2_convergence = sim$sigmag2_convergence, sigmad2_convergence = sim$sigmad2_convergence, sigmag2_bound = sim$sigmag2_bound, sigmad2_bound = sim$sigmad2_bound, is_seed_OK = sim$is_seed_OK)) %>%
        rbindlist(fill = TRUE)
    
    sims <- cbind(sim_params, sim_results) %>%
        relocate(method, pleiotropy, generate_seed, estimate_seed, .after = purpose)
    
    return(sims)
}

combine_simulation_data <- function(sim_files) {
    parameters <- parameter_mat(sub("^.*/", "", sim_files))
    
    simulation_data <- lapply(1:nrow(parameters), function(i) {
        print(i)
        params <- parameters[i, ]
        sim_filename <- sim_files[i]
        
        sim_data(params, sim_filename)
    }) %>%
        rbindlist(fill = TRUE)
    
    return(simulation_data)
}


simulations_wo_T1E <- list.files(paste0("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/results/urblauna/general/", c("Colaus", "Colaus-pleiotropy", "identity", "identity-pleiotropy", "randomExtreme", "randomExtreme-pleiotropy")), pattern = "*.RDS", full.names = TRUE) %>%
    combine_simulation_data() %>%
    filter(!(scenario %in% c("Colaus-pleiotropy", "identity") & purpose == "T1E"))
simulations_w_T1E <- list.files(paste0("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/results/", c("T1E/Colaus-pleiotropy", "T1E/identity", "Burgess/Colaus-pleiotropy", "Burgess/identity")), pattern = "*.RDS", full.names = TRUE) %>%
    combine_simulation_data()

dat_all <- rbind(simulations_wo_T1E, simulations_w_T1E)

dat <- distinct(dat_all) %>% # Everything should already be distinct but just in case
    mutate(is_cor_OK = ((total - alpha) / k_selected)^2 <= sigmag2 * sigmad2,
           MP = (total - alpha) / total,
           SE_MP = alpha^2/total^2 * (SE_alpha^2/alpha^2 + SE_total^2/total^2),
           P_MP = 2 * pnorm(-abs(MP / SE_MP)),
           SE_MP_naive = alpha^2/total^2 * (SE_alpha_naive^2/alpha^2 + SE_total_naive^2/total^2),
           P_MP_naive = 2 * pnorm(-abs(MP / SE_MP_naive)),
           method_label_succinct = factor(case_when(method == "naive" ~ "MR+MVMR",
                                                    method == "naive_zero" ~ "MR+MVMR (zeroed out)",
                                                    method == "Burgess" ~ "Zhu et al. 2022",
                                                    method == "original" ~ "LiMA", #"Likelihood method",
                                                    method == "original_diagonal" ~ "LiMA (identity)", #"Likelihood method",
                                                    method == "integrated_free" ~ "I-LiMA (free)", #"Integrated likelihood",
                                                    method == "integrated_just_alpha" ~ "I-LiMA (direct)", #"Integrated likelihood",
                                                    method == "integrated_fixed_naive" ~ "I-LiMA (naive)", #"Integrated likelihood",
                                                    method == "integrated_fixed_contained" ~ "I-LiMA"), #"Integrated likelihood",
                                          levels = c("MR+MVMR", "MR+MVMR (zeroed out)", "Zhu et al. 2022", "LiMA", "LiMA (identity)", "I-LiMA", "I-LiMA (direct)", "I-LiMA (naive)", "I-LiMA (free)"))) %>%
    #mutate(SE_MP = ifelse(is.na(SE_MP_naive), SE_MP, SE_MP_naive),
    #       P_MP = ifelse(is.na(P_MP_naive), P_MP, P_MP_naive)) %>%
    rbind(filter(., method %in% c("naive", "naive_zero", "Burgess")) %>% mutate(pleiotropy = TRUE)) %>%
    arrange(scenario, purpose, pleiotropy, generate_seed, estimate_seed, m, nX, nZ, nY, h2X, total_causal, p_direct, cor_mediated, var_explained_ZtoY, k, k_sig_perc, p_value_mediators, Sigma_type, true_sigmaC2, true_sigmac2, true_sigmab2, method)

dat_nZ <- mutate(dat[purpose == "nZ",], purpose_plotmath = "n[M]", x = nZ, default = (x == 10000), range = (x %in% c(1000, 3000, 10000, 30000)))
dat_nX <- mutate(dat[purpose == "nX",], purpose_plotmath = "n[X]", x = nX, default = (x == 300000), range = (x %in% c(5000, 10000, 50000, 300000)))
dat_nY <- mutate(dat[purpose == "nY",], purpose_plotmath = "n[Y]", x = nY, default = (x == 300000), range = (x %in% c(5000, 10000, 50000, 300000)))
dat_h2X <- mutate(dat[purpose == "h2X",], purpose_plotmath = "h[X]^{2}", x = h2X, default = (x == 0.35), range = (x %in% c(0.1, 0.2, 0.35, 0.5)))
dat_m <- mutate(dat[purpose == "m",], purpose_plotmath = "m", x = m, default = (x == 100), range = (x %in% c(100, 200, 400, 800)))
dat_p_value_mediators <- mutate(dat[purpose == "p-value-mediators",], purpose_plotmath = "P", x = p_value_mediators, default = (x == 0.0001), range = (x %in% c(0.00000005, 0.0001, 0.05, 0.2)))
dat_cor_mediated <- mutate(dat[purpose == "cor-mediated",], purpose_plotmath = "rho[gamma*list(,)*delta]^{2}", x = cor_mediated, default = (x == 0.15), range = (x %in% c(0.01, 0.15, 0.3, 0.5)))
dat_var_explained_ZtoY <- mutate(dat[purpose == "var-explained-ZtoY",], purpose_plotmath = "sigma[Y*list(,)*M]^{2}", x = var_explained_ZtoY, default = (x == 0.0192), range = (x %in% c(0.01, 0.0192, 0.04, 0.08)))
dat_k <- mutate(dat[purpose == "k",], purpose_plotmath = "k", x = k, default = (x == 500), range = (x %in% c(100, 300, 500, 1000)))
dat_k_sig_perc <- mutate(dat[purpose == "k-sig-perc",], purpose_plotmath = "p[k]", x = k_sig_perc, default = (x == 0.05), range = (x %in% c(0.01, 0.03, 0.05, 0.1)))
dat_total_causal <- mutate(dat[purpose == "total-causal",], purpose_plotmath = "E(theta)", x = total_causal, default = (x == 0.15), range = (x %in% c(0.05, 0.1, 0.15, 0.2)))
dat_p_indirect <- mutate(dat[purpose == "p-direct",], purpose_plotmath = "MP", x = 1 - p_direct, default = (p_direct == 0.85), range = (p_direct %in% c(0.8, 0.85, 0.9, 0.95)))
dat_T1E <- mutate(dat[purpose == "T1E",], purpose_plotmath = "T1E", x = k, default = (x == 1), range = (x %in% c(1, 5, 10, 25)))
dat_Truth <- mutate(dat[purpose == "Truth",], purpose_plotmath = "Tru", x = k, default = (x == 10), range = (x %in% c(1, 5, 10, 25)))

dat_not_filtered <- rbind(dat_nX, dat_nZ, dat_nY, dat_m, dat_h2X, dat_k, dat_k_sig_perc, dat_p_value_mediators, dat_cor_mediated, dat_var_explained_ZtoY, dat_total_causal, dat_p_indirect, dat_T1E, dat_Truth) %>%
    relocate(purpose_plotmath, x, default, range, .after = purpose) %>%
    mutate(purpose = factor(purpose, levels = unique(purpose)),
           purpose_plotmath = factor(purpose_plotmath, levels = unique(purpose_plotmath))) %>%
    arrange(x) %>%
    mutate(x = case_when(default & purpose == "p-value-mediators" ~ paste0("<b>", format(x, scientific = TRUE, digits = 1), "</b>"),
                         default & purpose != "p-value-mediators" ~ paste0("<b>", format(x, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "</b>"), 
                         TRUE ~ as.character(x)),
           x = factor(x, levels = unique(x)))

# We'll filter out results where the optimization didn't converge
dat_filtered <- filter(dat_not_filtered, convergence == 0 | method %in% c("naive", "naive_zero", "Burgess"),
                       sigmag2_convergence %in% c(0, NA),
                       sigmad2_convergence %in% c(0, NA))

##### Main figures: biases using the regression approach
dat_regression <- group_by(dat_filtered, scenario, purpose, purpose_plotmath, x, default, range, method, pleiotropy, m, nX, nZ, nY, h2X, total_causal, p_direct, cor_mediated, 
                           var_explained_ZtoY, k, k_sig_perc, p_value_mediators, Sigma_type, true_sigmaC2, true_sigmac2, true_sigmab2, method_label_succinct) %>%
    filter(row_number() <= 300) %>%
    summarize(SS_SE_total = sum(SE_total[!is.nan(SE_total) & !is.na(SE_total)]**2),
              SS_total = sum(total[!is.nan(SE_total) & !is.na(SE_total)]**2),
              lambda = sqrt(case_when(SS_total > SS_SE_total ~ 1 - SS_SE_total / SS_total, # Correction factor for regression dilution bias
                                      TRUE ~ 1)),
              slope = 1 - (t(total) %*% alpha / t(total) %*% total)[1],
              slope_corrected = 1 - (t(total) %*% alpha / t(total) %*% total)[1] / lambda,
              SE = summary(lm(alpha ~ 0 + total))$coef[1, 2],
              CI_lower = slope - 1.96 * SE,
              CI_upper = slope + 1.96 * SE,
              bias = slope - mean(1 - p_direct),
              MSE = mean((MP - (1 - p_direct))^2),
              Var = mean((MP - slope)^2),
              coverage = mean(between(1 - p_direct, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
              power = mean(!between(0, MP - 1.96*SE_MP, MP + 1.96*SE_MP)),
              total = mean(total),
              direct = mean(alpha),
              k_selected = round(mean(k_selected), 1),
              N = n(), .groups = "drop")

##### Figure 1 is the main figure about the simulations, showing the key parameters: bias, selected mediators, computation time etc
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 16
AXIS_TEXT_SIZE <- 14

methods <- c("integrated_fixed_contained", "original", "naive")
scenario <- "Colaus-pleiotropy"
pleiotropy <- TRUE
colors = c("#9C9C9C", "#E6A83A", "#D45626")

dat_Fig1 <- dat_regression %>%
    filter(method %in% methods, scenario == !!scenario, pleiotropy == !!pleiotropy) %>%
    arrange(method_label_succinct)

dat_time <- dat_filtered %>%
    filter(method %in% methods, scenario == !!scenario, pleiotropy == !!pleiotropy,
           m < 800, k_selected %in% c(1, 15, 30, 45)) %>%
    summarize(k_selected = unique(k_selected),
              N = n(),
              time_0 = min(time),
              time_25 = quantile(time, 0.25),
              time_50 = quantile(time, 0.5),
              time_75 = quantile(time, 0.75),
              time_100 = max(time),
              .by = c(method_label_succinct, k_selected))

gg_Fig1(dat_Fig1, dat_time, colors)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/Fig1_reordered.pdf", width = 15, height = 13)


##### Figure 2 is another main figure about the simulations, showing the key metrics: bias, variance, coverage, power, T1E
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 16
AXIS_TEXT_SIZE <- 14

methods <- c("integrated_fixed_contained", "original", "naive")
scenario <- "Colaus-pleiotropy"
pleiotropy <- TRUE
colors <- c("#9C9C9C", "#E6A83A", "#D45626")

dat_Fig2 <- dat_filtered %>%
    filter(method %in% methods, 
           scenario == !!scenario, 
           pleiotropy == !!pleiotropy,
           purpose %in% c("T1E", "Truth"), 
           k %in% c(1, 5, 10, 25))
dat_metrics_Fig2 <- calculate_metrics_Fig2(dat_Fig2)
dat_metrics_filtered_Fig2 <- calculate_metrics_Fig2(dat_Fig2, total_filter = 0.04)

gg_Fig2(dat_Fig2, dat_metrics_Fig2, dat_metrics_filtered_Fig2, colors)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/Fig2.pdf", width = 18, height = 15)


##### Supplementary figure about Sigma
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

scenarios <- c("Colaus-pleiotropy", "identity-pleiotropy", "randomExtreme-pleiotropy")
methods <- c("integrated_fixed_contained", "original", "naive")
colors <- c("#9C9C9C", "#E6A83A", "#D45626")

dat_sup <- dat_regression %>%
    filter(scenario %in% scenarios, method %in% methods, pleiotropy == TRUE, range == TRUE,
           !(purpose %in% c("T1E", "Truth"))) %>%
    mutate(scenario = case_when(scenario == "identity-pleiotropy" ~ "Identity",
                                scenario == "Colaus-pleiotropy" ~ "CoLaus",
                                scenario == "randomExtreme-pleiotropy" ~ "Extreme"),
           scenario = factor(scenario, levels = c("Identity", "CoLaus", "Extreme")))
highlight_scenario <- "CoLaus"
gg_sup(dat_sup, highlight_scenario, 0.32, colors)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/Sigma.pdf", width = 10, height = 13)


##### Supplementary figure about pleiotropy
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

scenarios <- c("Colaus", "Colaus-pleiotropy")
methods <- c("integrated_fixed_contained", "original", "naive")
colors <- c("#9C9C9C", "#E6A83A", "#D45626")

dat_sup <- filter(dat_regression, scenario %in% scenarios, method %in% methods, range == TRUE,
                  !(purpose %in% c("T1E", "Truth"))) %>%
    mutate(scenario = case_when(scenario == "Colaus" & pleiotropy ~ "generated: NO, estimated: YES",
                                scenario == "Colaus-pleiotropy" & pleiotropy ~ "generated: YES, estimated: YES",
                                scenario == "Colaus" & !pleiotropy ~ "generated: NO, estimated: NO",
                                scenario == "Colaus-pleiotropy" & !pleiotropy ~ "generated: YES, estimated: NO"),
           scenario = factor(scenario, levels = c("generated: NO, estimated: YES",
                                                  "generated: YES, estimated: YES",
                                                  "generated: NO, estimated: NO",
                                                  "generated: YES, estimated: NO")))
highlight_scenario <- "generated: YES, estimated: YES"
gg_sup(dat_sup, highlight_scenario, 0.46, colors)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/pleiotropy.pdf", width = 13, height = 14)


##### Supplementary figure about the zeroed out MR framework
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

methods <- c("naive", "naive_zero")
colors <- c("#9C9C9C", "#666666")

dat_sup <- filter(dat_regression, scenario == "Colaus-pleiotropy", method %in% methods, pleiotropy == TRUE, range == TRUE,
                  !(purpose %in% c("T1E", "Truth")))
highlight_scenario <- ""
gg_sup(dat_sup, highlight_scenario, 0.32, colors, TRUE)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/naive.pdf", width = 8, height = 8)


##### Supplementary figure about the original methods
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

methods <- c("original", "original_diagonal")
colors <- c("#E6A83A", "#986813")

dat_sup <- filter(dat_regression, scenario == "Colaus-pleiotropy", method %in% methods, pleiotropy == TRUE, range == TRUE,
                  !(purpose %in% c("T1E", "Truth")))
highlight_scenario <- ""
gg_sup(dat_sup, highlight_scenario, 0.32, colors, TRUE)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/original.pdf", width = 8, height = 8)


##### Supplementary figure about the integrated methods
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

methods <- c("integrated_fixed_contained", "integrated_just_alpha", "integrated_fixed_naive", "integrated_free")
colors <- c("#D45626", "#a6431d", "#7c3216", "#f8e1d8")

dat_sup <- filter(dat_regression, scenario == "Colaus-pleiotropy", method %in% methods, pleiotropy == TRUE, range == TRUE,
                  !(purpose %in% c("T1E", "Truth")))
highlight_scenario <- ""
gg_sup(dat_sup, highlight_scenario, 0.46, colors, TRUE)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/integrated.pdf", width = 12, height = 10)


##### Supplementary figures with Zhu et al. 2022
FONT_FAMILY <- "Helvetica"
PLOT_TAG_SIZE <- 20
AXIS_TITLE_SIZE <- 14
AXIS_TEXT_SIZE <- 12

methods <- c("integrated_fixed_contained", "Burgess", "naive")
pleiotropy <- TRUE
colors = c("#9C9C9C", "dodgerblue3", "#D45626")

gg_Zhu(dat_regression, dat_filtered, methods, "Colaus-pleiotropy", pleiotropy, colors, 0.21)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/Zhu_Colaus_pleiotropy.pdf", width = 14, height = 10)
gg_Zhu(dat_regression, dat_filtered, methods, "identity", pleiotropy, colors, 0.18)
#ggsave("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper/Zhu_identity_nopleiotropy.pdf", width = 14, height = 10)


##### Data for the supplementary tables
methods <- c("naive", "original", "integrated_fixed_contained")

create_dat_metrics <- function(dat_filtered, param = NULL, total_filter = 0) {
    dat_supp <- dat_filtered %>%
        filter(method %in% methods, 
               scenario == "Colaus-pleiotropy", 
               pleiotropy == TRUE, 
               range == TRUE,
               !(purpose %in% c("T1E", "Truth")),
               abs(total) > total_filter) %>%
        mutate(method = factor(method, levels = methods),
               x = factor(x, levels = rev(levels(x))))
    
    dat_metrics_default <- filter(dat_supp, default) %>%
        group_by(method) %>%
        calculate_metrics_sup() %>%
        mutate(purpose = "", x = "default", .before = 1)
    
    dat_metrics_all <- dat_supp %>%
        group_by(purpose, x, method) %>%
        calculate_metrics_sup() %>%
        arrange(purpose, x)
    
    dat_metrics <- rbind(dat_metrics_default, dat_metrics_all) %>%
        mutate(purpose = ifelse(duplicated(purpose), "", as.character(purpose)))
    
    if (is.null(param))
        dat_metrics <- select(dat_metrics, purpose, x, ends_with(c("converged", "bound")))
    else {
        dat_metrics <- select(dat_metrics, purpose, x, contains(param)) %>%
            select(!ends_with(paste0("bias_", param))) %>%
            rename_with(~ sub(paste0("_", param), "", .x), contains(param))
    }
    
    return(dat_metrics)
}

stylize_gtab <- function(gtab, nrow) {
    gtab_styled <- gtab %>%
        fmt_percent(columns = contains(c("coverage", "power", "prc", "converged", "bound")), decimals = 1, drop_trailing_zeros = TRUE) %>%
        fmt_number(columns = contains("variance"), n_sigfig = 1, drop_trailing_zeros = TRUE) %>%
        cols_width(starts_with(c("naive", "original", "integrated")) ~ px(63),
                   x ~ px(90),
                   purpose ~ px(40)) %>%
        cols_align(align = "center") %>%
        text_transform(locations = cells_body(columns = purpose),
                       fn = function(x) {
                           purrr::map(x, ~ gt::html(ifelse(.x == "default", "<span style='font-size:14pt; line-height:40px; overflow:visible'>default</span>",
                                                           paste0("<span style='font-size:14pt; line-height:10px; overflow:hidden;'>",
                                                                  case_when(.x == "nX" ~ "n<sub style='font-size:65%;'>X</sub>",
                                                                            .x == "nZ" ~ "n<sub style='font-size:65%;'>M</sub>",
                                                                            .x == "nY" ~ "n<sub style='font-size:65%;'>Y</sub>",
                                                                            .x == "m" ~ "m",
                                                                            .x == "h2X" ~ "h<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>X</sub>",
                                                                            .x == "k" ~ "k",
                                                                            .x == "k-sig-perc" ~ "p<sub style='font-size:65%;'>k</sub>",
                                                                            .x == "p-value-mediators" ~ "P",
                                                                            .x == "cor-mediated" ~ "&rho;<sub style='font-size:65%;'>&gamma;,&delta;</sub>",
                                                                            .x == "var-explained-ZtoY" ~ "&sigma;<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>Y,M</sub>",
                                                                            .x == "total-causal" ~ "E(&theta;)",
                                                                            .x == "p-direct" ~ "MP",
                                                                            TRUE ~ .x),
                                                                  "</span>"))))
                       }) %>%
        text_transform(locations = cells_body(columns = x),
                       fn = function(x) {
                           purrr::map_chr(x, ~ ifelse(.x == "default", "<span style='font-size:16pt; line-height:63px; overflow:visible'>default</span>",
                                                      paste0("<span style='font-size:12pt;'>", xml2::xml_text(xml2::read_html(paste0("<x>", .x, "</x>"))), "</span>")))
                       }) %>%
        tab_style(locations = cells_body(columns = purpose), style = cell_text(v_align = "bottom", align = "left")) %>%
        tab_style(locations = cells_body(rows = c(TRUE, 1:nrow %% 4 == 0)), style = cell_borders(sides = c("bottom"), weight = px(2))) %>%
        tab_style(locations = cells_body(rows = x == "default"), style = list(cell_text(weight = "bold", size = px(14)),
                                                                              cell_borders(sides = "bottom", weight = px(7), style = "double"))) %>%
        tab_style(locations = cells_column_labels(), style = cell_borders(sides = c("top", "bottom"), weight = px(2))) %>%
        tab_options(table.border.top.style = "hidden",
                    table.font.names = "Helvetica",
                    table.font.size = 10,
                    table_body.border.bottom.color = "black",
                    column_labels.border.bottom.color = "black",
                    column_labels.font.size = 14,
                    data_row.padding = 2)
    
    return(gtab_styled)
}

create_metrics_table <- function(dat_metrics, type = c("alpha", "total", "MP"), variance_level = 0.1) {
    type <- match.arg(type)
    
    gtab_metrics <- dat_metrics %>%
        gt() %>%
        cols_label(purpose = "",
                   x = "",
                   starts_with("naive") ~ "MR fmw",
                   starts_with("original") ~ "LiMA",
                   starts_with("integrated") ~ "I-LiMA") %>%
        tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>Bias</span>"), columns = contains("bias")) %>%
        tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>Variance</span>"), columns = contains("variance")) %>%
        tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>Coverage</span>"), columns = contains("coverage")) %>%
        tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>Power</span>"), columns = contains("power")) %>%
        tab_style(locations = cells_body(), style = cell_borders(weight = 0)) %>%
        tab_style(locations = cells_body(columns = c(x, starts_with("integrated"))), style = cell_borders(sides = "right", weight = px(2))) %>%
        tab_style(locations = cells_body(columns = contains("integrated") & contains(c("prc", "variance", "coverage"))), style = cell_borders(sides = "right", weight = px(7), style = "double")) %>%
        data_color(columns = contains(c("prc", "coverage", "power")), fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000), domain = c(0, 1), na.color = "#f1f1f1")) %>%
        data_color(columns = contains("variance"), fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000), domain = c(0, variance_level), na.color = "#f1f1f1"))
    
    for (column in colnames(dat_metrics)[grepl("bias_prc", colnames(dat_metrics))]) {
        if (any(dat_metrics[column] < 0))
            gtab_metrics <- data_color(gtab_metrics, columns = column, rows = dat_metrics[column] < 0, fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Blues")))(1000), reverse = TRUE, domain = c(-1, 0), na.color = "#f1f1f1"))
        
        if (any(dat_metrics[column] <= -1))
            gtab_metrics <- data_color(gtab_metrics, columns = column, rows = dat_metrics[column] <= -1, fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Blues")))(1000)[1000], reverse = TRUE, domain = c(-1e6, -1), na.color = "#f1f1f1"))
        
        if (any(dat_metrics[column] >= 1))
            gtab_metrics <- data_color(gtab_metrics, columns = column, rows = dat_metrics[column] >= 1, fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000)[1000], domain = c(1, 1e6), na.color = "#f1f1f1"))
    }
    
    for (column in colnames(dat_metrics)[grepl("variance", colnames(dat_metrics))]) {
        if (any(dat_metrics[column] >= variance_level))
            gtab_metrics <- data_color(gtab_metrics, columns = column, rows = dat_metrics[column] >= variance_level, fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000)[1000], domain = c(variance_level, 1e6), na.color = "#f1f1f1"))
    }
    
    gtab_metrics <- stylize_gtab(gtab_metrics, nrow = nrow(dat_metrics)-1)
    
    print(gtab_metrics)
    
    return(gtab_metrics)
}


##### Supplementary table about the metrics (bias, variance, coverage, power)
dat_metrics <- create_dat_metrics(dat_filtered, "MP")
dat_metrics_filtered <- create_dat_metrics(dat_filtered, "MP", total_filter = 0.04)
dat_metrics[, grepl("variance", colnames(dat_metrics))] <- dat_metrics_filtered[, grepl("variance", colnames(dat_metrics_filtered))]
gtab_metrics <- create_metrics_table(dat_metrics, variance_level = 0.5)
#gtsave(gtab_metrics, filename = "metrics_MP.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper", vwidth = 3822, vheight = 2846, zoom = 3)

dat_metrics <- create_dat_metrics(dat_filtered, "total")
dat_metrics_filtered <- create_dat_metrics(dat_filtered, "total", total_filter = 0.04)
dat_metrics[, grepl("variance", colnames(dat_metrics)) | grepl("bias", colnames(dat_metrics))] <- dat_metrics_filtered[, grepl("variance", colnames(dat_metrics_filtered)) | grepl("bias", colnames(dat_metrics_filtered))]
gtab_metrics <- create_metrics_table(dat_metrics, variance_level = 0.001)
#gtsave(gtab_metrics, filename = "metrics_total.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper", vwidth = 3822, vheight = 2846, zoom = 3)

dat_metrics <- create_dat_metrics(dat_filtered, "alpha")
dat_metrics_filtered <- create_dat_metrics(dat_filtered, "alpha", total_filter = 0.04)
dat_metrics[, grepl("variance", colnames(dat_metrics)) | grepl("bias", colnames(dat_metrics))] <- dat_metrics_filtered[, grepl("variance", colnames(dat_metrics_filtered)) | grepl("bias", colnames(dat_metrics_filtered))]
gtab_metrics <- create_metrics_table(dat_metrics, variance_level = 0.005)
#gtsave(gtab_metrics, filename = "metrics_alpha.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper", vwidth = 3822, vheight = 2846, zoom = 3)


##### Supplementary table about the failures
gtab_failures <- create_dat_metrics(dat_not_filtered) %>%
    select(purpose, x, original_not_converged, contains("integrated") & ends_with(c("bound", "converged"))) %>%
    gt() %>%
    cols_label(purpose = "",
               x = "",
               ends_with("not_converged") ~ "f",
               contains("sigmag2") ~ gt::html("&sigma;<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>&gamma;</sub>"),
               contains("sigmad2") ~ gt::html("&sigma;<sup style='font-size:65%;'>2</sup><sub style='position:relative; font-size:65%; left:-.5em;'>&delta;</sub>")) %>%
    tab_spanner(label = gt::html("I-LiMA"), id = "integrated_converged", columns = contains("integrated") & contains("converged")) %>%
    tab_spanner(label = gt::html("I-LiMA"), id = "integrated_bound", columns = contains("integrated") & contains("bound")) %>%
    tab_spanner(label = gt::html("LiMA"), columns = original_not_converged) %>%
    tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>Not converged</span>"), columns = contains("converged")) %>%
    tab_spanner(label = gt::html("<span style='font-size:14pt; font-weight:bold;'>On the bound</span>"), columns = contains("bound")) %>%
    tab_style(locations = cells_body(), style = cell_borders(weight = 0)) %>%
    tab_style(locations = cells_body(columns = c(x, original_not_converged)), style = cell_borders(sides = "right", weight = px(2))) %>%
    tab_style(locations = cells_body(columns = contains("integrated") & contains("sigmad2_not_converged")), style = cell_borders(sides = "right", weight = px(2))) %>%
    tab_style(locations = cells_body(columns = contains("integrated") & contains("sigmad2_not_converged")), style = cell_borders(sides = "right", weight = px(7), style = "double")) %>%
    data_color(columns = contains(c("converged", "bound")), fn = scales::col_numeric(colorRampPalette(c("white", brewer.pal(n = 3, name = "Reds")))(1000), domain = c(0, 1), na.color = "#f1f1f1")) %>%
    stylize_gtab(nrow = nrow(dat_metrics)-1)
print(gtab_failures)
#gtsave(gtab_failures, filename = "failures.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper", vwidth = 1554, vheight = 3855, zoom = 3)


##### Supplementary table about the distribution of mediator QTLs
FONT_FAMILY <- "Helvetica"

dat_eqtls <- fread("/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/data/mediator_QTLs.txt") %>%
    select(`Nr of eQTLs` = n_eqtls, Probability = p, `Cumulative probability` = cumulative_p)

gtab <- gt(dat_eqtls) %>%
    fmt_number(columns = c(Probability, `Cumulative probability`), decimals = 3, drop_trailing_zeros = TRUE) %>%
    cols_width(Probability ~ px(120)) %>%
    tab_style(locations = cells_body(columns = everything()), style = cell_borders(sides = c("top"), weight = px(0))) %>%
    tab_style(locations = cells_body(columns = everything()), style = cell_text(align = "center")) %>%
    tab_style(locations = cells_column_labels(), style = list(cell_text(align = "center", weight = "bold"), cell_borders(sides = c("top", "bottom"), color = "black", weight = px(2)))) %>%
    tab_style(locations = cells_body(rows = nrow(dat_eqtls)), style = cell_borders(sides = c("bottom"), weight = px(2))) %>%
    tab_options(table.font.names = FONT_FAMILY,
                table.font.size = 12,
                column_labels.font.size = 14,
                data_row.padding = 6)
print(gtab)
#gtsave(gtab, filename = "mediator_QTLs.png", path = "/Users/kaidolepik/Desktop/Work/PROJECTS_CH/ML_mediated/simulations/figures/paper", zoom = 6)
