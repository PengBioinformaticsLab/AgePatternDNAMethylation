library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)


# paths

project_dir <- getwd()

allTRUE_path    <- file.path(project_dir, "data", "matrix", "allTRUE.csv")
maleTRUE_path   <- file.path(project_dir, "data", "matrix", "maleTRUE.csv")
femaleTRUE_path <- file.path(project_dir, "data", "matrix", "femaleTRUE.csv")

rds_dir <- file.path(project_dir, "relevant_rds")
out_dir <- project_dir


windows_use <- c(
  "0-20","5-25","10-30","15-35","20-40",
  "25-45","30-50","35-55","40-60","45-65",
  "50-70","55-75","60-80"
)

young_windows <- c("0-20","5-25","10-30","15-35","20-40")


# helper functions

get_corr_col <- function(results, w) {
  cand <- c(paste0(w, "-corr_res"), paste0(w, "-corr"))
  hit <- cand[cand %in% colnames(results)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_abs_col <- function(results, w) {
  cand <- c(paste0(w, "-abs_res"), paste0(w, "-abs_beta"), paste0(w, "-abs"))
  hit <- cand[cand %in% colnames(results)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_sd_col <- function(results, w) {
  cand <- c(paste0(w, "-sd_res0.5"), paste0(w, "-sd_res"))
  hit <- cand[cand %in% colnames(results)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

read_processed_matrix <- function(path) {
  x <- fread(path)
  cpg <- as.character(x[[1]])
  mat <- as.data.frame(x[, -1, with = FALSE])
  rownames(mat) <- cpg
  mat[] <- lapply(mat, function(z) suppressWarnings(as.numeric(as.character(z))))
  mat
}

compute_window_means <- function(mat, windows) {
  stopifnot(all(grepl("^age_[0-9]+$", colnames(mat))))
  out <- vector("list", length(windows))
  names(out) <- windows
  for (w in windows) {
    bounds <- strsplit(w, "-", fixed = TRUE)[[1]]
    a <- as.integer(bounds[1])
    b <- as.integer(bounds[2])
    cols <- intersect(paste0("age_", a:b), colnames(mat))
    if (length(cols) == 0) { warning("No columns found for window: ", w); next }
    out[[w]] <- rowMeans(mat[, cols, drop = FALSE], na.rm = TRUE)
  }
  out
}

make_long_sd <- function(results, sex_label, windows) {
  out <- list()
  for (w in windows) {
    sd_col <- get_sd_col(results, w)
    if (is.na(sd_col)) next
    vals <- as.numeric(results[[sd_col]])
    df <- data.frame(CpG = rownames(results), window = w, value = vals,
                     sex = sex_label, stringsAsFactors = FALSE)
    df <- df[is.finite(df$value), , drop = FALSE]
    out[[w]] <- df
  }
  bind_rows(out)
}

make_long_abs <- function(results, sex_label, windows) {
  out <- list()
  for (w in windows) {
    abs_col <- get_abs_col(results, w)
    if (is.na(abs_col)) next
    vals <- as.numeric(results[[abs_col]])
    df <- data.frame(CpG = rownames(results), window = w, value = vals,
                     sex = sex_label, stringsAsFactors = FALSE)
    df <- df[is.finite(df$value), , drop = FALSE]
    out[[w]] <- df
  }
  bind_rows(out)
}

make_long_ratio_from_matrix <- function(results, mean_list, sex_label, windows) {
  out <- list()
  for (w in windows) {
    sd_col <- get_sd_col(results, w)
    if (is.na(sd_col)) next
    if (is.null(mean_list[[w]])) next
    sdv <- as.numeric(results[[sd_col]])
    mn <- as.numeric(mean_list[[w]][rownames(results)])
    denom <- sqrt(mn * (1 - mn))
    ratio <- sdv / denom
    ratio[!is.finite(ratio)] <- NA
    df <- data.frame(CpG = rownames(results), window = w, value = ratio,
                     sex = sex_label, stringsAsFactors = FALSE)
    df <- df[is.finite(df$value), , drop = FALSE]
    out[[w]] <- df
  }
  bind_rows(out)
}

make_density_plot <- function(df, title_txt, xlab_txt, vline = NULL, xlim_max = NULL) {
  if (is.null(df) || nrow(df) == 0) stop(paste("Empty dataframe for:", title_txt))
  df$window <- factor(df$window, levels = windows_use)
  p <- ggplot(df, aes(x = value, color = sex, fill = sex)) +
    geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1) +
    facet_wrap(~ window, scales = "free_y", ncol = 3) +
    labs(title = title_txt, x = xlab_txt, y = "Density") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  if (!is.null(vline)) p <- p + geom_vline(xintercept = vline, linetype = "dashed", linewidth = 0.7)
  if (!is.null(xlim_max)) p <- p + coord_cartesian(xlim = c(0, xlim_max))
  p
}

make_tail_density_plot <- function(df, title_txt, xlab_txt, cutoff, x_max = NULL, y_max = NULL) {
  if (is.null(df) || nrow(df) == 0) stop(paste("Empty dataframe for:", title_txt))
  df$window <- factor(df$window, levels = windows_use)
  
  p <- ggplot(df, aes(x = value, color = sex, fill = sex)) +
    geom_density(alpha = 0.15, linewidth = 0.8, adjust = 1) +
    facet_wrap(~ window, scales = "free_y", ncol = 3) +
    labs(title = title_txt, x = xlab_txt, y = "Density") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    ) +
    geom_vline(xintercept = cutoff, linetype = "dashed", linewidth = 0.7)
  
  if (!is.null(x_max) && !is.null(y_max)) {
    p <- p + coord_cartesian(xlim = c(cutoff, x_max), ylim = c(0, y_max))
  } else if (!is.null(x_max)) {
    p <- p + coord_cartesian(xlim = c(cutoff, x_max))
  } else if (!is.null(y_max)) {
    p <- p + coord_cartesian(ylim = c(0, y_max))
  }
  
  p
}


# load retained CpGs

allTRUE_dt <- fread(allTRUE_path)
keep_cpgs  <- unique(as.character(allTRUE_dt[[1]]))
keep_cpgs  <- keep_cpgs[!is.na(keep_cpgs) & nzchar(keep_cpgs)]

cat("allTRUE CpGs retained:", length(keep_cpgs), "\n")


# load/match result tables to CpGs

all_results    <- readRDS(file.path(rds_dir, "all_results.rds"))
male_results   <- readRDS(file.path(rds_dir, "male_results.rds"))
female_results <- readRDS(file.path(rds_dir, "female_results.rds"))

all_results_match    <- all_results[intersect(rownames(all_results), keep_cpgs), , drop = FALSE]
male_results_match   <- male_results[intersect(rownames(male_results), keep_cpgs), , drop = FALSE]
female_results_match <- female_results[intersect(rownames(female_results), keep_cpgs), , drop = FALSE]

cat("all_results_match:", nrow(all_results_match), "\n")
cat("male_results_match:", nrow(male_results_match), "\n")
cat("female_results_match:", nrow(female_results_match), "\n\n")


# load mean matrices, restrict to allTRUE

male_mat   <- read_processed_matrix(maleTRUE_path)
female_mat <- read_processed_matrix(femaleTRUE_path)

male_mat   <- male_mat[intersect(rownames(male_mat), keep_cpgs), , drop = FALSE]
female_mat <- female_mat[intersect(rownames(female_mat), keep_cpgs), , drop = FALSE]

cat("nrow(male_mat):", nrow(male_mat), "\n")
cat("nrow(female_mat):", nrow(female_mat), "\n\n")

male_mean_list   <- compute_window_means(male_mat, windows_use)
female_mean_list <- compute_window_means(female_mat, windows_use)


# aaCpG filtering by window

count_forced_top100 <- function(results, label,
                                cor_cut = 0.5, abs_cut_young = 0.2,
                                abs_cut_old = 0.2, sd_cut = 2, top_n = 100) {
  out <- list()
  for (w in windows_use) {
    corr_col <- get_corr_col(results, w)
    abs_col  <- get_abs_col(results, w)
    sd_col   <- get_sd_col(results, w)
    if (is.na(corr_col) || is.na(abs_col) || is.na(sd_col)) next
    abs_cut <- if (w %in% young_windows) abs_cut_young else abs_cut_old
    df <- data.frame(
      CpG        = rownames(results),
      corr       = as.numeric(results[[corr_col]]),
      abs_change = as.numeric(results[[abs_col]]),
      sd_metric  = as.numeric(results[[sd_col]]),
      stringsAsFactors = FALSE
    )
    df <- df[is.finite(df$corr) & is.finite(df$abs_change) & is.finite(df$sd_metric), , drop = FALSE]
    if (nrow(df) == 0) next
    df$fail_corr <- abs(df$corr) < cor_cut
    df$fail_abs  <- df$abs_change < abs_cut
    df$fail_sd   <- df$sd_metric >= sd_cut
    df$pass_all  <- (!df$fail_corr) & (!df$fail_abs) & (!df$fail_sd)
    pass_df <- df[df$pass_all, , drop = FALSE]
    if (nrow(pass_df) >= top_n) {
      selected <- pass_df %>% arrange(desc(abs_change)) %>% slice(1:top_n)
    } else {
      need_n  <- top_n - nrow(pass_df)
      fill_df <- df[!df$pass_all, , drop = FALSE] %>% arrange(desc(abs_change)) %>% head(need_n)
      selected <- bind_rows(pass_df, fill_df) %>% arrange(desc(abs_change)) %>% slice(1:top_n)
    }
    out[[w]] <- data.frame(
      group = label, window = w,
      total_selected = nrow(selected),
      corr_fail = sum(selected$fail_corr, na.rm = TRUE),
      abs_fail  = sum(selected$fail_abs,  na.rm = TRUE),
      sd_fail   = sum(selected$fail_sd,   na.rm = TRUE),
      corr_pass = sum(!selected$fail_corr, na.rm = TRUE),
      abs_pass  = sum(!selected$fail_abs,  na.rm = TRUE),
      sd_pass   = sum(!selected$fail_sd,   na.rm = TRUE),
      pass_all  = sum(selected$pass_all,   na.rm = TRUE),
      forced_in = sum(!selected$pass_all,  na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
  bind_rows(out)
}

all_top100_match    <- count_forced_top100(all_results_match, "All")
male_top100_match   <- count_forced_top100(male_results_match, "Male")
female_top100_match <- count_forced_top100(female_results_match, "Female")

top100_counts <- bind_rows(all_top100_match, male_top100_match, female_top100_match)

cat("\naaCpG filtering summary for the final selected set\n")
print(top100_counts)


# aaCpG summary table

get_pass_sites_by_window <- function(results,
                                     cor_cut = 0.5,
                                     abs_cut_young = 0.2,
                                     abs_cut_old = 0.2,
                                     sd_cut = 2) {
  out <- list()
  
  for (w in windows_use) {
    corr_col <- get_corr_col(results, w)
    abs_col  <- get_abs_col(results, w)
    sd_col   <- get_sd_col(results, w)
    if (is.na(corr_col) || is.na(abs_col) || is.na(sd_col)) next
    
    abs_cut <- if (w %in% young_windows) abs_cut_young else abs_cut_old
    
    df <- data.frame(
      CpG        = rownames(results),
      corr       = as.numeric(results[[corr_col]]),
      abs_change = as.numeric(results[[abs_col]]),
      sd_metric  = as.numeric(results[[sd_col]]),
      stringsAsFactors = FALSE
    )
    
    df <- df[is.finite(df$corr) & is.finite(df$abs_change) & is.finite(df$sd_metric), , drop = FALSE]
    
    if (nrow(df) == 0) {
      out[[w]] <- character(0)
      next
    }
    
    df$pass_all <- abs(df$corr) >= cor_cut &
      df$abs_change >= abs_cut &
      df$sd_metric < sd_cut
    
    pass_df <- df[df$pass_all, , drop = FALSE] %>%
      arrange(desc(abs_change))
    
    out[[w]] <- pass_df$CpG
  }
  
  out
}

all_pass_list    <- get_pass_sites_by_window(all_results_match)
male_pass_list   <- get_pass_sites_by_window(male_results_match)
female_pass_list <- get_pass_sites_by_window(female_results_match)

aaCpG_final_table <- data.frame(
  `Age Window` = windows_use,
  All = sapply(windows_use, function(w) {
    n_pass <- length(unique(all_pass_list[[w]]))
    ifelse(n_pass >= 100, n_pass, 100)
  }),
  Male = sapply(windows_use, function(w) {
    n_pass <- length(unique(male_pass_list[[w]]))
    ifelse(n_pass >= 100, n_pass, 100)
  }),
  Female = sapply(windows_use, function(w) {
    n_pass <- length(unique(female_pass_list[[w]]))
    ifelse(n_pass >= 100, n_pass, 100)
  }),
  `Male Female Overlap` = sapply(windows_use, function(w) {
    length(intersect(male_pass_list[[w]], female_pass_list[[w]]))
  }),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

all_unique_sites    <- unique(unlist(all_pass_list, use.names = FALSE))
male_unique_sites   <- unique(unlist(male_pass_list, use.names = FALSE))
female_unique_sites <- unique(unlist(female_pass_list, use.names = FALSE))
male_female_overlap_unique <- intersect(male_unique_sites, female_unique_sites)

total_row <- data.frame(
  `Age Window` = "Total Unique Sites",
  All = length(all_unique_sites),
  Male = length(male_unique_sites),
  Female = length(female_unique_sites),
  `Male Female Overlap` = length(male_female_overlap_unique),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

aaCpG_final_table <- bind_rows(aaCpG_final_table, total_row)

cat("\naaCpG summary table by age window\n")
print(aaCpG_final_table)


# counts/percentages by each cutoff

cut_by_cutoff_table <- top100_counts %>%
  transmute(
    group, window, total_selected,
    corr_cut_n   = corr_fail,
    abs_cut_n    = abs_fail,
    sd_cut_n     = sd_fail,
    corr_cut_pct = round(100 * corr_fail / total_selected, 2),
    abs_cut_pct  = round(100 * abs_fail  / total_selected, 2),
    sd_cut_pct   = round(100 * sd_fail   / total_selected, 2),
    pass_all_n   = pass_all,
    forced_in_n  = forced_in
  ) %>%
  arrange(group, factor(window, levels = windows_use))

cat("\n# of selected CpGs flagged by each cutoff\n")
print(cut_by_cutoff_table)


# male vs female fail counts 

make_fail_sex_table <- function(male_results, female_results, windows,
                                metric = c("corr", "abs", "sd"),
                                corr_cut = 0.5, abs_cut = 0.2, sd_cut = 2) {
  metric <- match.arg(metric)
  out <- list()
  for (w in windows) {
    if (metric == "corr") {
      colname <- paste0(w, "-corr_res")
      male_vals   <- as.numeric(male_results[[colname]])
      female_vals <- as.numeric(female_results[[colname]])
      male_fail   <- abs(male_vals) < corr_cut
      female_fail <- abs(female_vals) < corr_cut
    }
    if (metric == "abs") {
      colname <- paste0(w, "-abs_res")
      male_vals   <- as.numeric(male_results[[colname]])
      female_vals <- as.numeric(female_results[[colname]])
      male_fail   <- male_vals < abs_cut
      female_fail <- female_vals < abs_cut
    }
    if (metric == "sd") {
      colname <- paste0(w, "-sd_res0.5")
      male_vals   <- as.numeric(male_results[[colname]])
      female_vals <- as.numeric(female_results[[colname]])
      male_fail   <- male_vals >= sd_cut
      female_fail <- female_vals >= sd_cut
    }
    male_valid   <- is.finite(male_vals)
    female_valid <- is.finite(female_vals)
    male_total   <- sum(male_valid)
    female_total <- sum(female_valid)
    fail_n_male   <- sum(male_fail[male_valid],    na.rm = TRUE)
    fail_n_female <- sum(female_fail[female_valid], na.rm = TRUE)
    out[[w]] <- data.frame(
      window          = w,
      fail_n_Male     = fail_n_male,
      fail_n_Female   = fail_n_female,
      fail_n_diff     = fail_n_male - fail_n_female,
      fail_pct_Male   = round(100 * fail_n_male   / male_total,   2),
      fail_pct_Female = round(100 * fail_n_female / female_total, 2),
      fail_pct_diff   = round(100 * fail_n_male / male_total - 100 * fail_n_female / female_total, 2),
      stringsAsFactors = FALSE
    )
  }
  as_tibble(bind_rows(out))
}

corr_fail_table <- make_fail_sex_table(male_results_match, female_results_match, windows_use, metric = "corr")
abs_fail_table  <- make_fail_sex_table(male_results_match, female_results_match, windows_use, metric = "abs")
sd_fail_table   <- make_fail_sex_table(male_results_match, female_results_match, windows_use, metric = "sd")

cat("\nmale vs female fail table: correlation cutoff\n")
print(corr_fail_table)

cat("\nmale vs female fail table: absolute change cutoff\n")
print(abs_fail_table)

cat("\nmale vs female fail table: SD cutoff\n")
print(sd_fail_table)


# build long-format data for density

male_sd_df    <- make_long_sd(male_results_match, "Male", windows_use)
female_sd_df  <- make_long_sd(female_results_match, "Female", windows_use)
sd_df         <- bind_rows(male_sd_df, female_sd_df)

male_abs_df   <- make_long_abs(male_results_match, "Male", windows_use)
female_abs_df <- make_long_abs(female_results_match, "Female", windows_use)
abs_df        <- bind_rows(male_abs_df, female_abs_df)

male_ratio_df   <- make_long_ratio_from_matrix(male_results_match, male_mean_list, "Male", windows_use)
female_ratio_df <- make_long_ratio_from_matrix(female_results_match, female_mean_list, "Female", windows_use)
ratio_df        <- bind_rows(male_ratio_df, female_ratio_df)

cat("\ndensity plot data sizes\n")
cat("nrow(sd_df):", nrow(sd_df), "\n")
cat("nrow(abs_df):", nrow(abs_df), "\n")
cat("nrow(ratio_df):", nrow(ratio_df), "\n")


# main density plots (full distrib)

p_sd <- make_density_plot(
  sd_df,
  "Male vs Female: SD(beta)",
  "SD(beta)",
  vline = 0.2,
  xlim_max = 0.35
)

p_abs <- make_density_plot(
  abs_df,
  "Male vs Female: |Δβ|",
  "|Δβ|",
  vline = 0.2,
  xlim_max = 0.35
)

p_ratio <- make_density_plot(
  ratio_df,
  "Male vs Female: SD(beta) / sqrt(mean(beta) * (1 - mean(beta)))",
  "SD(beta) / sqrt(mean(beta) * (1 - mean(beta)))",
  vline = 0.2,
  xlim_max = 1.0
)

dev.new(width = 12, height = 10)
print(p_sd)

dev.new(width = 12, height = 10)
print(p_abs)

dev.new(width = 12, height = 10)
print(p_ratio)


# tail zoom plots

p_sd_tail <- make_tail_density_plot(
  sd_df,
  "Male vs Female: SD(beta) Tail Zoom",
  "SD(beta)",
  cutoff = 0.3,
  x_max  = 0.6
)

p_ratio_tail <- make_tail_density_plot(
  ratio_df,
  "Male vs Female: SD(beta) / sqrt(mean(beta) * (1 - mean(beta))) Tail Zoom",
  "SD(beta) / sqrt(mean(beta) * (1 - mean(beta)))",
  cutoff = 0.75,
  x_max  = 2.0
)

p_abs_tail <- make_tail_density_plot(
  abs_df,
  "Male vs Female: |Δβ| Tail Zoom",
  "|Δβ|",
  cutoff = 0.2,
  x_max  = 0.5,
  y_max  = 2
)

dev.new(width = 12, height = 10)
print(p_sd_tail)

dev.new(width = 12, height = 10)
print(p_ratio_tail)

dev.new(width = 12, height = 10)
print(p_abs_tail)


# counts of CpGs in the tail per window per sex

abs_tail_counts <- abs_df %>%
  filter(is.finite(value)) %>%
  mutate(window = factor(as.character(window), levels = windows_use)) %>%
  group_by(sex, window) %>%
  summarise(total = n(), n_gt_0.2  = sum(value > 0.2,  na.rm = TRUE),
            pct_gt_0.2  = round(100 * n_gt_0.2  / total, 2), .groups = "drop") %>%
  arrange(window, sex)

sd_tail_counts <- sd_df %>%
  filter(is.finite(value)) %>%
  mutate(window = factor(as.character(window), levels = windows_use)) %>%
  group_by(sex, window) %>%
  summarise(total = n(), n_gt_0.3  = sum(value > 0.3,  na.rm = TRUE),
            pct_gt_0.3  = round(100 * n_gt_0.3  / total, 2), .groups = "drop") %>%
  arrange(window, sex)

ratio_tail_counts <- ratio_df %>%
  filter(is.finite(value)) %>%
  mutate(window = factor(as.character(window), levels = windows_use)) %>%
  group_by(sex, window) %>%
  summarise(total = n(), n_gt_0.75 = sum(value > 0.75, na.rm = TRUE),
            pct_gt_0.75 = round(100 * n_gt_0.75 / total, 2), .groups = "drop") %>%
  arrange(window, sex)

cat("\n|Δβ| > 0.2: counts per window per sex\n")
print(abs_tail_counts)

cat("\nSD(beta) > 0.3: counts per window per sex\n")
print(sd_tail_counts)

cat("\nRatio > 0.75: counts per window per sex\n")
print(ratio_tail_counts)

