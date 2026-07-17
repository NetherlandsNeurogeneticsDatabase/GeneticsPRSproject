# ==============================================================================
# 1. PREDEFINED VARIABLES & CONFIGURATION
# ==============================================================================

# --- Library Imports ---
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(broom)
  library(forcats)
  library(parallel)
})

# --- File Paths ---
INPUT_PATH  <- "/home/../coombes/coombes_df.csv"
OUTPUT_DIR  <- "output"

# --- Analysis Configurations ---
PRS_VARS <- c(
  'AD_APOE', 'MDD', 'BP', 'SCZ', 'risk_taking', 
  'intelligence', 'insomnia_schoeler_2023', 'sociability', 'neuroticism_hu'
)

ATTRIBUTES_OF_INTEREST <- c(
  "Depressed_mood", "Mania", "Paranoia_suspiciousness", "Hallucinations",
  "Suicidal_ideation", "Aggressive_behavior", "Delusions", "Delirium", 
  "Psychiatric_admissions"
)

# Named colors mapped directly to attributes
CUSTOM_COLORS <- setNames(
  c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"),
  ATTRIBUTES_OF_INTEREST
)

N_PERMUTATIONS <- 10000  # Set to lower (e.g., 50) for fast iterative testing

# ==============================================================================
# 2. REQUIRED DATA (INPUT)
# ==============================================================================

if (!file.exists(INPUT_PATH)) {
  stop(paste("CRITICAL ERROR: Input data file missing at", INPUT_PATH))
}

df <- read_csv(INPUT_PATH, show_col_types = FALSE)

# Build a unique baseline matrix of donors
donor_df <- df %>% 
  distinct(donor_id, .keep_all = TRUE) %>%
  select(donor_id, nd_simple, age, most_likely_sex, PC1, PC2, PC3, all_of(PRS_VARS))


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

build_df_att <- function(att_of_int, disorder = NULL, disorder_exclude = NULL, u = NULL) {
  # Generate boolean targets indicating symptom presence
  att_bool_df <- df %>%
    filter(
      (is.null(disorder) | nd_simple %in% disorder) &
      (is.null(disorder_exclude) | !nd_simple %in% disorder_exclude)
    ) %>%
    group_by(donor_id) %>%
    summarize(
      att_bool = if (is.null(u)) {
        as.integer(any(attribute %in% att_of_int))
      } else {
        as.integer(any(attribute %in% att_of_int & age_at_symptom < u, na.rm = TRUE))
      },
      .groups = "drop"
    )
  
  donor_df_filtered <- donor_df %>%
    filter(
      (is.null(disorder) | nd_simple %in% disorder) &
      (is.null(disorder_exclude) | !nd_simple %in% disorder_exclude)
    )
  
  df_att <- donor_df_filtered %>%
    left_join(att_bool_df, by = "donor_id")
  
  return(df_att)
}

run_prs_models <- function(df_att, label) {
  covars <- c("PC1", "PC2", "PC3")
  
  results <- lapply(PRS_VARS, function(prs) {
    formula_obj <- as.formula(
      paste("att_bool ~", prs, "+", paste(covars, collapse = " + "))
    )
    
    model <- glm(formula_obj, data = df_att, family = binomial)
    
    tidy(model) %>%
      filter(term == prs) %>%
      mutate(
        PRS = prs,
        CI_lower = estimate - 1.96 * std.error,
        CI_upper = estimate + 1.96 * std.error,
        Selection = label
      ) %>%
      select(PRS, estimate, CI_lower, CI_upper, p.value, Selection)
  }) %>%
    bind_rows()
  
  return(results)
}

run_permutation_test <- function(df_att, prs_vars, n_perms = 100) {
  # Multi-threaded cluster setup
  n_cores <- max(1, parallel::detectCores() - 1)
  covars_rhs <- "+ PC1 + PC2 + PC3"
  
  perm_min_pvalues <- unlist(
    parallel::mclapply(1:n_perms, mc.cores = n_cores, FUN = function(perm) {
      df_perm <- df_att
      df_perm$att_bool <- sample(df_perm$att_bool)
      
      perm_pvalues <- sapply(prs_vars, function(prs) {
        formula_str <- paste("att_bool ~", prs, covars_rhs)
        mod_perm <- glm(as.formula(formula_str), data = df_perm, family = binomial())
        coef_table <- summary(mod_perm)$coefficients
        
        if (prs %in% rownames(coef_table)) {
          coef_table[prs, 4] 
        } else {
          NA
        }
      })
      min(perm_pvalues, na.rm = TRUE)
    })
  )
  return(perm_min_pvalues)
}


# ==============================================================================
# 4. PROCESSES
# ==============================================================================

# Configure target clinical subsets for profiling
disorder_sets <- list(
  PSYCH = list(disorder = c("MDD", "Other PSYCH", "BP"), exclude = NULL)
)

all_results <- list()

for (dname in names(disorder_sets)) {
  cat("\nRunning analysis pipeline for disorder cohort group:", dname, "\n")
  
  disorder1 <- disorder_sets[[dname]]$disorder
  disorder_exclude1 <- disorder_sets[[dname]]$exclude
  
  results_list <- list()
  perm_thresholds <- c()
  attr_counts <- c()
  n_total_disorder <- NA
  
  for (att in ATTRIBUTES_OF_INTEREST) {
    df_att <- build_df_att(
      att_of_int = att,
      disorder = disorder1,
      disorder_exclude = disorder_exclude1,
      u = NULL
    )
    
    n_total <- length(unique(df_att$donor_id))
    n_attr  <- sum(df_att$att_bool == 1, na.rm = TRUE)
    
    if (is.na(n_total_disorder)) {
      n_total_disorder <- n_total
    }
    
    attr_counts[att] <- n_attr
    
    # Generate permutation significance thresholds
    perm <- run_permutation_test(df_att, PRS_VARS, n_perms = N_PERMUTATIONS)
    thr <- quantile(perm, 0.05)
    perm_thresholds[att] <- thr
    
    # Run multivariate regressions
    results <- run_prs_models(df_att, att)
    results_list[[att]] <- results
  }
  
  results_all <- bind_rows(results_list)
  
  # Format plot categories
  attr_labels <- paste0(names(attr_counts), " (n=", attr_counts, ")")
  names(attr_labels) <- names(attr_counts)
  
  perm_df <- data.frame(
    Selection = names(perm_thresholds),
    perm_thr = as.numeric(perm_thresholds),
    stringsAsFactors = FALSE
  )
  
  results_all_stats <- results_all %>%
    left_join(perm_df, by = "Selection") %>%
    mutate(
      Selection_label = attr_labels[Selection],  
      perm_sig = p.value < perm_thr,
      disorder = dname,
      n_total = n_total_disorder
    )
  
  all_results[[dname]] <- results_all_stats
}


# ==============================================================================
# 5. OUTPUT: Save & Save Plots
# ==============================================================================

# Save continuous stats model dataset to disk
all_results_df <- dplyr::bind_rows(all_results)
readr::write_csv(all_results_df, "coombes_all_results_psych_classic.csv")

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

# Generate forest dot plot layouts
pos <- position_dodge(width = 0.9)

for (dname in unique(all_results_df$disorder)) {
  results_all_stats <- all_results_df %>% filter(disorder == dname)
  
  label_df <- results_all_stats %>% distinct(Selection, Selection_label)
  label_map <- setNames(label_df$Selection_label, label_df$Selection)
  
  results_all_stats <- results_all_stats %>%
    mutate(
      PRS = factor(PRS, levels = rev(PRS_VARS)),
      Selection = factor(Selection, levels = names(CUSTOM_COLORS))
    )
  
  plot_prs <- ggplot(
    results_all_stats,
    aes(x = estimate, y = PRS, color = Selection, group = Selection)
  ) +
    geom_errorbarh(
      aes(xmin = CI_lower, xmax = CI_upper),
      height = 0.35,
      position = pos
    ) +
    geom_point(
      aes(fill = perm_sig),
      shape = 21,
      size = 3.5,
      stroke = 1,
      position = pos
    ) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(
      yintercept = seq(1.5, length(PRS_VARS) - 0.5, by = 1),
      color = "grey85",
      linewidth = 1
    ) +
    scale_fill_manual(
      values = c("TRUE" = "black", "FALSE" = "white"),
      guide = "none"
    ) +
    scale_color_manual(
      values = CUSTOM_COLORS,
      labels = label_map,
      breaks = names(label_map)
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.1))) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0(dname, " (n = ", results_all_stats$n_total[1], " donors)"),
      x = "Change in log(odds) per SD increase in PRS",
      y = NULL,
      color = "Attribute"
    )
  
  # Save vector PDF graphics
  output_plot_path <- file.path(OUTPUT_DIR, paste0("PRS_plot_con_", dname, ".pdf"))
  ggsave(
    filename = output_plot_path,
    plot = plot_prs,
    width = 9,
    height = 0.9 * length(PRS_VARS) + 5,
    dpi = 300,
    device = cairo_pdf
  )
  
  print(plot_prs)
}