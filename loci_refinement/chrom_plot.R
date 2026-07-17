# ==============================================================================
# 1. PREDEFINED VARIABLES & CONFIGURATION
# ==============================================================================

# --- Selected Locus File Index ---
# Explicitly isolates target CSV index 34 to generate the desired layout
CHOSEN_LOCUS_INDEX   <- 34
ID_COLUMN_KEY        <- "ID_ad_public"

# --- Plot Aesthetics & Scaling Parameters ---
GENOME_BUILD         <- 37
SIGNIFICANCE_THRESH  <- 0.01
POINT_SIZE_VAL       <- 2
TRACKS_GENEBAR       <- 2
SCALE_ASSOC_PLOT     <- 20
HEIGHT_PLOT_PARAM    <- 44
X_AXIS_BREAKS        <- 6

# --- Plot Dimension Configurations ---
GG_SAVE_WIDTH        <- 8
GG_SAVE_HEIGHT       <- 10

# --- High-Throughput Stack Trait Metadata ---
STACKED_TRAITS       <- c("Public AD GWAS", "NND AD GWAS", "NND FTD GWAS", "NND PD GWAS")
HIGHLIGHT_COLOURS    <- c("blue", "yellow", "blue")

# --- File Paths ---
FOLDER_PATH          <- "/home/jupyter-n.mekkes@gmail.com-f6d87/thesis_ch4_genetics/loci_of_interest/relaxed_ad_loci2"
TARGET_OUTPUT_PLOT   <- "figures/chrom19zoomplot.pdf"

# --- Package Dependency Layer ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(geni.plots)
  library(LDlinkR)
  library(stringr)
  library(ggplot2)
})


# ==============================================================================
# 2. REQUIRED DATA (INPUT)
# ==============================================================================
"""
REQUIRED INPUT DATA:
1. `csv_files`: List of localized variant mapping collections matching '\\AD_schwartzentruber.csv$'.
2. The data structure within each file must feature the following columns:
   - 'CHR'              : Chromosomal mapping index.
   - 'POS'              : Base pair genomic coordinate location.
   - 'ID_ad_public'     : Lead public reference identifier mapping name.
   - 'P_ad_public'      : Association p-value derived from public reference cohorts.
   - 'P_AD'             : Cleaned dataset association p-value matching AD tracking.
   - 'P_FTD'            : Cleaned dataset association p-value matching FTD tracking.
   - 'P_PD'             : Cleaned dataset association p-value matching PD tracking.
   - 'chrom_coordinate' : Strict location string used by the downstream LDlinkR API context.
"""

# Fetch locus tracking records from directory path
if (!dir.exists(FOLDER_PATH)) {
  stop(paste("CRITICAL ERROR: Configuration pathway directory missing at location:", FOLDER_PATH))
}

csv_files <- list.files(path = FOLDER_PATH, pattern = "\\AD_schwartzentruber.csv$", full.names = TRUE)
target_csv_file <- csv_files[CHOSEN_LOCUS_INDEX]


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

process_loci <- function(file_path, id_column = "ID_ad_public") {
  # Flags, standardizes, and populates null variants with uniform tracking IDs
  loci_of_interest <- read.csv(file_path, header = TRUE, sep = ",")
  na_indices <- which(is.na(loci_of_interest[[id_column]]) | loci_of_interest[[id_column]] == "")
  new_names <- paste0("nameless", seq_along(na_indices))
  loci_of_interest[[id_column]][na_indices] <- new_names
  return(loci_of_interest)
}

get_ld_matrix_safe <- function(chrom_coordinates, max_retries = 5) {
  # Performs structural API matrix extraction queries with a failure protection layer
  for (attempt in 1:max_retries) {
    tryCatch({
      result <- LDmatrix(
        snps = chrom_coordinates,
        pop = c("CEU", "IBS", "FIN", "GBR", "TSI"),
        r2d = "r2",
        token = "23edcb3e3ad7",
        file = FALSE,
        genome_build = "grch37",
        api_root = "https://ldlink.nih.gov/LDlinkRest"
      )
      return(result)
    }, error = function(e) {
      if (attempt < max_retries) {
        Sys.sleep(5)
      } else {
        stop("API connection request failure limit reached after multiple attempts.")
      }
    })
  }
}

process_ld_matrix <- function(loci_of_interest, loci_interest_ldmatrix, id_column = "ID_ad_public") {
  # Synchronizes dimension alignments between the summary framework and LD matrix outputs
  loci_interest_ldmatrix_allrows <- loci_of_interest %>%
    select(all_of(id_column)) %>%
    left_join(loci_interest_ldmatrix, by = setNames("RS_number", id_column))
  
  markers <- loci_interest_ldmatrix_allrows[[id_column]]
  row.names(loci_interest_ldmatrix_allrows) <- loci_interest_ldmatrix_allrows[[1]]
  loci_interest_ldmatrix_allrows <- loci_interest_ldmatrix_allrows[-1]
  
  missing_columns <- setdiff(markers, names(loci_interest_ldmatrix_allrows))
  for (col in missing_columns) {
    loci_interest_ldmatrix_allrows[[col]] <- NA
  }
  
  loci_interest_ldmatrix_allrows_allcols <- loci_interest_ldmatrix_allrows[, markers]
  loci_interest_ldmatrix_allrows_allcols <- as.matrix(loci_interest_ldmatrix_allrows_allcols)
  
  sqrt_matrix <- matrix(sqrt(loci_interest_ldmatrix_allrows_allcols),
                        nrow = nrow(loci_interest_ldmatrix_allrows_allcols))
  rownames(sqrt_matrix) <- rownames(loci_interest_ldmatrix_allrows_allcols)
  colnames(sqrt_matrix) <- colnames(loci_interest_ldmatrix_allrows_allcols)
  
  return(sqrt_matrix)
}


# ==============================================================================
# 4. PROCESSES
# ==============================================================================

# --- Process 1: Framing File Metrics & Genome Coordinates ---
filename <- basename(target_csv_file)
parts <- str_split(filename, "_")[[1]]

chromosome_of_interest <- str_remove(parts[1], "chr")
lead_snp      <- parts[4]
pos_interest  <- as.numeric(parts[7])
start_snp     <- pos_interest - 25000
stop_snp      <- pos_interest + 25000

# Establish local naming metrics for asset tracking
assoc_filename <- sub("\\.csv$", "_assoc.csv", target_csv_file)
corr_filename  <- sub("\\.csv$", "_corr.csv", target_csv_file)

# --- Process 2: Coordinate Matrix Verification & LD Matrix Synchronization ---
if (!file.exists(assoc_filename) || !file.exists(corr_filename)) {
  loci_of_interest <- process_loci(target_csv_file, id_column = ID_COLUMN_KEY)
  
  # Trim tracking records strictly inside boundary constraints
  loci_of_interest <- loci_of_interest[
    loci_of_interest$CHR == chromosome_of_interest &
    loci_of_interest$POS > start_snp &
    loci_of_interest$POS < stop_snp, 
  ]
  
  assoc <- loci_of_interest
  corr  <- get_ld_matrix_safe(assoc$chrom_coordinate)
  corr  <- process_ld_matrix(assoc, corr, id_column = ID_COLUMN_KEY)
  
  write.csv(assoc, assoc_filename, row.names = FALSE)
  write.csv(corr, corr_filename, row.names = TRUE)
} else {
  assoc <- read.csv(assoc_filename, header = TRUE, sep = ",")
  corr  <- read.csv(corr_filename, header = TRUE, row.names = 1, sep = ",")
}

# --- Process 3: Standardizing Variable Labels for Stacked Regional Plots ---
merged_df_filtered <- assoc %>%
  rename(
    chr      = CHR,
    marker   = ID_ad_public,
    pvalue_1 = P_ad_public,
    pvalue_2 = P_AD,
    pvalue_3 = P_FTD,  
    pvalue_4 = P_PD,
    pos      = POS
  )

# Process adjustment statistics metrics
merged_df_filtered$qvalue    <- p.adjust(merged_df_filtered$pvalue_3, method = "BH")
merged_df_filtered$qvalue_AD <- p.adjust(merged_df_filtered$pvalue_2, method = "BH")
merged_df_filtered$qvalue_PD <- p.adjust(merged_df_filtered$pvalue_4, method = "BH")


# ==============================================================================
# 5. OUTPUT: Visualizations & Saving File
# ==============================================================================

# Ensure layout directories exist on disk before writing assets
if (!dir.exists(dirname(TARGET_OUTPUT_PLOT))) {
  dir.create(dirname(TARGET_OUTPUT_PLOT), recursive = TRUE)
}

# Identify regional min-peak indicators across cohorts
list_of_markers <- c(
  merged_df_filtered$marker[which.min(merged_df_filtered$pvalue_1)],
  merged_df_filtered$marker[which.min(merged_df_filtered$pvalue_2)],
  merged_df_filtered$marker[which.min(merged_df_filtered$pvalue_3)],
  merged_df_filtered$marker[which.min(merged_df_filtered$pvalue_4)]
)

top_marker_to_use <- if (lead_snp %in% merged_df_filtered$marker) lead_snp else NULL

# Construct multi-trait high-resolution canvas map
p <- fig_region_stack(
  data               = merged_df_filtered,
  traits             = STACKED_TRAITS,
  corr               = corr,
  build              = GENOME_BUILD,
  top_marker         = top_marker_to_use, 
  highlights         = list_of_markers,
  thresh             = SIGNIFICANCE_THRESH,
  point_size         = POINT_SIZE_VAL,  
  highlights_colours = HIGHLIGHT_COLOURS,  
  title_center       = TRUE,
  genebar_ntracks    = TRACKS_GENEBAR,
  assoc_plot_size    = SCALE_ASSOC_PLOT,
  plot_height        = HEIGHT_PLOT_PARAM
)

# Re-scale linear nucleotide ranges to mega-basepair units
p <- p + scale_x_continuous(
  labels   = function(x) sprintf("%.2f", x / 1e6),
  n.breaks = X_AXIS_BREAKS
)

# Export layout safely to the designated disk directory path
ggsave(TARGET_OUTPUT_PLOT, plot = p, width = GG_SAVE_WIDTH, height = GG_SAVE_HEIGHT)
print(p)