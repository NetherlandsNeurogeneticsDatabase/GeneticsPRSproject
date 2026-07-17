# ==============================================================================
# 1. PREDEFINED VARIABLES & CONFIGURATION
# ==============================================================================

# --- Analysis & Filtering Parameters ---
P_VALUE_THRESHOLD  <- 1e-10
PP4_AD_THRESH      <- 0.5
PP4_FTD_THRESH     <- 0.4
PP4_PD_THRESH      <- 0.5

# --- Visual & Layout Settings ---
PLOT_WIDTH_PDF     <- 7
PLOT_HEIGHT_PDF    <- 20
PLOT_WIDTH_PNG     <- 900
PLOT_HEIGHT_PNG    <- 2000
PLOT_RES_PNG       <- 150

CELL_WIDTH         <- 20
CELL_HEIGHT        <- 15
FONT_SIZE_ROW      <- 12
FONT_SIZE_COL      <- 14
FONT_SIZE_GLOBAL   <- 10

# --- Labels & Color Schemes ---
PLOT_MAIN_TITLE    <- "Coloc PP4 Schwartzentruber (AD) vs. NND AD, NND FTD, NND PD"
COLUMN_LABELS      <- c("PP4 AD NND", "PP4 FTD NND", "PP4 PD NND")
COLOR_PALETTE_REV  <- "RdYlBu"
COLOR_STEPS        <- 100

# --- File Paths ---
INPUT_EXCEL_PATH   <- "/home/jupyter-n.mekkes@gmail.com-f6d87/thesis_ch4_genetics/loci_of_interest/EMS118040-supplement-Supplementary_Tables_1_14.xlsx"
OUTPUT_DIR         <- "/home/jupyter-n.mekkes@gmail.com-f6d87/thesis_ch4_genetics/loci_of_interest/coloc"
OUTPUT_FILE_PREFIX <- "heatmap_output_AD"

PDF_OUTPUT_PATH    <- file.path(OUTPUT_DIR, paste0(OUTPUT_FILE_PREFIX, ".pdf"))
PNG_OUTPUT_PATH    <- file.path(OUTPUT_DIR, paste0(OUTPUT_FILE_PREFIX, ".png"))


# ==============================================================================
# 2. REQUIRED DATA (INPUT)
# ==============================================================================
"""
REQUIRED INPUT DATA:
1. `ad_loci`: Master summary matrix DataFrame imported from the supplementary spreadsheet.
2. The sheet structure must contain the following input columns:
   - 'Chr'                  : Chromosome location descriptor.
   - 'Lead p'               : Association p-value metric.
   - 'Lead SNP pos'         : Genome basepair coordinate index.
   - 'PP4_ad_public_ad_nnd' : Posterior probability value for AD coloc.
   - 'PP4_ad_public_ftd_nnd': Posterior probability value for FTD coloc.
   - 'PP4_ad_public_pd_nnd' : Posterior probability value for PD coloc.
"""

# --- Package Dependency Layer ---
suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
})

# --- Verify and Read Input Matrix ---
if (!file.exists(INPUT_EXCEL_PATH)) {
  stop(paste("CRITICAL ERROR: Input dataset missing at coordinate path:", INPUT_EXCEL_PATH))
}

ad_loci <- read_excel(INPUT_EXCEL_PATH, sheet = '1-AD loci')


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================

save_pheatmap_pdf <- function(x, filename, width = 7, height = 20) {
  # Directs the pheatmap grid output to an editable vector PDF canvas
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width = 900, height = 2000, res = 150) {
  # Directs the pheatmap grid output to a pixel-scaled high-DPI PNG format canvas
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


# ==============================================================================
# 4. PROCESSES
# ==============================================================================

# --- Process 1: Statistical Threshold Subsetting ---
ad_loci_fitered <- ad_loci %>%
  filter(`Lead p` <= P_VALUE_THRESHOLD)

ad_loci_heatmap <- ad_loci_fitered %>%
  filter(
    (PP4_ad_public_ad_nnd > PP4_AD_THRESH) |
    (PP4_ad_public_ftd_nnd > PP4_FTD_THRESH) |
    (PP4_ad_public_pd_nnd > PP4_PD_THRESH)
  )

# Normalize column strings to match downstream vector transformations
if ("Lead SNP pos" %in% names(ad_loci_heatmap)) {
  names(ad_loci_heatmap)[names(ad_loci_heatmap) == "Lead SNP pos"] <- "lead_snp_pos"
}


# --- Process 2: Framing Heatmap Input Vectors ---
heatmap_df <- data.frame(
  ad_loci_heatmap[, c(
    "lead_snp_pos",
    "Chr",
    "PP4_ad_public_ad_nnd",
    "PP4_ad_public_ftd_nnd",
    "PP4_ad_public_pd_nnd"
  )],
  stringsAsFactors = FALSE
)

# Convert row addresses to target [Chr_Position] character matrices
rownames(heatmap_df) <- paste0(
  as.integer(heatmap_df$Chr), "_",
  as.integer(heatmap_df$lead_snp_pos)
)

# Sort genome rows chronologically across chromosome boundaries
heatmap_df <- heatmap_df[order(heatmap_df$Chr, heatmap_df$lead_snp_pos), ]

# Drop auxiliary sorting variables to isolate matrix profiles exclusively
heatmap_df <- heatmap_df[, !(names(heatmap_df) %in% c("Chr", "lead_snp_pos"))]


# ==============================================================================
# 5. OUTPUT: Formatting, Rendering & File Export
# ==============================================================================

# Ensure target folder exists on disk
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# Configure color palette scales
my_colors <- colorRampPalette(rev(brewer.pal(7, COLOR_PALETTE_REV)))(COLOR_STEPS)
my_breaks <- seq(0, 1, length.out = COLOR_STEPS + 1)

# Render core pheatmap structure
p <- pheatmap(
  heatmap_df,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth    = CELL_WIDTH,
  cellheight   = CELL_HEIGHT,
  main         = PLOT_MAIN_TITLE,
  color        = my_colors,
  breaks       = my_breaks,    
  fontsize_row = FONT_SIZE_ROW,
  fontsize_col = FONT_SIZE_COL,
  fontsize     = FONT_SIZE_GLOBAL,
  legend       = TRUE,
  labels_col   = COLUMN_LABELS,
  angle_col    = 90
)

# Save structured canvas allocations to files
save_pheatmap_pdf(p, PDF_OUTPUT_PATH, width = PLOT_WIDTH_PDF, height = PLOT_HEIGHT_PDF)
save_pheatmap_png(p, PNG_OUTPUT_PATH, width = PLOT_WIDTH_PNG, height = PLOT_HEIGHT_PNG, res = PLOT_RES_PNG)

cat("Process completed successfully.\n")