# ==============================================================================
# 1. PREDEFINED VARIABLES & CONFIGURATION
# ==============================================================================

# --- Cohort Selection Control Switch ---
# Options: "AD", "FTD", "PD"
COHORT_SELECTION     <- "AD" 

# --- Graphical Output Layout Constants ---
SIGNIFICANCE_LINE    <- 5e-8
THIN_THRESH_MANHATTAN <- 0.2
THIN_THRESH_QQ       <- 2.5
POINTS_PER_BIN_VAL   <- 30
LINE_WIDTH_VAL       <- 0.5

# --- Image Output Dimensions ---
DIM_WIDTH_PDF        <- 20
DIM_HEIGHT_PDF       <- 8
DIM_WIDTH_PNG        <- 1700
DIM_HEIGHT_PNG       <- 1600
RESOLUTION_DPI       <- 300

# --- File Directories & Input Data Path Map ---
FIGURES_OUT_DIR      <- "output"

DATA_PATH_MAP <- list(
  "AD"  = "/home/../NNDGWAS/ADCON2/REDUCED_GWAS_AD_CON_AGE_SEX_PC123.txt",
  "FTD" = "/home/../NNDGWAS/FTDCON2/REDUCED_GWAS_FTD_CON_AGE_SEX_PC123.txt",
  "PD"  = "/home/..t/NNDGWAS/PDCON/REDUCED_GWAS_PD_CON_AGE_SEX_PC123.txt"
)

TITLE_MAP <- list(
  "AD"  = "AD CON",
  "FTD" = "FTD CON",
  "PD"  = "PD CON"
)

# --- Package Dependency Layer ---
suppressPackageStartupMessages({
  library(GWASTools)
})


# ==============================================================================
# 2. REQUIRED DATA (INPUT)
# ==============================================================================
"""
REQUIRED INPUT DATA:
1. `gwas_raw`: Tabulated space/tab-delimited text summary statistics matrix.
2. The file structure must contain the following input columns:
   - 'X.CHROM' or 'CHR' : Chromosome identification mapping values.
   - 'P'                : Unadjusted asymptotic association p-value column.
"""

INPUT_FILE_PATH <- DATA_PATH_MAP[[COHORT_SELECTION]]
TARGET_LABEL    <- TITLE_MAP[[COHORT_SELECTION]]

if (!file.exists(INPUT_FILE_PATH)) {
  stop(paste("CRITICAL ERROR: Selected summary statistics file path missing:", INPUT_FILE_PATH))
}

gwas_raw <- read.table(INPUT_FILE_PATH, header = TRUE, comment.char = "")


# ==============================================================================
# 3. HELPER FUNCTIONS
# ==============================================================================
# (None required. Pipeline utilizes natively optimized GWASTools core bindings.)


# ==============================================================================
# 4. PROCESSES
# ==============================================================================

# --- Process 1: Standardizing Table Metadata Headers ---
if ("X.CHROM" %in% colnames(gwas_raw)) {
  colnames(gwas_raw)[colnames(gwas_raw) == "X.CHROM"] <- "CHR"
}

# --- Process 2: Cleaning Missing Row Vectors (QC) ---
gwas_data <- gwas_raw[complete.cases(gwas_raw[, c("P", "CHR")]), ]


# ==============================================================================
# 5. OUTPUT: Visualizations & Saving File
# ==============================================================================

if (!dir.exists(FIGURES_OUT_DIR)) {
  dir.create(FIGURES_OUT_DIR, recursive = TRUE)
}

# --- Set Up Alternating Chromosome Color Vectors ---
chr_colors <- rep(c("cadetblue2", "darkslategrey"), length.out = length(unique(gwas_data$CHR)))
col_vec    <- chr_colors[as.numeric(as.factor(gwas_data$CHR)) %% 2 + 1]

# ------------------------------------------------------------------------------
# Endpoint 1: Manhattan Plot Execution (PDF File Format)
# ------------------------------------------------------------------------------
manhattan_out_path <- file.path(FIGURES_OUT_DIR, paste0(tolower(COHORT_SELECTION), "gwas.pdf"))

pdf(manhattan_out_path, width = DIM_WIDTH_PDF, height = DIM_HEIGHT_PDF)

par(
  mar      = c(5, 5, 4, 2),
  cex.main = 8,
  cex.lab  = 6,
  cex.axis = 6
)
par(cex = 0.6, ps = 4)

manhattanPlot(
  p             = gwas_data$P,
  chromosome    = gwas_data$CHR,
  main          = paste("Association Results for", TARGET_LABEL),
  thinThreshold = THIN_THRESH_MANHATTAN,
  pointsPerBin  = POINTS_PER_BIN_VAL,
  signif        = SIGNIFICANCE_LINE,
  col           = col_vec
)

dev.off()


# ------------------------------------------------------------------------------
# Endpoint 2: Quantile-Quantile (QQ) Plot Execution (PNG File Format)
# ------------------------------------------------------------------------------
qq_out_path <- file.path(FIGURES_OUT_DIR, paste0("qqplotpc_", tolower(COHORT_SELECTION), ".png"))

png(qq_out_path, width = DIM_WIDTH_PNG, height = DIM_HEIGHT_PNG, res = RESOLUTION_DPI)

par(cex = 0.6)

qqPlot(
  pval          = gwas_data$P,
  thinThreshold = THIN_THRESH_QQ,
  main          = paste("qq plot Results for", TARGET_LABEL),
  lwd           = LINE_WIDTH_VAL
)

dev.off()


# --- Console Logging Confirmations ---
cat(paste0("✅ Pipeline operations completed for cohort: ", COHORT_SELECTION, "\n"))
cat(paste0(" -> Manhattan Plot: ", manhattan_out_path, "\n"))
cat(paste0(" -> QQ Plot:        ", qq_out_path, "\n"))