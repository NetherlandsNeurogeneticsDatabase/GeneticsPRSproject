#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# ==============================================================================
# 1. USER CONFIGURATION (INDICATIVE PATHS - REPLACE WITH YOUR ACTUAL PATHS)
# ==============================================================================

# Pipeline Root and Workspace Directories
PROJECT_ROOT_DIR="/path/to/your/project/prs_metabolite"
PREPROCESSED_DATA_DIR="${PROJECT_ROOT_DIR}/prep"
FINAL_OUTPUT_SCORE_DIR="${PROJECT_ROOT_DIR}/score"

# Input Summary Statistics (GWAS Data)
RAW_GWAS_SUMMARY_STATS="${PROJECT_ROOT_DIR}/input_gwas_summary_statistics.tsv"
SORTED_GWAS_SUMMARY_STATS="${PREPROCESSED_DATA_DIR}/gwas_summary_statistics_sorted.tsv"
STANDARDIZED_MA_FILE="${PREPROCESSED_DATA_DIR}/gwas_converted_format.ma"

# SBayesRC Reference Resources
SBAYESRC_SNP_WHITE_LIST="/path/to/reference/LDfolder/SNPsForSBayesRC.txt"
SBAYESRC_SNP_METADATA_INFO="/path/to/reference/LDfolder/snp.info"
SBAYESRC_LD_REFERENCE_DIR="/path/to/reference/LDfolder/"
SBAYESRC_FUNCTIONAL_ANNOTATION_FILE="/path/to/reference/annotfolder/annot_baseline2.2.txt"

# SBayesRC Analysis Parameters
GWAS_SAMPLE_SIZE="300139"
GWAS_TRAIT_TYPE="continuous" # Options: continuous, binary

# Target Cohort Genotype Data (PLINK2 Format)
TARGET_GENOTYPE_PGEN="/path/to/cohort/genotypes/imputed_dataset.pgen"
TARGET_GENOTYPE_FAM="/path/to/cohort/genotypes/imputed_dataset.fam"
TARGET_GENOTYPE_BIM="/path/to/cohort/genotypes/imputed_dataset.bim"

# Executable Tool Binaries
PLINK2_BINARY_PATH="/path/to/bioinformatics_tools/plink2/plink2"

# Output Prefix Targets
SBAYESRC_MODEL_OUTPUT_PREFIX="${PREPROCESSED_DATA_DIR}/sbayesrc_adjusted_effects"
FINAL_PRS_SCORE_PREFIX="${FINAL_OUTPUT_SCORE_DIR}/calculated_prs_profile"

# Ensure runtime directories exist
mkdir -p "$PREPROCESSED_DATA_DIR" "$FINAL_OUTPUT_SCORE_DIR"

# ==============================================================================
# 2. PIPELINE EXECUTION STEPS
# ==============================================================================

echo "=== Step 0: Sorting raw summary statistics by Genomic Coordinates ==="
{ 
  head -n 1 "$RAW_GWAS_SUMMARY_STATS"; 
  tail -n +2 "$RAW_GWAS_SUMMARY_STATS" | sort -k1,1n -k2,2n; 
} > "$SORTED_GWAS_SUMMARY_STATS"


echo "=== Step 1: Standardizing Summary Statistics format for SBayesRC ==="
time Rscript step1_prepareSummaryStatistics.R \
    "$SBAYESRC_SNP_WHITE_LIST" \
    "$SORTED_GWAS_SUMMARY_STATS" \
    "$STANDARDIZED_MA_FILE" \
    variant_id chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error p_value \
    "$GWAS_SAMPLE_SIZE" \
    "$GWAS_TRAIT_TYPE"


echo "=== Step 2: Harmonizing Variant Identifiers with Reference Panel ==="
time Rscript renameSNPIDs.R "$STANDARDIZED_MA_FILE" "$SBAYESRC_SNP_METADATA_INFO" 1 TRUE


echo "=== Step 3: Running SBayesRC - Quality Control / Tidy ==="
time Rscript -e "SBayesRC::tidy(mafile='$STANDARDIZED_MA_FILE', LDdir='$SBAYESRC_LD_REFERENCE_DIR', output='${PREPROCESSED_DATA_DIR}/dataset_tidy.ma', log2file=TRUE)"


echo "=== Step 4: Running SBayesRC - Missing Variant Imputation ==="
time Rscript -e "SBayesRC::impute(mafile='${PREPROCESSED_DATA_DIR}/dataset_tidy.ma', LDdir='$SBAYESRC_LD_REFERENCE_DIR', output='${PREPROCESSED_DATA_DIR}/dataset_imp.ma', log2file=TRUE)"


echo "=== Step 5: Running SBayesRC - Core Risk Effect Size Modeling ==="
time Rscript -e "SBayesRC::sbayesrc(mafile='${PREPROCESSED_DATA_DIR}/dataset_imp.ma', LDdir='$SBAYESRC_LD_REFERENCE_DIR', outPrefix='$SBAYESRC_MODEL_OUTPUT_PREFIX', annot='$SBAYESRC_FUNCTIONAL_ANNOTATION_FILE', log2file=TRUE)"


echo "=== Step 6: Computing Individual Profile Polygenic Risk Scores (PLINK2) ==="
time "$PLINK2_BINARY_PATH" \
    --pgen "$TARGET_GENOTYPE_PGEN" \
    --fam "$TARGET_GENOTYPE_FAM" \
    --bim "$TARGET_GENOTYPE_BIM" \
    --score "${SBAYESRC_MODEL_OUTPUT_PREFIX}.txt" 1 2 3 header center \
    --out "$FINAL_PRS_SCORE_PREFIX"

echo "=== Pipeline Workflow Finished Successfully ==="