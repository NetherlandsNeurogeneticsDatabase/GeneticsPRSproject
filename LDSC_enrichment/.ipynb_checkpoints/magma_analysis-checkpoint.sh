#!/usr/bin/env bash
# ==============================================================================
# 1. PREDEFINED VARIABLES & CONFIGURATION
# ==============================================================================

# --- Executable Binaries ---
MAGMA_BIN="./magma_v1.10_mac/magma"

# --- Genome Reference Data Configuration ---
GENE_LOC_REF="data/NCBI37.3/NCBI37.3.gene.loc"
REFERENCE_BFILE="data/g1000_eur/g1000_eur"

# --- Input Summary Statistics & Metadata ---
GWAS_SNP_LOC="data/sumstats/GCST90012877.snploc"
GWAS_PVAL_DATA="data/sumstats/GCST90012877.pval"
GWAS_SAMPLE_SIZE=472868

# --- Gene Set / Pathway Mapping Data ---
RAW_GMT_PATH="data/gene_annot/c5.all.v2025.1.Hs.entrez.gmt"
CLEAN_GMT_PATH="data/gene_annot/5.all.v2025.1.Hs.entrez.clean.gmt"
FINAL_SET_ANNOT_PATH="data/gene_annot/5.all.v2025.1.Hs.entrez.clean.tab.gmt"

# --- Output Pipeline Directories ---
RESULTS_DIR="results"
ANNOT_OUT_PREFIX="${RESULTS_DIR}/AD_zwarschentruber_annotation"
GENE_OUT_PREFIX="${RESULTS_DIR}/my_gwas_gene_analysis"
PATHWAY_OUT_PREFIX="${RESULTS_DIR}/AD_zwarschentruber_gene_analysis"

# ==============================================================================
# 2. REQUIRED DATA (INPUT)
# ==============================================================================
# Ensure the following prerequisite input structures are prepared:
#
# 1. GWAS SNP Location File (GWAS_SNP_LOC): Tab-separated, no-header format containing [rs_id, chr, pos].
# 2. GWAS P-Value Summary File (GWAS_PVAL_DATA): Structured document containing [rs_id, P-value].
# 3. Download link for the GSEA Molecular Signatures Gene-Set Collection:
#    https://www.gsea-msigdb.org/gsea/msigdb/download_file.jsp?filePath=/msigdb/release/2025.1.Hs/c5.all.v2025.1.Hs.entrez.gmt
# 4. Download link for MAGMA reference panels (REFERENCE_BFILE): https://cncr.nl/research/magma/

# ==============================================================================
# 3. PROCESSES
# ==============================================================================

# Ensure directory pathways exist
mkdir -p "${RESULTS_DIR}"
mkdir -p "$(dirname "${CLEAN_GMT_PATH}")"

# ------------------------------------------------------------------------------
# Process 1: Annotation Mapping
# ------------------------------------------------------------------------------
# Reads raw variant locations from the snp-loc file, associates them against 
# official Entrez coordinate models, and produces structured variant annotations.
echo "Running Process 1: Variant-to-Gene Mapping Annotation..."
"${MAGMA_BIN}" \
  --annotate \
  --snp-loc "${GWAS_SNP_LOC}" \
  --gene-loc "${GENE_LOC_REF}" \
  --out "${ANNOT_OUT_PREFIX}"


# ------------------------------------------------------------------------------
# Process 2: Gene-Level Analysis
# ------------------------------------------------------------------------------
# Aggregates raw SNP-level association profiles into explicit gene-level units 
# via standard mean-based projection models utilizing an external reference panel.
echo "Running Process 2: Genic Statistical Aggregation..."
"${MAGMA_BIN}" \
  --bfile "${REFERENCE_BFILE}" \
  --pval "${GWAS_PVAL_DATA}" N="${GWAS_SAMPLE_SIZE}" use=SNP,P \
  --gene-annot "${ANNOT_OUT_PREFIX}.genes.annot" \
  --gene-model snp-wise=mean \
  --out "${GENE_OUT_PREFIX}"


# ------------------------------------------------------------------------------
# Process 3: Functional Pathway Cleansing & Transformation
# ------------------------------------------------------------------------------
# Strips hyper-linked indexing and standardizes Molecular Signatures Database (MSigDB) 
# parsing matrices to strictly output tab-separated pathway definitions.
echo "Running Process 3: Cleansing Pathway Definitions Database..."

awk -F'\t' '{
    # Extract structural pathway name, then iterate sequentially across functional genes
    printf "%s", $1
    for (i=3; i<=NF; i++) {
        printf "\t%s", $i
    }
    printf "\n"
}' "${RAW_GMT_PATH}" > "${CLEAN_GMT_PATH}"

# Force space-delimited text models to strict tabbed format tables
awk '{$1=$1; print}' OFS='\t' "${CLEAN_GMT_PATH}" > "${FINAL_SET_ANNOT_PATH}"


# ------------------------------------------------------------------------------
# Process 4: Gene-Set Pathway Analysis
# ------------------------------------------------------------------------------
# Performs generalized linear regressions to evaluate whether specific 
# biological pathways carry enriched risk signals compared to background genomes.
echo "Running Process 4: Assessing Functional Over-Representation Analytics..."
"${MAGMA_BIN}" \
  --gene-results "${GENE_OUT_PREFIX}.genes.raw" \
  --set-annot "${FINAL_SET_ANNOT_PATH}" \
  --out "${PATHWAY_OUT_PREFIX}"

# ==============================================================================
# 4. OUTPUT
# ==============================================================================
# The following final summary data elements are written to disk:
#
# - ${PATHWAY_OUT_PREFIX}.gsa.out   : Analytical pathway / gene-set statistical p-values.
# - ${PATHWAY_OUT_PREFIX}.gsa.raw   : Accompanying pathway covariate tracking residuals matrix.
# - ${GENE_OUT_PREFIX}.genes.out    : Core gene-centric target testing values.

echo "✅ Pipeline successfully executed. Results safely routed to '${RESULTS_DIR}/' directory."