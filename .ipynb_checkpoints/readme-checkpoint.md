# Neurogenomics Analysis Pipeline

Analysis pipeline accompanying the neurogenomics study integrating genotype, neuropathology, and clinical data from the Netherlands Neurogenomics Database (NND).

---

# 1. Data Integration & LLM Pipeline

## Neuropathology Feature Extraction
Converts unstructured neuropathology text reports into structured quantitative phenotypes using Large Language Models (LLMs).

* **Input:** Raw clinical/neuropathological text reports from the NND.
* **Function/Tool:** LLM pipeline script (`llm_pipeline`).
* **Output:** Structured pathology feature tables (e.g., Braak stage, CERAD score, regional pathology burden).

---

# 2. Genomic Quality Control & Stratification

## SNP Genotype Quality Control
Preprocesses and filters raw genomic data to remove systemic artifacts and low-quality samples/variants.

* **Input:** Raw genotype dataset (VCF/PLINK format).
* **Function/Tool:** `PLINK2`, `zcall`, `shape-it`, and `minimac4`.
* **Key Filters Applied:** * Excess heterozygosity removal (Hardy-Weinberg equilibrium threshold $1 \times 10^{-5}$)
    * Low call rate removal ($<97.5\%$ for variants and samples via `zcall`)
    * Sex mismatch removal (`shape-it` & `minimac4`)
    * Cryptic relatedness filtering (removal of closely related donors)
    * Pre-association filters (Imputation score $>0.8$, MAF $>0.01$, HWE)
* **Output:** Cleaned and imputed genotype dataset.

## PCA & Population Stratification
Controls for ancestry-driven population structure across regional cohorts.

* **Input:** Cleaned genotype dataset.
* **Function/Tool:** `PLINK2` (MAF $\ge 0.01$, LD pruning $r^2 = 0.2$, PCA calculation), `Kruskal-Wallis` test, and `geopandas` visualization.
* **Process:** Computes principal components (PCs) to strip out the first four ancestry components. Runs Kruskal-Wallis tests to detect unique PCs matching approximate geographical addresses across Dutch provinces, correcting via False Discovery Rate (FDR $<0.05$).
* **Output:** Sig-PC clusters and geographical distribution heatmaps mapped by normalized donor counts per 100k.

---

# 3. Public GWAS Summary Statistics Refinement

## Preprocessing & Standardization
Standardizes heterogeneous external public GWAS files into uniform layouts before genetic risk estimation.

* **Input:** Public GWAS summary statistics across selected use cases (FTLD-MND, Neuropsychiatric symptoms, Neuropathologically defined diagnosis).
* **Function/Tool:** Custom processing script.
* **Process:** Standardizes data schemas to GRCh37 coordinates, sorts by chromosome/position, and infers missing parameters (e.g., sample size per SNP derived from total sample size; Effect Allele Frequency (EAF) mapped via SbayesRC LD references; Standard Error (SE) estimated using the natural log of odds ratios).
* **Output:** Harmonized, sorted GWAS summary statistic tables containing standard genomic headers.

---

# 4. Genetic Correlation & Enrichment

## Matrix-Based Co-regulation Analysis
Quantifies global genetic overlaps and regional functional enrichment across traits.

* **Input:** Cleaned NND genotype data and preprocessed public GWAS summary statistics.
* **Function/Tool:** `LDSC` (Linkage Disequilibrium Score Regression) and `MAGMA`.
* **Process:** 1. Runs global pairwise genetic correlation ($r_g$) via LDSC.
    2. Runs gene/pathway-level enrichment mapping across specific Gene Ontology blocks (GOBP/GOCC/GOMF) via MAGMA.
    3. Transforms long-form statistical records into structured matrices.
* **Output:** Genetic correlation tables, MAGMA enrichment files, and matrix-based $-\log_{10}(\text{P-value})$ heatmaps.

---

# 5. Polygenic Risk Scores (PRS) Execution

## PRS Calculation & Normalization
Generates personalized genomic liability metrics across disease profiles.

* **Input:** Cleaned genotype data and standardized public GWAS summary statistics.
* **Function/Tool:** `SbayesRC` (SNP risk effect estimation) and `PLINK2` (profile score multiplication).
* **Process:** Computes raw individual risk burdens, standardizes dosage scores (centered to zero), and applies Z-score transformations.
* **APOE Splitting Feature:** When isolated profiling is required, the APOE gene boundary (Chr 19: 44,409,039-) is trimmed away to record independent background PRSs, while tracking the standalone $e3/e4$ risk haplotype dose metrics independently.
* **Output:** Normalized individual PRS vector profiles.

---

# 6. Downstream Main Analysis & Visualizations

## PRS Association Engine
Executes statistical tests mapping genetic boundaries to endophenotypes, symptoms, and clinical labels.

* **Input:** Individual PRS files, structured LLM-derived pathology data, and clinical diagnostic tables.
* **Function/Tool:** Multi-test regression models, Permutation-based resampling mechanism (10k shuffles for empirical p-values), and Benjamini-Hochberg FDR correction.
* **Output:** Publication-ready visual assets saved into targeted outputs:
    * **Figure: Dot Plots** / **Figure: Scatter Plots** (Case vs. control mean PRS differences; PRS vs Neuropathology load maps)
    * **Figure: Heatmap - Feature Association** (Clustered disease/shade intensity configurations)

> **Note:** All script paths and file directories within these workflows are indicative. Please update configuration paths to match your unique deployment footprint.

---

# Software

## Genetics
* PLINK 2
* bcftools
* IMPUTE5
* VCFtools
* LDSC
* MAGMA

## Python
* pandas
* NumPy
* SciPy
* scikit-learn
* statsmodels
* geopandas
* matplotlib / seaborn

---

# Citation

```bibtex
@article{mekkes2026dynamic,
  title={Dynamic interplay of polygenic risk across brain disorders, neuropathological endophenotypes, and neuropsychiatric symptoms},
  author={Mekkes, Nienke Jacobine and Kumar, Shivam and Hoekstra, Eric and Marmolejo Garza, Alejandro and Dagkesamanskaia, Ekaterina and Wever, Dennis and Groot, Minke and Kreft, Karim L and Rajicic, Ana and Seelaar, Harro and others},
  journal={medRxiv},
  pages={2026--06},
  year={2026},
  publisher={Cold Spring Harbor Laboratory Press}
}