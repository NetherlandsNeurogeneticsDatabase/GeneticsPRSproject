# Neurogenomics Analysis Pipeline

Analysis pipeline accompanying the neurogenomics study integrating genotype, neuropathology, and clinical data from the Netherlands Neurogenomics Database (NND).

---

## Directory Structure

```text
├── GWAS_QC/
│   └── QC_plot.R
├── LDSC_enrichment/
│   ├── LDSC_correlation.ipynb
│   ├── MAGMA_visualization.ipynb
│   └── magma_analysis.sh
├── Neuropath_analysis/
│   ├── Generate_plot.ipynb
│   ├── Perform_analysis.ipynb
│   ├── input/
│   │   ├── config.json
│   │   └── prompt_abc.json
│   ├── llm_extract.py
│   ├── output/
│   └── utils.py
├── PCA/
│   ├── PCA_cluster_on_map.ipynb
│   └── output/
├── PRS_calculation/
│   ├── Data_prep.ipynb
│   ├── PRS_analysis_workflow.ipynb
│   ├── Symptoms_plot.R
│   ├── output/
│   ├── prs_calculation.sh
│   └── utils.py
├── loci_refinement/
│   ├── chrom_plot.R
│   └── heatmap_risk_loci.R
└── readme.md
```

---

## 1. Neuropathology Feature Extraction (`Neuropath_analysis/`)

Converts unstructured neuropathology text reports into structured quantitative phenotypes using Large Language Models (LLMs).

* **Scripts:** 
  * `llm_extract.py` — Core feature extraction engine.
  * `Perform_analysis.ipynb` — Downstream statistical analysis.
  * `Generate_plot.ipynb` — Visualizations for neuropathological features.
* **Input Data:** 
  * Raw clinical/neuropathological text reports from the NND.
  * Configuration files (`input/config.json`, `input/prompt_abc.json`).
  * *Protocol:* Files are formatted as a matrix containing an `id` column, an unstructured `sentences` column, and an empty feature name column.
* **Output Data:** 
  * Structured pathology feature tables exported into `output/` with columns populated based on explicit text criteria.
* **Associated Paper Figure:**
  * **Figure 5:** Neuropathology group analysis and structural phenotype mapping (contained within this separate modular folder workflow).

---

## 2. Genomic Quality Control & Stratification (`GWAS_QC/` & `PCA/`)

Preprocesses and filters genomic datasets to eliminate artifacts, while evaluating geographic-driven population stratifications.

### SNP Genotype Quality Control (`GWAS_QC/`)
* **Script:** `QC_plot.R`
* **Input Data:** Raw genotype datasets (VCF/PLINK format).
* **Key Filters Applied:** Excess heterozygosity removal (HWE threshold $1 \times 10^{-5}$), low call rates ($<97.5\%$ via `zcall`), sex mismatches (`shape-it` & `minimac4`), and cryptic relatedness filtering.
* **Output Data:** Cleaned, filtered, and imputed genomic matrices.
* **Associated Paper Figure:**
  * **Figure 1b:** Data distribution with a violin plot tracking donor age profiles.

### PCA & Population Stratification (`PCA/`)
* **Script:** `PCA_cluster_on_map.ipynb`
* **Input Data:** Cleaned genotype dataset from the QC workflow.
* **Output Data:** Ancestry PC components and distribution records saved inside `output/`.
* **Associated Paper Figure:**
  * **Figure 1c:** PCA cluster projection mapped onto geographic spatial plots of the Netherlands.

---

## 3. Public GWAS Summary Statistics Refinement (`loci_refinement/`)

Standardizes heterogeneous external public GWAS findings before down-stream hazard modeling and fine-mapping risk loci.

* **Scripts:** `chrom_plot.R`, `heatmap_risk_loci.R`
* **Input Data:** External public summary statistics across specific cohorts (FTLD-MND, neuropsychiatric profiles).
* **Output Data:** Unified summary statistics layouts containing standard genomic metrics.
* **Associated Paper Figure:**
  * **Figure 4:** Loci refinement diagrams, risk heatmaps, and target region definitions.

---

## 4. Genetic Correlation & Enrichment (`LDSC_enrichment/`)

Tracks global genetic overlaps and regional functional enrichment metrics across target endophenotypes.

* **Scripts:** `magma_analysis.sh`, `LDSC_correlation.ipynb`, `MAGMA_visualization.ipynb`
* **Input Data:** Cleaned NND genotype data and preprocessed public GWAS summaries.
* **Output Data:** Clustered statistical grids and matrix logs.
* **Associated Paper Figures:**
  * **Figure 2e & 2f:** Linkage Disequilibrium Score Regression (LDSC) correlation profiles and trait heatmap matrices.

---

## 5. Polygenic Risk Score Pipeline (`PRS_calculation/`)

Generates individualized genetic liability metrics across disease profiles to map clinical symptoms and mutations.

* **Scripts:** 
  * `prs_calculation.sh` & `Data_prep.ipynb` — Scoring setup and SbayesRC array profile calculations.
  * `PRS_analysis_workflow.ipynb` — Regression mapping, case-control distributions, and group analysis.
  * `Symptoms_plot.R` — Target visualization script isolating neuropsychiatric indicators.
* **Input Data:** Cleaned internal genotypes, diagnostic registries, and external summary statistics weights.
* **Output Data:** Normalized individual PRS vector tracks and regression diagnostics inside `output/`.
* **Associated Paper Figures:**
  * **Figure 2a, b, c:** Three distinct violin plots assessing specific polygenic risk score distributions.
  * **Figure 2d:** Dot plot showing associations between specific polygenic risk scores and diagnostic states.
  * **Figure 3a:** Dot plot mapping risk trends across Frontotemporal Dementia (FTLD) subgroups.
  * **Figure 3b, c, d, e:** Violin plots outlining polygenic scores stratified by explicit genetic mutation status (e.g., *C9orf72* status).
  * **Figure 6:** Symptoms analysis and dot plots connecting polygenic risk metrics with neuropsychiatric manifestations.

---

## Prerequisites & Installation

### Core Software Environment
Ensure your local compute node has access to both **Python (3.8+)** and **R (4.0.0+)** environments.

#### Required Command Line Tools:
* `PLINK 2`
* `bcftools`
* `IMPUTE5`
* `VCFtools`
* `LDSC`
* `MAGMA`

#### Python Library setup

```bash
pip install pandas numpy scipy scikit-learn statsmodels geopandas matplotlib seaborn


#### R Library Setup:
```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Biological data scaling packages
devtools::install_github("privefl/bigsnpr")
```

---

## Citation

```bibtex
@article{mekkes2026dynamic,
  title={Dynamic interplay of polygenic risk across brain disorders, neuropathological endophenotypes, and neuropsychiatric symptoms},
  author={Mekkes, Nienke Jacobine and Kumar, Shivam and Hoekstra, Eric and Marmolejo Garza, Alejandro and Dagkesamanskaia, Ekaterina and Wever, Dennis and Groot, Minke and Kreft, Karim L and Rajicic, Ana and Seelaar, Harro and others},
  journal={medRxiv},
  pages={2026--06},
  year={2026},
  publisher={Cold Spring Harbor Laboratory Press}
}
```
