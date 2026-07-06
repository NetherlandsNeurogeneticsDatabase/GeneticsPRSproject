# Neurogenomics Analysis Pipeline

Analysis pipeline accompanying the neurogenomics study integrating genotype, neuropathology, and clinical data from the Netherlands Neurogenomics Database (NND).

---

# 1. Data Integration

The project begins with harmonization of three complementary data sources:

- Genotype data
- Neuropathology reports
- Clinical phenotypes

---

# 2. Genomic Analysis

Genome-wide genotype data are processed through standard genetic analysis pipelines.

## Quality Control

Sample and variant quality control includes:


Outputs:

- Clean genotype dataset
- QC reports

---

## Principal Component Analysis (PCA)

Population structure is estimated using genome-wide variants.

Outputs:

- Principal components
- Population clustering

---

## Genome-Wide Association Study (GWAS)

Genome-wide association analyses identify variants associated with disease phenotypes.

Outputs:

- Manhattan plots
- QQ plots
- Summary statistics
- Significant loci

---

## Polygenic Risk Scores (PRS)

Polygenic risk scores are calculated using external GWAS summary statistics.

Workflow:

- PRS calculation
- Score normalization

Outputs:

- Individual PRS
- Distribution plots
- Case-control comparisons

---

## Locus Refinement

Associated genomic regions are further refined through:

- LD analysis
- Regional association plots

Outputs:

- Prioritized variants
- Candidate genes

---

# 3. LLM-Based Neuropathology Extraction

Neuropathology reports are converted into structured quantitative phenotypes using Large Language Models (LLMs).


Outputs:

- Structured pathology tables
- Quality-control reports
- Validation metrics

---

# 4. Integrative Analysis

The final stage combines genetics, pathology, and clinical phenotypes.

## PRS Enrichment

Evaluate enrichment of disease-specific PRS across neurodegenerative disorders.

---

## PRS–Clinical Associations

Association analyses between PRS and:

- diagnosis
- symptoms
- disease progression
- clinical outcomes

---

## PRS–Neuropathology Associations

Evaluate relationships between genetic risk and neuropathological hallmarks.

Examples:

- PRS vs Braak stage
- PRS vs CERAD score
- PRS vs regional pathology burden

---

## Locus Reallocation Assessment

Compare association signals across:

- Public GWAS
- NND GWAS
- Integrated analyses

to refine candidate loci and prioritize variants for follow-up.

---

# Software

## Genetics

- PLINK 2
- bcftools
- IMPUTE5
- VCFtools

## Python

- pandas
- NumPy
- SciPy
- scikit-learn
- statsmodels
- matplotlib


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

