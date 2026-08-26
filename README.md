# CardioMet

This repository contains the source code used in the study:

**Data-driven cardiometabolic phenotyping reveals distinct risk patterns in coronary artery disease**

The repository includes scripts for data preprocessing, feature selection, unsupervised clustering, subtype assignment, external validation, prognostic analyses, medication analyses, and figure generation.

---

## Repository contents

### Clustering analysis in the BIPASS cohort.R

This script performs the complete clustering workflow in the BIPASS discovery cohort, including:

- data preparation and quality control
- exclusion of clinically implausible outliers
- age- and sex-adjustment of candidate metabolic variables
- selection of the optimal number of clusters
- K-means clustering
- cluster stability assessment using bootstrap Jaccard similarity
- visualization of cardiometabolic subtypes (pie chart, UMAP, boxplots, radar plot, and pairwise comparison heatmap)
- development of the multinomial ridge regression classifier for subtype assignment in external cohorts

### External validation in the external cohorts.R

This script provides the common external-validation workflow applied to the IMPROVE, Mount Sinai Data Warehouse (MSDW), and UK Biobank (UKB) cohorts, including:

- exclusion of clinically implausible values
- age- and sex-adjustment and standardization of the prespecified clustering variables
- independent K-means clustering using the six-variable feature set and four-cluster solution derived from the BIPASS discovery framework
- assignment and characterization of the four cardiometabolic subtypes
- descriptive summaries and visualization of subtype distributions and metabolic profiles
- application of the multinomial ridge regression classifier developed in the BIPASS discovery cohort
- individual-level subtype assignment without refitting the BIPASS-derived classifier
- evaluation of classifier performance using one-versus-rest ROC analyses and subtype-specific AUC estimates with 95% confidence intervals

The same analytical workflow was applied across the three external validation cohorts after cohort-specific data preprocessing and harmonization. Cohort-specific raw-data extraction, variable construction, and harmonization are therefore not included in this common validation script.

---

Additional scripts for external validation, prognostic analyses, medication analyses, and the web-based subtype assignment tool will be added to this repository.
