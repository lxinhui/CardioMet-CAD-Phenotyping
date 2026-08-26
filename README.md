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

### Prognostic analysis across study cohorts.R

This script provides the common prognostic-analysis workflow applied across the BIPASS discovery cohort and the IMPROVE, Mount Sinai Data Warehouse (MSDW), and UK Biobank (UKB) validation cohorts, including:

- construction of competing-risk outcomes for major adverse cardiovascular events (MACE), with non-cardiovascular death treated as a competing event
- estimation of subtype-specific cumulative incidence of MACE using the Aalen–Johansen estimator
- calculation of MACE event counts and crude event proportions across cardiometabolic subtypes
- overall comparison of cumulative incidence functions using Gray's test
- pairwise Gray's tests between cardiometabolic subtypes with adjustment for multiple comparisons
- Fine–Gray competing-risk regression to estimate age- and sex-adjusted subdistribution hazard ratios for MACE
- estimation of follow-up duration using the reverse Kaplan–Meier method

The same prognostic-analysis framework was applied across the four study cohorts after cohort-specific outcome construction and data harmonization. Cohort-specific outcome definitions, follow-up periods, and variable construction may differ according to study design and data availability and are therefore handled before application of this common analytical workflow. The metabolically healthy (MH) subtype was used as the reference group in the Fine–Gray regression analyses.

