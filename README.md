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

---

Additional scripts for external validation, prognostic analyses, medication analyses, and the web-based subtype assignment tool will be added to this repository.
