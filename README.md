# CardioMet

This repository contains the source code used in the study:

**Unsupervised cardiometabolic phenotyping reveals prognostic and treatment-associated heterogeneity in coronary artery disease**

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

### Medication analysis in the MSDW cohort.R

This script performs the medication analyses in the Mount Sinai Data Warehouse (MSDW) cohort to evaluate subtype-specific patterns of medication use and medication-associated differences in cardiometabolic biomarkers and major adverse cardiovascular events (MACE), including:

- construction of medication exposure variables for glucose-lowering, antihypertensive, and lipid-lowering medications using prespecified prescription windows around the index CAD diagnosis
- identification of prior medication use and concomitant medications within the same therapeutic class
- characterization and visualization of medication use patterns across the four cardiometabolic subtypes
- extraction of post-treatment cardiometabolic biomarker measurements using medication-specific lag periods and follow-up windows
- calculation of changes in relevant biomarkers from baseline and estimation of medication-associated biomarker changes using inverse probability of treatment weighting (IPTW)
- construction of medication-specific MACE outcomes and time-to-event variables using prespecified treatment-effect lag periods
- propensity score matching (PSM) of medication users and non-users within each cardiometabolic subtype, with adjustment for demographic characteristics, baseline cardiometabolic features, prior heart failure, prior medication use, and concomitant therapy as appropriate
- assessment of covariate balance after PSM using standardized mean differences
- estimation of medication-associated hazard ratios for MACE using Cox proportional hazards regression with robust standard errors in the overall cohort and within individual cardiometabolic subtypes
- visualization of cumulative MACE incidence among medication users across cardiometabolic subtypes
- comparison of medication-associated hazard ratios across cardiometabolic subtypes using ratios of hazard ratios (RHRs), including pairwise subtype comparisons with 95% confidence intervals and corresponding P values
- generation of figures and tabulated results for medication use patterns, biomarker changes, propensity-score balance, cumulative incidence, subtype-specific medication–MACE associations, and pairwise RHR analyses

Medication analyses were observational and were designed to characterize subtype-specific heterogeneity in medication-associated biomarker and cardiovascular outcomes rather than to establish causal treatment effects. Medication-specific eligibility criteria, exposure definitions, lag periods, covariates, and outcome construction were prespecified according to medication class and the availability of longitudinal electronic health record data.
