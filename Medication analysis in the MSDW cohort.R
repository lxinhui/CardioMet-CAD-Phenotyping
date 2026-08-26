#===============================================================================
# CAD Cardiometabolic Phenotyping and Medication Analysis
#===============================================================================
#
# Purpose:
# This script performs cardiometabolic phenotyping and evaluates
# subtype-specific associations of cardiometabolic medications with
# intermediate biomarkers and major adverse cardiovascular events (MACE)
# among patients with coronary artery disease (CAD).
#
# Main analytical workflow:
#
# 1. Data preparation and cohort construction
#    - Import CAD, metabolic biomarker, medication, and outcome datasets.
#    - Exclude patients with acute heart failure at the index CAD diagnosis.
#    - Exclude individuals with unknown sex.
#    - Link cardiometabolic measurements to the index CAD diagnosis using
#      prespecified measurement windows.
#    - Harmonize measurement units and calculate estimated glomerular
#      filtration rate (eGFR).
#    - Select one index CAD record per patient and handle missing values
#      using multiple imputation.
#
# 2. Cardiometabolic subtype identification
#    - Use six cardiometabolic features:
#        BMI, SBP, DBP, FBG, TG, and eGFR.
#    - Remove physiologically implausible values.
#    - Adjust each clustering feature for age and sex using linear regression.
#    - Standardize the resulting residuals.
#    - Perform K-means clustering (K = 4).
#    - Map clusters to four clinically interpretable CAD subtypes:
#        MH = metabolically healthy
#        HT = hypertensive
#        HG = hyperglycemic
#        RI = renal impairment
#
# 3. Medication exposure definition
#    - Evaluate glucose-lowering, antihypertensive, and lipid-lowering
#      medications.
#    - Define medication exposure around the index CAD diagnosis using
#      prespecified prescription windows.
#    - Construct indicators for prior medication use and concomitant
#      medications within the same therapeutic class.
#
# 4. Medication use patterns
#    - Describe medication use and non-use across the four CAD subtypes.
#    - Visualize subtype-specific medication utilization patterns.
#
# 5. Medication-associated changes in cardiometabolic biomarkers
#    - Extract post-treatment biomarker measurements within prespecified
#      drug-specific lag periods and follow-up windows.
#    - Calculate absolute changes from baseline in relevant biomarkers.
#    - Estimate medication-associated biomarker changes using inverse
#      probability of treatment weighting (IPTW).
#    - Perform analyses in the overall cohort and separately within each
#      cardiometabolic subtype.
#
# 6. Medication-associated MACE outcomes
#    - Define drug-specific follow-up windows accounting for treatment
#      effect-lag periods.
#    - Construct drug-specific MACE outcomes and time-to-event variables.
#    - Perform propensity score matching (PSM) within each CAD subtype.
#    - Assess post-matching covariate balance using standardized mean
#      differences (SMDs).
#    - Estimate medication-associated MACE hazard ratios using Cox
#      proportional hazards models with robust standard errors.
#
# 7. Subtype-specific heterogeneity in medication associations
#    - Compare medication-associated hazard ratios across CAD subtypes using
#      ratios of hazard ratios (RHRs).
#    - Calculate all pairwise subtype comparisons:
#
#          RHR = HR(comparison subtype) / HR(reference subtype)
#
#    - Estimate 95% confidence intervals and P values on the log-HR scale.
#    - Visualize pairwise RHRs using heatmaps to characterize heterogeneity
#      in medication-associated MACE risk across cardiometabolic subtypes.
#
# Main outputs:
#    - Medication use pattern figure
#    - IPTW estimates of medication-associated biomarker changes
#    - PSM covariate balance plots
#    - Cumulative incidence curves
#    - Subtype-specific medication-MACE HR estimates
#    - Pairwise RHR heatmap
#    - Pairwise RHR results table
#
# Notes:
#    - Medication analyses are observational and should be interpreted as
#      medication-associated rather than causal treatment effects.
#    - RHR calculations based on separately estimated subtype-specific HRs
#      assume independence of the corresponding log-HR estimates.
#
#===============================================================================

library(arrow)
library(dplyr)
library(data.table)
library(stringr)
library(readr)

ukbdata_meta <- read_parquet("cad_metabolic_combined_20260331_demo.parquet")
ukbdata_cad <- read_parquet("cad_cvd_outcomes_20260331_demo.parquet")
ukbdata_drug <- read_parquet("cad_drug_exposure_20260331_demo.parquet")
ukbdata_outcome <- read_parquet("cad_outcome_events_combined_20260331_demo.parquet")
ukbdata_outcome <- read_parquet("cad_outcome_events_combined_20260331_demo.parquet")

length(unique(ukbdata_cad$patient_id))
ukbdata_cad <- unique(ukbdata_cad)

#===============================================================
####-----Clustering Variable Selection and Preparation------####
#===============================================================

# 1. Exclude individuals with acute heart failure on the date of CAD diagnosis
ukbdata_cad <- ukbdata_cad %>%
  filter(
    is.na(HF_DATE) |
      abs(
        as.numeric(
          difftime(HF_DATE,
                   CAD_diagnosis_date,
                   units = "hours")
        )
      ) > 24
  )

# 2. Exclude individuals with unknown sex
unknown_gender_ids <- ukbdata_cad %>%
  group_by(patient_id) %>%
  summarise(
    all_not_male_female = all(!gender %in% c("Female", "Male")),
    .groups = "drop"
  ) %>%
  filter(all_not_male_female) %>%
  pull(patient_id)

ukbdata_cad <- subset(ukbdata_cad,!patient_id %in% unknown_gender_ids)
length(unique(ukbdata_cad$patient_id))

# Sex assignment for each participant
gender_by_person <- ukbdata_cad %>%
  group_by(patient_id) %>%
  summarise(
    gender_final = {
      g <- unique(gender[gender %in% c("Female", "Male")])
      if (length(g) == 1) {
        g
      } else {
        "Unknown"
      }
    },
    .groups = "drop"
  )

head(gender_by_person)

table(gender_by_person$gender_final)
table(is.na(gender_by_person$gender_final))
ukbdata_cad <- merge(ukbdata_cad,gender_by_person,all.x = T, by = "patient_id") 

# Review of biomarker measurement units
linshi <- subset(ukbdata_meta,metabolite_name == "BMI")
table(linshi$unit)
linshi <- subset(ukbdata_meta,metabolite_name == "SBP")
table(linshi$unit)
linshi <- subset(ukbdata_meta,metabolite_name == "DBP")
table(linshi$unit)
linshi <- subset(ukbdata_meta,metabolite_name == "LDL")#MG/DL
table(linshi$unit)
ukbdata_meta <- ukbdata_meta %>%
  mutate(
    value = case_when(
      metabolite_name == "LDL" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_real_,
      TRUE ~ value
    ),
    unit = case_when(
      metabolite_name == "LDL" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_character_,
      TRUE ~ unit
    )
  )

linshi <- subset(ukbdata_meta,metabolite_name == "TG")#MG/DL
table(linshi$unit)
ukbdata_meta <- ukbdata_meta %>%
  mutate(
    value = case_when(
      metabolite_name == "TG" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_real_,
      TRUE ~ value
    ),
    unit = case_when(
      metabolite_name == "TG" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_character_,
      TRUE ~ unit
    )
  )

linshi <- subset(ukbdata_meta,metabolite_name == "FBG")#MG/DL
table(linshi$unit)
ukbdata_meta <- ukbdata_meta %>%
  mutate(
    value = case_when(
      metabolite_name == "FBG" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_real_,
      TRUE ~ value
    ),
    unit = case_when(
      metabolite_name == "FBG" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL")) ~ NA_character_,
      TRUE ~ unit
    )
  )

linshi <- subset(ukbdata_meta,metabolite_name == "Scr")#MG/DL
table(linshi$unit)
ukbdata_meta <- ukbdata_meta %>%
  mutate(
    value = case_when(
      metabolite_name == "Scr" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL", "mg.dl")) ~ NA_real_,
      TRUE ~ value
    ),
    unit = case_when(
      metabolite_name == "Scr" & 
        !(unit %in% c("mg/dl", "mg/dL", "MG/DL", "mg.dl")) ~ NA_character_,
      TRUE ~ unit
    )
  )


# 2. Define the participant–CAD diagnosis pair as the unit of analysis.
# Cardiometabolic measurements were linked to the CAD diagnosis if obtained within:
# - ±30 days: FBG, TG, LDL, HDL, TC, Scr
# - ±90 days: SBP, DBP
# - ±365 days: BMI (or height and weight)
# Measurements within these windows were considered index-hospitalization values.

unique(ukbdata_meta$metabolite_name)

library(data.table)

setDT(ukbdata_cad)
setDT(ukbdata_meta)

ukbdata_cad[, CAD_diagnosis_date := as.POSIXct(CAD_diagnosis_date)]
ukbdata_meta[, measurement_date := as.POSIXct(measurement_date)]

window_list <- list(
  FBG = 30, TG = 30, LDL = 30, Scr = 30,
  SBP = 90, DBP = 90,
  BMI = 365
)

extract_nearest_measure <- function(meta_dt, cad_dt, metabolite, window_days) {
  
  meta_sub <- meta_dt[metabolite_name == metabolite]
  
  dt <- meta_sub[cad_dt, on = .(patient_id), allow.cartesian = TRUE]
  
  dt[, diff_days := as.numeric(difftime(measurement_date, CAD_diagnosis_date, units = "days"))]
  
  dt <- dt[abs(diff_days) <= window_days]
  
  dt <- dt[order(abs(diff_days))]
  dt <- dt[, .SD[1], by = .(patient_id, CAD_diagnosis_date)]

  dt <- dt[, .(patient_id, CAD_diagnosis_date, value)]
  setnames(dt, "value", metabolite)
  
  return(dt)
}

result_list <- list()

for (met in names(window_list)) {
  cat("Processing:", met, "\n")
  
  result_list[[met]] <- extract_nearest_measure(
    ukbdata_meta,
    ukbdata_cad,
    metabolite = met,
    window_days = window_list[[met]]
  )
}

final_dt <- copy(ukbdata_cad)

for (met in names(result_list)) {
  final_dt <- merge(
    final_dt,
    result_list[[met]],
    by = c("patient_id", "CAD_diagnosis_date"),
    all.x = TRUE
  )
}

ukbdata_cad <- final_dt
colnames(ukbdata_cad)

# Common unit conversion factors (mg/dL → mmol/L)
ukbdata_cad <- ukbdata_cad %>%
  mutate(
    LDL  = LDL * 0.02586,
    TG   = TG  * 0.01129,
    FBG  = FBG * 0.0555
  )

# eGFR calculation
ukbdata_cad$k1 <- ifelse(ukbdata_cad$gender_final == "Female",0.7,0.9)
ukbdata_cad$a1 <- ifelse(ukbdata_cad$gender_final == "Female",-0.329,-0.411)
ukbdata_cad$b <- ifelse(ukbdata_cad$gender_final == "Female",1.018,1)
ukbdata_cad$sex <- ifelse(ukbdata_cad$gender_final == "Female", 0, 1)

ukbdata_cad <- ukbdata_cad %>%
  mutate(
    eGFR = apply(., 1, function(row) {
      sex <- as.numeric(row["sex"])
      age <- as.numeric(row["age_at_cad"])
      Scr <- as.numeric(row["Scr"])
      k1 <- as.numeric(row["k1"])
      a1 <- as.numeric(row["a1"])
      b <- as.numeric(row["b"])
      
      141 * ((min(Scr/k1,1))^a1) * ((max(Scr/k1,1))^(-1.209)) * (0.993^age) * b
      
    })
  )

# Exclude observations with more than four missing values among the six clustering variables
cluster_vars <- c("BMI", "DBP", "FBG", "SBP", "TG", "eGFR")

ukbdata_clean <- ukbdata_cad %>%
  mutate(
    miss_n = rowSums(is.na(select(., all_of(cluster_vars))))
  ) %>%
  filter(miss_n <= 2) %>%   
  select(-miss_n)         

length(unique(ukbdata_cad$patient_id))
length(unique(ukbdata_clean$patient_id))

# Rule 1 (primary):
# Within each patient_id, retain the record with the fewest missing values
# among the six key clustering variables.

# Rule 2 (secondary):
# If multiple records have the same minimum number of missing values,
# retain the record with the earliest CAD diagnosis date.
vars8 <-  c("BMI", "DBP", "FBG", "SBP", "TG", "eGFR")

ukbdata_final <- ukbdata_clean %>%

  rowwise() %>%
  mutate(
    missing_n = sum(is.na(c_across(all_of(vars8))))
  ) %>%
  ungroup() %>%
  
  arrange(
    patient_id,
    missing_n,                   
    CAD_diagnosis_date      
  ) %>%
  
  group_by(patient_id) %>%
  slice(1) %>%
  ungroup() %>%
  
  select(-missing_n)

# Missingness assessment of variables
missing_df <- ukbdata_final %>%
  summarise(
    across(
      all_of(vars8),
      list(
        missing_n = ~ sum(is.na(.)),
        missing_prop = ~ mean(is.na(.))
      ),
      .names = "{.col}_{.fn}"
    )
  ) %>%
  tidyr::pivot_longer(
    cols = everything(),
    names_to = c("variable", ".value"),
    names_pattern = "^(.*)_(missing_n|missing_prop)$"
  ) %>%
  arrange(desc(missing_prop))

missing_df

# Multiple imputation of missing values
library(mice)

vars_imp <- c(
  "BMI", "DBP", "FBG",
  "LDL", "SBP", "TG",
  "eGFR", "age_at_cad"
)

data_imp <- ukbdata_final[, vars_imp]
meth <- make.method(data_imp)
meth[vars_imp] <- "pmm" 
pred <- make.predictorMatrix(data_imp)
pred[,] <- 1
diag(pred) <- 0
set.seed(20260126)

imp <- mice(
  data_imp,
  m = 10,
  method = meth,
  predictorMatrix = pred,
  maxit = 10,
  printFlag = TRUE
)

imp_list <- lapply(1:10, function(i) complete(imp, i))
imp_mean <- Reduce(
  "+",
  imp_list
) / length(imp_list)

ukbdata_final[, vars_imp] <- imp_mean

# Convert sex coding in the new dataset from 0 = female / 1 = male
# to 1 = male / 2 = female to match the training dataset
ukbdata_final$sex <- ifelse(
  ukbdata_final$sex == 1, 1,
  ifelse(ukbdata_final$sex == 0, 2, NA)
)

datafixfix <- ukbdata_final
data_mean_imputed <- ukbdata_final
ukbdata <- ukbdata_final

#==========================================================================
###############-------------Clustering analysis--------------###############
#==========================================================================
library(cluster)
library(factoextra)
library(fpc)
library(tidyverse)

vars <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI", "age_at_cad","sex")#5类
cluster_vars_numeric <- vars

data_cluster <- data_mean_imputed[, c("patient_id",cluster_vars_numeric)]
data_cluster <- na.omit(data_cluster)


id_col <- data_cluster["patient_id"]  
data_values <- data_cluster[, cluster_vars_numeric] 
abnormal_limits <- list(
  FBG  = c(low = 2,  high = 35),
  SBP  = c(low = 60,   high = 280),
  DBP  = c(low = 30,   high = 150),
  eGFR = c(low = 1,    high = 180),
  TC   = c(low = 1,    high = 15),
  TG   = c(low = 0.1,  high = 40),
  LDL  = c(low = 0.2,    high = 12), 
  HDL  = c(low = 0.2,  high = 5), 
  BMI  = c(low = 12,   high = 70),
  age_at_cad  = c(low = 15,   high = 85)
)

detect_medical_outlier <- function(df, limits_list) {
  
  limits_use <- limits_list[names(limits_list) %in% names(df)]
  
  limits_use <- limits_use[match(names(df), names(limits_use))]
  
  extreme_mat <- mapply(function(x, limits) {
    if (is.null(limits)) {
      return(rep(FALSE, length(x)))
    }
    x < limits["low"] | x > limits["high"]
  }, df, limits_use)
  
 
  extreme_rows <- apply(extreme_mat, 1, any)
  return(extreme_rows)
}

extreme_rows <- detect_medical_outlier(data_values, abnormal_limits)
data_values_clean <- data_values[!extreme_rows, ]
id_col_clean <- id_col[!extreme_rows, , drop = FALSE]

# for each metabolic biomarker: regression---Biomarker = β0 + β1 × Age + β2 × Gender + ε
# ε: Centered and scaled to mean = 0, SD = 1
data_values_clean$sex <- as.factor(data_values_clean$sex)

cluster_vars <- setdiff(colnames(data_values_clean), c("age_at_cad", "sex"))

data_scaled <- data.frame(row.names = rownames(data_values_clean))

for (v in cluster_vars) {
  
  fit <- lm(
    formula = as.formula(paste(v, "~ age_at_cad + sex")),
    data = data_values_clean
  )
  
  res <- resid(fit)
  
  res_scaled <- as.numeric(scale(res))
  
  data_scaled[[v]] <- res_scaled
}

str(data_scaled)

data_cluster_scaled <- cbind(id = id_col_clean, as.data.frame(data_scaled))

ukbdata <- subset(ukbdata_final,patient_id %in% data_cluster_scaled$patient_id)
datafixfix <- subset(datafixfix,patient_id %in% data_cluster_scaled$patient_id)#####
data_cluster <- subset(data_cluster,patient_id %in% data_cluster_scaled$patient_id)#####

# statistical description
library(dplyr)
library(gtsummary)

vars1 <- c(
  "age_at_cad", "sex", "BMI", "SBP", "DBP",
  "TG", "FBG", "eGFR"
)

table1 <- data_cluster %>%
  select(all_of(vars1)) %>%
  tbl_summary(
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    )
  ) %>%
  add_n() %>% 
  modify_header(label ~ "**Variables**") %>%
  bold_labels()

table1

# Clustering analysis
set.seed(20250814)

km_res <- kmeans(data_scaled, centers = 4, nstart = 25)

data_cluster$cluster <- factor(km_res$cluster)
datafixfix$clusterpamage <- factor(km_res$cluster)
ukbdata$clusterpamage <- factor(km_res$cluster)

recode_cluster <- function(x) {
  dplyr::recode(x,
                `1` = 3,
                `2` = 4,
                `3` = 1,
                `4` = 2)
}

data_cluster$cluster        <- recode_cluster(data_cluster$cluster)
datafixfix$clusterpamage    <- recode_cluster(datafixfix$clusterpamage)
ukbdata$clusterpamage       <- recode_cluster(ukbdata$clusterpamage)


colnames(data_cluster)[1] <- "id"
colnames(ukbdata)[1] <- "id"

#==========================================================
######-------------Medication Analysis---------------######
#==========================================================
# Restrict the medication analysis dataset to individuals with assigned cluster labels
ukbdata_drug <- subset(ukbdata_drug,patient_id %in% data_cluster$id)
colnames(ukbdata_drug)[1] <- "id"

# Restrict the metabolic dataset to individuals with assigned cluster labels
ukbdata_meta <- subset(ukbdata_meta,patient_id %in% data_cluster$id)
colnames(ukbdata_meta)[1] <- "id"

ukbdata_drug <- ukbdata_drug %>%
  mutate(
    drug_category = case_when(
      drug_category %in% c("ACEI", "ARBs") ~ "ACEIARBs",
      TRUE ~ drug_category
    )
  )

table(ukbdata_drug$drug_category)

# Definition of major categories
drug_superclass <- list(
  
  glucose_lowering = c(
    "metformin",
    "sulfonylureas",
    "SGLT2"
  ),
  
  antihypertensive = c(
    "ACEIARBs",
    "beta_blocker",
    "CCB"
  ),
  
  lipid_lowering = "statin"
  
)

# 1. Create binary indicators for medication use within each medication category
create_drug_exposure_by_category <- function(
    ukbdata,
    ukbdata_drug,
    drug_superclass,
    window_before = 60,
    window_after  = 30,
    prior_start   = 365,
    prior_end     = 60
){
  
  library(dplyr)
  
  ukbdata <- ukbdata %>%
    mutate(CAD_diagnosis_date = as.Date(CAD_diagnosis_date))
  
  ukbdata_drug <- ukbdata_drug %>%
    mutate(prescription_date = as.Date(prescription_date))
  
  drug_with_cad <- ukbdata_drug %>%
    left_join(
      ukbdata %>% select(id, CAD_diagnosis_date),
      by = "id"
    ) %>%
    mutate(
      time_diff = as.numeric(prescription_date - CAD_diagnosis_date)
    )
  
  drug_window <- drug_with_cad %>%
    filter(!is.na(time_diff),
           time_diff >= -window_before,
           time_diff <= window_after)
  
  drug_prior <- drug_with_cad %>%
    filter(!is.na(time_diff),
           time_diff >= -prior_start,
           time_diff < -prior_end)
  
  result <- ukbdata %>%
    select(id)
  
  all_drugs <- unlist(drug_superclass)
  
  for (drug in all_drugs) {
    
    superclass <- names(drug_superclass)[
      sapply(drug_superclass, function(x) drug %in% x)
    ]
    
    superclass_drugs <- drug_superclass[[superclass]]
    
    # ----------------------------
    # exposure
    # ----------------------------
    
    exposed_df <- drug_window %>%
      filter(drug_category == drug) %>%
      group_by(id) %>%
      summarise(
        drug_start_date = min(prescription_date),
        .groups = "drop"
      )
    
    exposed_ids <- exposed_df$id
    
    # ----------------------------
    # Other medications within the same class (excluding itself)
    # ----------------------------
    
    concomitant_ids <- drug_window %>%
      filter(drug_category %in% superclass_drugs,
             drug_category != drug) %>%
      pull(id) %>%
      unique()
    
    # ----------------------------
    # clean control
    # ----------------------------
    
    window_flags <- drug_window %>%
      group_by(id) %>%
      summarise(
        has_target = any(drug_category == drug),
        has_super  = any(drug_category %in% superclass_drugs),
        .groups = "drop"
      )
    
    excluded_ids <- window_flags %>%
      filter(!has_target & has_super) %>%
      pull(id)
    
    all_ids <- result$id
    
    clean_control <- setdiff(all_ids,
                             union(exposed_ids, excluded_ids))
    
    # ----------------------------
    # exposure
    # ----------------------------
    
    result[[drug]] <- case_when(
      result$id %in% exposed_ids ~ 1,
      result$id %in% clean_control ~ 0,
      TRUE ~ NA_real_
    )
    
    # ----------------------------
    # concomitant
    # ----------------------------
    
    result[[paste0(drug,"_concomitant")]] <- case_when(
      result$id %in% exposed_ids &
        result$id %in% concomitant_ids ~ 1,
      TRUE ~ 0
    )
    
    # ----------------------------
    # start time
    # ----------------------------
    
    result <- result %>%
      left_join(
        exposed_df %>%
          rename(!!paste0(drug,"_start_date") := drug_start_date),
        by = "id"
      )
    
    # ----------------------------
    # prior
    # ----------------------------
    
    prior_ids <- drug_prior %>%
      filter(drug_category == drug) %>%
      pull(id) %>%
      unique()
    
    result[[paste0("prior_", drug)]] <- ifelse(
      result$id %in% prior_ids, 1, 0
    )
    
  }
  
  return(result)
}

result <- create_drug_exposure_by_category(
  ukbdata = ukbdata,
  ukbdata_drug = ukbdata_drug,
  drug_superclass = drug_superclass
)

ukbdata <- ukbdata %>%
  left_join(result, by = "id")

data_cluster <- data_cluster %>%
  left_join(
    result,
    by = "id"
  )

#========================================================
# Count the number of users and non-users for each medication
#========================================================

library(dplyr)
library(tidyr)

drug_vars <- c(
  "metformin",
  "sulfonylureas",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

df_long <- data_cluster %>%
  select(cluster, all_of(drug_vars)) %>%
  pivot_longer(cols = all_of(drug_vars),
               names_to = "drug",
               values_to = "use")

table_drug <- df_long %>%
  
  group_by(cluster, drug) %>%
  
  summarise(
    
    N_use = sum(use == 1, na.rm = TRUE),
    
    N_non = sum(use == 0, na.rm = TRUE),
    
    N_total = N_use + N_non, 
    
    pct_use = round(N_use / N_total * 100, 1),
    
    .groups = "drop"
    
  )

table_drug <- table_drug %>%
  mutate(
    Use = paste0(N_use, " (", pct_use, "%)"),
    Non_use = paste0(N_non, " (", 100 - pct_use, "%)")
  )

library(ggplot2)
library(dplyr)

cluster_map <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)

custom_colors <- c(
  "MH" = "#4E79A7",
  "HG" = "#E15759",
  "HT" = "#B07AA1",
  "RI" = "#F28E2B"
)

keep_drugs <- c(
  "ACEIARBs",
  "CCB",
  "beta_blocker",
  "metformin",
  "SGLT2",
  "statin"
)

plot_df <- df_long %>%
  filter(drug %in% keep_drugs) %>%
  group_by(cluster, drug) %>%
  summarise(
    N_use = sum(use == 1, na.rm = TRUE),
    N_non = sum(use == 0, na.rm = TRUE),
    N_total = N_use + N_non,
    pct = if_else(
      N_total > 0,
      N_use / N_total * 100,
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  
  mutate(
    
    cluster = factor(
      cluster_map[as.character(cluster)],
      levels = c("MH", "HT", "HG", "RI")
    ),
    
    drug = factor(
      drug,
      levels = c(
        "CCB",
        "ACEIARBs",
        "beta_blocker",
        "SGLT2",
        "metformin",
        "statin"
      )
    ),

    label = paste0(
      N_use,
      "\n(",
      sprintf("%.1f", pct),
      "%)"
    )
  )

p <- ggplot(
  plot_df,
  aes(
    x = drug,
    y = pct,
    fill = cluster
  )
) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.72),
    vjust = -0.25,
    size = 2.9,
    lineheight = 0.9
  ) +
  scale_fill_manual(values = custom_colors) +
  scale_y_continuous(
    limits = c(0, max(plot_df$pct) + 14),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Medication",
    y = "Use (%)",
    fill = "Cluster"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    ),
    axis.title = element_text(
      face = "bold",
      size = 13
    ),
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(
      size = 11
    ),
    panel.grid = element_blank()
  )

print(p)

ggsave(
  filename = "Medication_use_patterns.pdf",
  width = 12,
  height = 7,
  device = cairo_pdf  
)

#========================================================
# Generate follow-up biomarker levels
#========================================================
extract_postdrug_biomarker_from_ukb <- function(
    ukbdata,
    ukbmeta,
    drug_vars,
    drug_lag_days,
    drug_biomarker_map,
    start_date_suffix = "_start_date",
    followup_days = 365
){
  library(dplyr)
  library(lubridate)
  library(tidyr)

  ukbdata <- ukbdata %>%
    mutate(
      CAD_diagnosis_date = as.Date(CAD_diagnosis_date)
    )
  
  ukbmeta <- ukbmeta %>%
    mutate(
      measurement_date = as.Date(measurement_date)  
    )
  
  result <- ukbdata %>%
    select(id)
  
  for(drug in drug_vars){
    
    cat("Processing:", drug, "\n")
    
    lag        <- drug_lag_days[[drug]]
    biomarkers <- drug_biomarker_map[[drug]]
    start_col  <- paste0(drug, start_date_suffix)
    
    # 校验
    if(!(drug %in% names(ukbdata))){
      cat("Missing variables:", drug, "\n")
      next
    }
    if(!(start_col %in% names(ukbdata))){
      cat("Missing start date:", start_col, "\n")
      next
    }
    
    tmp <- ukbdata %>%
      select(id, CAD_diagnosis_date, all_of(drug), all_of(start_col)) %>%
      mutate(
        drug_start_date = as.Date(.data[[start_col]]),
        exposure = .data[[drug]],
        
        index_date = case_when(
          exposure == 1 ~ drug_start_date,
          exposure == 0 ~ CAD_diagnosis_date,
          TRUE ~ as.Date(NA)
        ),
        
        window_start = index_date + days(lag),
        window_end   = index_date + days(followup_days)
      ) %>%
      filter(exposure %in% c(0,1) & !is.na(index_date))
    
    if(nrow(tmp) == 0){
      cat("No available samples:", drug, "\n")
      next
    }
    
    biomarker_window <- ukbmeta %>%
      filter(metabolite_name %in% biomarkers) %>%  
      inner_join(
        tmp %>% select(id, window_start, window_end),
        by = "id"
      ) %>%
      filter(
        measurement_date >= window_start &        
          measurement_date <= window_end
      ) %>%
      arrange(id, metabolite_name, measurement_date)  
    
    if(nrow(biomarker_window) == 0){
      cat("No biomarker data available:", drug, "\n")
      next
    }
    
    biomarker_summary <- biomarker_window %>%
      group_by(id, metabolite_name) %>%       
      
      # Use the first available measurement
      slice(1) %>% 
      ungroup() %>%
      select(id, metabolite_name, value) %>%
      pivot_wider(
        names_from  = metabolite_name,          
        values_from = value
      )
    
    new_names <- paste0(drug, "_", colnames(biomarker_summary)[-1])
    colnames(biomarker_summary)[-1] <- new_names
    
    result <- result %>%
      left_join(biomarker_summary, by = "id")
  }
  
  return(result)
}

drug_vars <- c(
  "metformin",
  "sulfonylureas",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

drug_biomarker_map <- list(
  metformin      = c("FBG","SBP","DBP","TG","BMI","Scr"),
  sulfonylureas  = c("FBG","SBP","DBP","TG","BMI","Scr"),
  SGLT2          = c("FBG","SBP","DBP","TG","BMI","Scr"),
  ACEIARBs       = c("FBG","SBP","DBP","TG","BMI","Scr"),
  beta_blocker   = c("FBG","SBP","DBP","TG","BMI","Scr"),
  CCB            = c("FBG","SBP","DBP","TG","BMI","Scr"),
  statin         = c("FBG","SBP","DBP","TG","BMI","Scr")
)

drug_lag_days <- list(
  metformin      = 90,
  sulfonylureas  = 90,
  SGLT2          = 90,
  ACEIARBs       = 60,
  beta_blocker   = 60,
  CCB            = 60,
  statin         = 90
)

drug_followup_biomarker <- extract_postdrug_biomarker_from_ukb(
  ukbdata = ukbdata,
  ukbmeta = ukbdata_meta,
  drug_vars = drug_vars,
  drug_lag_days = drug_lag_days,
  drug_biomarker_map = drug_biomarker_map,
  start_date_suffix = "_start_date",
  followup_days = 365
)

ukbdata <- ukbdata %>%
  left_join(drug_followup_biomarker, by = "id")

library(dplyr)
library(stringr)

scr_vars <- names(ukbdata)[grepl("_Scr$", names(ukbdata))]
scr_vars

ukbdata <- ukbdata %>%
  mutate(
    k1 = ifelse(sex == "2", 0.7, 0.9),
    a1 = ifelse(sex == "2", -0.329, -0.411),
    b  = ifelse(sex == "2", 1.018, 1)
  ) %>%
  mutate(
    across(
      ends_with("_Scr"),
      ~ 141 * (pmin(.x / k1, 1)^a1) *
        (pmax(.x / k1, 1)^(-1.209)) *
        (0.993^age_at_cad) * b,
      .names = "{str_remove(.col, '_Scr$')}_eGFR"
    )
  )


# Convert units using standard conversion factors (mg/dL to mmol/L)
ukbdata <- ukbdata %>%
  mutate(
    metformin_TG   = metformin_TG  * 0.01129,
    sulfonylureas_TG   = sulfonylureas_TG  * 0.01129,
    SGLT2_TG   = SGLT2_TG  * 0.01129,
    ACEIARBs_TG   = ACEIARBs_TG  * 0.01129,
    beta_blocker_TG   = beta_blocker_TG  * 0.01129,
    CCB_TG   = CCB_TG  * 0.01129,
    statin_TG   = statin_TG  * 0.01129,
    
    metformin_FBG  = metformin_FBG * 0.0555,
    sulfonylureas_FBG  = sulfonylureas_FBG * 0.0555,
    SGLT2_FBG  = SGLT2_FBG * 0.0555,
    ACEIARBs_FBG  = ACEIARBs_FBG * 0.0555,
    beta_blocker_FBG  = beta_blocker_FBG * 0.0555,
    CCB_FBG  = CCB_FBG * 0.0555,
    statin_FBG  = statin_FBG * 0.0555
  )

#========================================================
# Generate History of Heart Failure Variable
#========================================================

colnames(ukbdata_outcome)[1] <- "id"

library(dplyr)
library(lubridate)

ukbdata <- ukbdata %>%
  mutate(
    CAD_diagnosis_date = as.POSIXct(CAD_diagnosis_date)
  )

ukbdata_outcome <- ukbdata_outcome %>%
  mutate(
    event_date = as.POSIXct(event_date)
  )

hf_data <- ukbdata_outcome %>%
  filter(event_type == "HF") %>%
  select(id, event_date)

drug_vars <- c(
  "metformin",
  "sulfonylureas",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

for (drug in drug_vars) {
  
  cat("Processing:", drug, "\n")
  
  start_col <- paste0(drug, "_start_date")
  
  if (!(drug %in% colnames(ukbdata))) {
    cat("missing variable:", drug, "\n")
    next
  }
  
  if (!(start_col %in% colnames(ukbdata))) {
    cat("missing start time:", start_col, "\n")
    next
  }
  
  ukbdata[[start_col]] <- as.POSIXct(ukbdata[[start_col]])
  
  tmp_index <- ukbdata %>%
    select(id, CAD_diagnosis_date, all_of(drug), all_of(start_col)) %>%
    mutate(
      drug_start_date = .data[[start_col]],
      exposure = .data[[drug]],
      
      index_date = case_when(
        exposure == 1 ~ drug_start_date,
        exposure == 0 ~ CAD_diagnosis_date,
        TRUE ~ as.POSIXct(NA)
      )
    ) %>%
    filter(!is.na(index_date)) %>%
    select(id, index_date)
  
  prior_hf_tmp <- hf_data %>%
    inner_join(tmp_index, by = "id") %>%
    filter(event_date < index_date) %>%
    distinct(id) %>%
    mutate(prior_HF = 1)
  
  new_var <- paste0("prior_HF_", drug)
  
  ukbdata <- ukbdata %>%
    left_join(prior_hf_tmp %>% rename(!!new_var := prior_HF), by = "id") %>%
    mutate(
      !!new_var := ifelse(is.na(.data[[new_var]]), 0, .data[[new_var]])
    )
}

#========================================================
# Calculate absolute changes
#========================================================

drug_biomarker_map <- list(
  metformin      = c("FBG"),
  sulfonylureas  = c("FBG"),
  SGLT2          = c("FBG","SBP","DBP","TG","BMI","eGFR"),
  ACEIARBs       = c("SBP","DBP"),
  beta_blocker   = c("SBP","DBP"),
  CCB            = c("SBP","DBP"),
  statin         = c("TG")
)

create_delta_biomarkers <- function(
    data,
    drug_biomarker_map
){
  result <- data
  
  for(drug in names(drug_biomarker_map)){
    
    biomarkers <- drug_biomarker_map[[drug]]
    
    for(bio in biomarkers){
      
      post_var <- paste0(drug, "_", bio)
      base_var <- bio
      delta_var <- paste0(drug, "_", bio, "_delta")
      
      if(!(post_var %in% colnames(data))) next
      if(!(base_var %in% colnames(data))) next
      
      result[[delta_var]] <- result[[post_var]] - result[[base_var]]
      
      cat("✔ Created:", delta_var, "\n")
    }
  }
  
  return(result)
}

ukbdata <- create_delta_biomarkers(
  data = ukbdata,
  drug_biomarker_map = drug_biomarker_map
)

#========================================================
# Estimate Medication Effects Using IPTW
#========================================================

covariate_map <- list(
  
  # Glucose-lowering medications → FBG
  
  metformin_FBG = c(
    "age_at_cad",
    "sex",
    "FBG",
    "BMI",
    "eGFR",
    "metformin_concomitant",
    "prior_metformin",
    "prior_sulfonylureas",
    "prior_SGLT2"
  ),
  
  SGLT2_FBG = c(
    "age_at_cad",
    "sex",
    "FBG",
    "BMI",
    "eGFR",
    "SGLT2_concomitant",
    "prior_metformin",
    "prior_sulfonylureas",
    "prior_SGLT2"
  ),
  
  # SGLT2i → Blood pressure
  
  SGLT2_SBP = c(
    "age_at_cad",
    "sex",
    "SBP",
    "FBG",
    "BMI",
    "eGFR",
    "SGLT2_concomitant",
    "prior_metformin",
    "prior_sulfonylureas",
    "prior_SGLT2"
  ),
  
  SGLT2_DBP = c(
    "age_at_cad",
    "sex",
    "DBP",
    "FBG",
    "BMI",
    "eGFR",
    "SGLT2_concomitant",
    "prior_metformin",
    "prior_sulfonylureas",
    "prior_SGLT2"
  ),
  
  # SGLT2i → Kidney function
  SGLT2_eGFR = c(
    "age_at_cad",
    "sex",
    "eGFR",
    "FBG",
    "BMI",
    "SBP",
    "SGLT2_concomitant",
    "prior_metformin",
    "prior_sulfonylureas",
    "prior_SGLT2"
  ),
  
  # ACEI/ARB → Blood pressure
  
  ACEIARBs_SBP = c(
    "age_at_cad",
    "sex",
    "SBP",
    "BMI",
    "eGFR",
    "ACEIARBs_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  ),
  
  ACEIARBs_DBP = c(
    "age_at_cad",
    "sex",
    "DBP",
    "BMI",
    "eGFR",
    "ACEIARBs_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  ),
  
  # beta-blocker → Blood pressure
  
  beta_blocker_SBP = c(
    "age_at_cad",
    "sex",
    "SBP",
    "BMI",
    "eGFR",
    "beta_blocker_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  ),
  
  beta_blocker_DBP = c(
    "age_at_cad",
    "sex",
    "DBP",
    "BMI",
    "eGFR",
    "beta_blocker_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  ),
  
  # CCB → Blood pressure
  
  CCB_SBP = c(
    "age_at_cad",
    "sex",
    "SBP",
    "BMI",
    "eGFR",
    "CCB_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  ),
  
  CCB_DBP = c(
    "age_at_cad",
    "sex",
    "DBP",
    "BMI",
    "eGFR",
    "CCB_concomitant",
    "prior_CCB",
    "prior_beta_blocker",
    "prior_ACEIARBs"
  )
)

# IPTW function
run_iptw_lm <- function(
    data,
    drug_name,
    biomarker,
    subgroup_name = "Overall"
){
  
  cat("\n---", drug_name, "|", biomarker, "|", subgroup_name, "---\n")
  
  #======================================================
  # outcome
  #======================================================
  
  outcome_var <- paste0(
    drug_name, "_", biomarker, "_delta"
  )
  
  #======================================================
  # covariates
  #======================================================
  
  covar_key <- paste0(
    drug_name, "_", biomarker
  )
  
  ps_covars <- covariate_map[[covar_key]]
  
  if(is.null(ps_covars)){
    
    cat("Undefined covariates", covar_key, "\n")
    
    return(NULL)
  }

  ps_covars <- unique(ps_covars)
  
  vars_needed <- unique(c(
    drug_name,
    outcome_var,
    ps_covars
  ))
  
  missing_vars <- vars_needed[
    !vars_needed %in% colnames(data)
  ]
  
  if(length(missing_vars) > 0){
    
    cat(
      "missing variable",
      paste(missing_vars, collapse=", "),
      "\n"
    )
    
    return(NULL)
  }
  
  d <- data %>%
    
    select(all_of(vars_needed)) %>%
    
    rename(
      drug = all_of(drug_name),
      outcome = all_of(outcome_var)
    ) %>%
    
    filter(
      !is.na(drug),
      !is.na(outcome)
    )
  
  if(nrow(d) < 30) return(NULL)
  
  if(dplyr::n_distinct(d$drug) < 2)
    return(NULL)
  
  n_treat <- sum(d$drug == 1)
  
  n_ctrl <- sum(d$drug == 0)
  
  if(n_treat < 10 | n_ctrl < 10)
    return(NULL)
  
  #======================================================
  # PS model
  #======================================================
  
  formula_ps <- as.formula(
    paste(
      "drug ~",
      paste(ps_covars, collapse = " + ")
    )
  )
  
  cat("\nPS model:\n")
  
  print(formula_ps)
  
  ps_model <- glm(
    formula_ps,
    data = d,
    family = binomial
  )
  
  d$ps <- predict(
    ps_model,
    type = "response"
  )
  
  #======================================================
  # IPTW weight（ATE）
  #======================================================
  
  d$w <- ifelse(
    d$drug == 1,
    1 / d$ps,
    1 / (1 - d$ps)
  )
  
  w_cap <- quantile(
    d$w,
    0.99,
    na.rm = TRUE
  )
  
  d$w <- pmin(d$w, w_cap)
  
  #======================================================
  # weighted linear regression
  #======================================================
  
  fit <- lm(
    outcome ~ drug,
    data = d,
    weights = w
  )
  
  res <- broom::tidy(
    fit,
    conf.int = TRUE
  ) %>%
    
    filter(term == "drug") %>%
    
    mutate(
      drug = drug_name,
      biomarker = biomarker,
      clusterpamage = as.character(subgroup_name),
      n_total = nrow(d),
      n_treat = n_treat,
      n_ctrl = n_ctrl,
      n_weighted = sum(d$w)
    ) %>%
    
    select(
      drug,
      biomarker,
      clusterpamage,
      estimate,
      conf.low,
      conf.high,
      p.value,
      n_total,
      n_weighted,
      n_treat,
      n_ctrl
    )
  
  return(res)
}


pairs_df <- tibble::tibble(
  drug = c(
    "metformin","SGLT2",
    "SGLT2","SGLT2","SGLT2",
    "ACEIARBs","ACEIARBs",
    "beta_blocker","beta_blocker",
    "CCB","CCB"
  ),
  biomarker = c(
    "FBG","FBG",
    "SBP","DBP","eGFR",
    "SBP","DBP",
    "SBP","DBP",
    "SBP","DBP"
  )
)

drug_superclass_map <- list(
  glucose_lowering = c("metformin","sulfonylureas","SGLT2"),
  antihypertensive = c("ACEIARBs","beta_blocker","CCB"),
  lipid_lowering   = c("statin")
)

results_list <- list()

for(i in seq_len(nrow(pairs_df))){
  
  drug      <- pairs_df$drug[i]
  biomarker <- pairs_df$biomarker[i]
  
  cat("\n========================\n")
  cat("Processing:", drug, "-", biomarker, "\n")
  
  # ===== cluster =====
  cluster_res <- lapply(unique(ukbdata$clusterpamage), function(cl){
    
    dat_sub <- ukbdata %>% filter(clusterpamage == cl)
    
    run_iptw_lm(
      data = dat_sub,
      drug_name = drug,
      biomarker = biomarker,
      subgroup_name = cl
    )
  })
  
  cluster_res <- bind_rows(cluster_res)
  
  # ===== overall =====
  overall_res <- run_iptw_lm(
    data = ukbdata,
    drug_name = drug,
    biomarker = biomarker,
    subgroup_name = "Overall"
  )
  
  results_list[[paste0(drug,"_",biomarker)]] <- 
    bind_rows(cluster_res, overall_res)
}

results <- bind_rows(results_list)

write.csv(
  results,
  file = "drug_biomarkerchange_results-main.csv",
  row.names = FALSE
)

# See the separate script for plotting



# 5. Estimate the effects of medications on clinical outcomes
colnames(ukbdata_outcome)[1] <- "id"

library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(survival)
library(broom)
library(ggplot2)

drug_vars <- c(
  "metformin",
  "sulfonylureas",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

drug_lag_days <- list(
  metformin      = 90,
  sulfonylureas  = 90,
  SGLT2          = 90,
  ACEIARBs       = 60,
  beta_blocker   = 60,
  CCB            = 60,
  statin         = 90
)

# Define baseline MACE status (0/1)
library(dplyr)
library(lubridate)

mace_types <- c("HF", "IS", "MI", "CVD_Death")

ukbdata <- ukbdata %>%
  mutate(
    CAD_diagnosis_date = as.Date(CAD_diagnosis_date)
  )

ukbdata_outcome <- ukbdata_outcome %>%
  mutate(
    event_date = as.Date(event_date)
  )

# 3. Construct the baseline time window
baseline_mace_ids <- ukbdata_outcome %>%
  filter(event_type %in% mace_types) %>%
  
  inner_join(
    ukbdata %>% select(id, CAD_diagnosis_date),
    by = "id"
  ) %>%
  
  filter(
    event_date >= CAD_diagnosis_date - days(60),
    event_date <= CAD_diagnosis_date + days(30)
  ) %>%
  
  distinct(id) %>%
  mutate(mace_baseline = 1)

ukbdata <- ukbdata %>%
  left_join(baseline_mace_ids, by = "id") %>%
  mutate(
    mace_baseline = ifelse(is.na(mace_baseline), 0, mace_baseline)
  )

# Define follow-up MACE events and time-to-event
library(dplyr)
library(lubridate)

create_drug_specific_mace <- function(
    ukbdata,
    ukbdata_outcome,
    drug_vars,
    drug_lag_days,
    followup_days = 365
){
  
  mace_types <- c("HF", "IS", "MI", "CVD_Death")
  
  ukbdata <- ukbdata %>%
    mutate(
      CAD_diagnosis_date = as.Date(CAD_diagnosis_date)
    )
  
  ukbdata_outcome <- ukbdata_outcome %>%
    mutate(
      event_date = as.Date(event_date)
    )
  
  result <- ukbdata
  
  for (drug in drug_vars) {
    
    cat("Processing:", drug, "\n")
    
    start_col <- paste0(drug, "_start_date")
    lag_days  <- drug_lag_days[[drug]]
    
    if (!(drug %in% names(ukbdata))) next
    if (!(start_col %in% names(ukbdata))) next
    
    # =========================
    # Step 1: Construct the index date
    # =========================
    tmp <- ukbdata %>%
      select(id, all_of(drug), all_of(start_col), CAD_diagnosis_date) %>%
      mutate(
        drug_start_date = as.Date(.data[[start_col]]),
        exposure = .data[[drug]],
        
        index_date = case_when(
          exposure == 1 ~ drug_start_date + days(lag_days),
          exposure == 0 ~ CAD_diagnosis_date + days(lag_days),
          TRUE ~ as.Date(NA)
        ),
        
        window_end = case_when(
          exposure == 1 ~ drug_start_date + days(followup_days),
          exposure == 0 ~ CAD_diagnosis_date + days(followup_days),
          TRUE ~ as.Date(NA)
        )
      ) %>%
      filter(exposure %in% c(0,1) & !is.na(index_date))
    
    if (nrow(tmp) == 0) next
    
    # =========================
    # Step 2: Exclude individuals with baseline MACE
    # =========================
    prior_ids <- ukbdata_outcome %>%
      filter(event_type %in% mace_types) %>%
      inner_join(
        ukbdata %>% select(id, CAD_diagnosis_date),
        by = "id"
      ) %>%
      filter(
        event_date >= CAD_diagnosis_date - days(365),
        event_date <  CAD_diagnosis_date - days(60)
      ) %>%
      distinct(id)
    
    tmp <- tmp %>%
      filter(!id %in% prior_ids$id)
    
    # =========================
    # Step 3: Identify events occurring during follow-up
    # =========================
    tmp_event <- ukbdata_outcome %>%
      filter(event_type %in% mace_types) %>%
      inner_join(tmp %>% select(id, index_date, window_end), by = "id") %>%
      filter(
        event_date >= index_date,
        event_date <= window_end
      )
    
    # =========================
    # Step 4: Identify the earliest event
    # =========================
    mace_df <- tmp_event %>%
      group_by(id) %>%
      summarise(
        mace_date = min(event_date),
        .groups = "drop"
      )
    
    # =========================
    # Step 5: Define the outcome
    # =========================
    one_drug <- tmp %>%
      left_join(mace_df, by = "id") %>%
      mutate(
        MACE = ifelse(!is.na(mace_date), 1, 0),
        
        MACE_TIME = case_when(
          !is.na(mace_date) ~ as.numeric(mace_date - index_date),
          TRUE ~ as.numeric(window_end - index_date)
        )
      ) %>%
      transmute(
        id,
        !!paste0(drug, "_MACE") := MACE,
        !!paste0(drug, "_MACE_TIME") := MACE_TIME
      )
    
    # =========================
    # Step 6：merge
    # =========================
    result <- result %>%
      left_join(one_drug, by = "id")
  }
  
  return(result)
}

ukbdata <- create_drug_specific_mace(
  ukbdata = ukbdata,
  ukbdata_outcome = ukbdata_outcome,
  drug_vars = drug_vars,
  drug_lag_days = drug_lag_days
)



library(dplyr)
library(survival)
library(broom)
library(ggplot2)
library(patchwork)
library(MatchIt)

drug_labels <- c(
  CCB           = "CCB",
  ACEIARBs      = "ACEI/ARBs",
  beta_blocker  = "Beta-blockers",
  metformin     = "Metformin",
  SGLT2         = "SGLT2 inhibitors",
  statin         = "Statin"
)

cluster_labels <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)

drug_colors <- c(
  
  # CCB
  "CCB | Use"              = "#BAB0AC",
  "CCB | Non-use"          = "#D2C9C4",
  
  # ACEI/ARB
  "ACEIARBs | Use"         = "#6D597A",
  "ACEIARBs | Non-use"     = "#B7AABD",
  
  # beta-blocker
  "beta_blocker | Use"     = "#9C755F",
  "beta_blocker | Non-use" = "#C7AA99",
  
  # metformin
  "metformin | Use"        = "#F28E2B",
  "metformin | Non-use"    = "#F4B266",
  
  # SGLT2
  "SGLT2 | Use"            = "#76B7B2",
  "SGLT2 | Non-use"        = "#A5D1CD"
)

line_types <- c(
  "Use" = "solid",
  "Non-use" = "22"
)

drug_superclass_map <- list(
  
  glucose_lowering = c(
    "metformin",
    "sulfonylureas",
    "SGLT2"
  ),
  
  antihypertensive = c(
    "ACEIARBs",
    "beta_blocker",
    "CCB"
  ),
  
  lipid_lowering = c(
    "statin"
  )
)

get_psm_data_cluster <- function(
    data,
    drug,
    cluster_id
){
  
  event_col <- paste0(drug, "_MACE")
  time_col  <- paste0(drug, "_MACE_TIME")
  
  prior_var <- paste0("prior_HF_", drug)
  
  concomitant_var <- paste0(drug, "_concomitant")
  
  superclass_name <- NULL
  
  for(cls in names(drug_superclass_map)){
    
    if(drug %in% drug_superclass_map[[cls]]){
      superclass_name <- cls
    }
  }
  
  if(is.null(superclass_name))
    return(NULL)
  
  #======================================================
  # covariates
  #======================================================
  
  if(superclass_name == "antihypertensive"){
    
    ps_covars <- c(
      "mace_baseline",
      prior_var,
      "age_at_cad",
      "sex",
      "SBP",
      "BMI",
      "eGFR",
      "FBG",
      "TG",
      concomitant_var,
      "prior_CCB",
      "prior_beta_blocker",
      "prior_ACEIARBs"
    )
  }
  
  if(superclass_name == "glucose_lowering"){
    
    ps_covars <- c(
      "mace_baseline",
      prior_var,
      "age_at_cad",
      "sex",
      "FBG",
      "BMI",
      "eGFR",
      "SBP",
      "TG",
      concomitant_var,
      "prior_metformin",
      "prior_sulfonylureas",
      "prior_SGLT2"
    )
  }
  
  if(superclass_name == "lipid_lowering"){
    
    ps_covars <- c(
      "mace_baseline",
      prior_var,
      "age_at_cad",
      "sex",
      "FBG",
      "BMI",
      "eGFR",
      "SBP",
      "TG",
      concomitant_var,
      "prior_statin"
    )
  }
  
  ps_covars <- unique(ps_covars)
  
  ps_covars <- ps_covars[
    ps_covars %in% names(data)
  ]
  
  needed_vars <- unique(c(
    drug,
    event_col,
    time_col,
    ps_covars,
    "clusterpamage"
  ))
  
  dat_sub <- data %>%
    
    filter(clusterpamage == cluster_id) %>%
    
    select(all_of(needed_vars)) %>%
    
    filter(.data[[drug]] %in% c(0,1)) %>%
    
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]])
    ) %>%
    
    drop_na(all_of(ps_covars))
  
  if(nrow(dat_sub) < 30)
    return(NULL)
  
  n_treat <- sum(dat_sub[[drug]] == 1)
  n_ctrl  <- sum(dat_sub[[drug]] == 0)
  
  if(n_treat < 5 | n_ctrl < 5)
    return(NULL)
  
  formula_ps <- as.formula(
    paste(
      drug,
      "~",
      paste(ps_covars, collapse = " + ")
    )
  )
  
  m.out <- MatchIt::matchit(
    formula = formula_ps,
    data = dat_sub,
    method = "nearest",
    ratio = 2,
    caliper = 0.3
  )
  
  matched_data <- MatchIt::match.data(m.out)
  
  return(
    list(
      matched_data = matched_data,
      match_object = m.out,
      ps_covars = ps_covars,
      drug = drug,
      cluster_id = cluster_id
    )
  )
}

# covariate balance assessment
library(cobalt)
library(dplyr)
library(ggplot2)

variable_labels <- c(
  
  age_at_cad = "Age",
  sex_2 = "Male sex",
  
  BMI = "Body mass index",
  SBP = "Systolic blood pressure",
  DBP = "Diastolic blood pressure",
  FBG = "Fasting blood glucose",
  TG = "Triglycerides",
  eGFR = "Estimated glomerular filtration rate",
  
  mace_baseline = "Baseline MACE",
  
  prior_HF_SGLT2 = "Previous heart failure",
  prior_HF_metformin = "Previous heart failure",
  prior_HF_ACEIARBs = "Previous heart failure",
  prior_HF_beta_blocker = "Previous heart failure",
  prior_HF_CCB = "Previous heart failure",
  
  prior_metformin = "Previous metformin use",
  prior_sulfonylureas = "Previous sulfonylurea use",
  prior_SGLT2 = "Previous SGLT2 inhibitor use",
  
  prior_CCB = "Previous CCB use",
  prior_beta_blocker = "Previous β-blocker use",
  prior_ACEIARBs = "Previous ACEI/ARB use",
  
  SGLT2_concomitant = "Concomitant glucose-lowering therapy",
  metformin_concomitant = "Concomitant glucose-lowering therapy",
  
  ACEIARBs_concomitant = "Concomitant antihypertensive therapy",
  beta_blocker_concomitant = "Concomitant antihypertensive therapy",
  CCB_concomitant = "Concomitant antihypertensive therapy",
  
  distance = "Propensity score"
)

plot_psm_balance <- function(
    psm_object,
    drug_label = NULL,
    smd_threshold = 0.10
){
  
  library(cobalt)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  
  m.out <- psm_object$match_object
  
  bal <- bal.tab(
    m.out,
    un = TRUE,
    binary = "std"
  )
  
  bal_df <- bal$Balance %>%
    as.data.frame() %>%
    rownames_to_column("Variable") %>%
    transmute(
      Variable,
      Before = abs(Diff.Un),
      After  = abs(Diff.Adj)
    )
  
  bal_df <- bal_df %>%
    filter(Variable != "distance")
  
  bal_df$Variable <- ifelse(
    bal_df$Variable %in% names(variable_labels),
    variable_labels[bal_df$Variable],
    bal_df$Variable
  )
  
  bal_df <- bal_df %>%
    arrange(desc(Before))
  
  bal_df$Variable <- factor(
    bal_df$Variable,
    levels = rev(bal_df$Variable)
  )
  
  bal_df <- pivot_longer(
    bal_df,
    cols = c(Before, After),
    names_to = "Sample",
    values_to = "SMD"
  )
  
  p <- ggplot(
    bal_df,
    aes(
      SMD,
      Variable,
      color = Sample,
      shape = Sample
    )
  ) +
    
    geom_vline(
      xintercept = smd_threshold,
      linetype = "dashed",
      colour = "grey40",
      linewidth = 0.6
    ) +
    
    geom_point(
      size = 2.8
    ) +
    
    scale_color_manual(
      values = c(
        Before = "#B0B0B0",
        After = "#0072B2"
      )
    ) +
    
    scale_shape_manual(
      values = c(
        Before = 16,
        After = 17
      )
    ) +
    
    scale_x_continuous(
      limits = c(0,0.5),
      breaks = seq(0,0.5,0.1)
    ) +
    
    labs(
      title = drug_label,
      x = "Absolute standardized mean difference",
      y = NULL,
      color = NULL,
      shape = NULL
    ) +
    
    theme_classic(base_size = 12) +
    
    theme(
      
      plot.title = element_text(
        face="bold",
        hjust=0.5,
        size=13
      ),
      
      axis.title.x = element_text(
        face="bold"
      ),
      
      axis.text.y = element_text(
        size=10
      ),
      
      legend.position = "top",
      
      legend.text = element_text(size=10)
    )
  
  return(p)
  
}


#========================================================
# Generate one 4-panel balance plot for each medication
#========================================================

library(patchwork)

cluster_labels <- c(
  "1" = "MH subtype",
  "2" = "HT subtype",
  "3" = "HG subtype",
  "4" = "RI subtype"
)

drug_list <- c(
  "metformin",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

drug_titles <- c(
  metformin = "Metformin",
  SGLT2 = "SGLT2 inhibitors",
  ACEIARBs = "ACEI/ARBs",
  beta_blocker = "β-blockers",
  CCB = "CCB",
  statin = "statin"
)

balance_figures <- list()

for(drug in drug_list){
  
  plot_list <- list()
  
  for(cl in c(1, 2, 3, 4)){
    
    psm <- get_psm_data_cluster(
      data = ukbdata,
      drug = drug,
      cluster_id = cl
    )
    
    if(is.null(psm)){
      next
    }
    
    p <- plot_psm_balance(
      psm_object = psm,
      drug_label = cluster_labels[as.character(cl)]
    ) +
      theme(
        plot.title = element_text(
          face = "bold",
          hjust = 0.5,
          size = 12
        ),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_text(size = 9),
        legend.position = "top"
      )
    
    plot_list[[as.character(cl)]] <- p
  }
  
  final_p <- wrap_plots(
    plot_list,
    ncol = 2,
    guides = "collect"
  ) +
    plot_annotation(
      title = paste0(
        "Covariate balance after propensity score matching: ",
        drug_titles[[drug]]
      )
    ) &
    theme(
      legend.position = "top",
      plot.title = element_text(
        face = "bold",
        size = 15,
        hjust = 0.5
      )
    )
  
  balance_figures[[drug]] <- final_p
  
  print(final_p)
  
  ggsave(
    filename = paste0("PSM_balance_", drug, "_by_subtype.pdf"),
    plot = final_p,
    width = 10,
    height = 8.5,
    device = cairo_pdf
  )
}


# km curve
#========================================================
# Plot cumulative incidence curves
# Panel = Drug
# Line = CAD subtype (Use group only)
#========================================================

library(dplyr)
library(tidyr)
library(survival)
library(broom)
library(ggplot2)
library(scales)

#========================================================
# 0. Settings
#========================================================

cluster_order <- c("MH", "HT", "HG", "RI")

drug_vars_plot <- c(
  "CCB",
  "ACEIARBs",
  "beta_blocker",
  "metformin",
  "SGLT2",
  "statin"
)

drug_titles <- c(
  CCB = "CCB",
  ACEIARBs = "ACEI/ARBs",
  beta_blocker = "\u03B2-blocker",
  metformin = "Metformin",
  SGLT2 = "SGLT2 inhibitor",
  statin = "Statin"
)

cluster_labels <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)

my_colors <- c(
  MH = "#4E79A7",
  HT = "#B07AA1",
  HG = "#E15759",
  RI = "#F28E2B"
)

# max_followup_days <- 300
# days_per_month <- 30.4375
# max_followup_months <- max_followup_days / days_per_month

# Drug-specific effect-lag periods
drug_lag_days <- c(
  CCB = 60,
  ACEIARBs = 60,
  beta_blocker = 60,
  metformin = 90,
  SGLT2 = 90,
  statin = 90
)

max_followup_days <- 300
days_per_month <- 30.4375

# Total plotted time includes the lag period plus post-lag follow-up
max_plot_days <- max_followup_days + max(drug_lag_days)
max_plot_months <- max_plot_days / days_per_month

#========================================================
# 1. Extract cumulative incidence data for one drug
#========================================================

get_single_drug_surv_data <- function(
    data,
    drug_name,
    max_followup_days = 300
) {
  
  all_surv <- list()
  
  lag_days <- unname(drug_lag_days[drug_name])
  
  if (length(lag_days) == 0 || is.na(lag_days)) {
    stop("No lag period specified for drug: ", drug_name)
  }
  
  lag_months <- lag_days / days_per_month
  
  for (cl in 1:4) {
    
    psm <- get_psm_data_cluster(
      data = data,
      drug = drug_name,
      cluster_id = cl
    )
    
    if (is.null(psm)) {
      next
    }
    
    df <- psm$matched_data
    
    outcome_var <- paste0(drug_name, "_MACE")
    time_var <- paste0(drug_name, "_MACE_TIME")
    
    required_vars <- c(drug_name, outcome_var, time_var)
    
    if (!all(required_vars %in% names(df))) {
      warning(
        paste0(
          "Variables missing for ", drug_name,
          ", cluster ", cl, ": ",
          paste(setdiff(required_vars, names(df)), collapse = ", ")
        )
      )
      next
    }
    
    # Only medication users
    df <- df %>%
      filter(
        .data[[drug_name]] == 1,
        !is.na(.data[[outcome_var]]),
        !is.na(.data[[time_var]])
      ) %>%
      mutate(
        subtype = factor(
          cluster_labels[as.character(cl)],
          levels = cluster_order
        ),
        
        # time_raw is assumed to be measured from the end of the lag period
        time_raw = as.numeric(.data[[time_var]]),
        
        # Administrative censoring during the post-lag follow-up period
        postlag_time_days = pmin(
          time_raw,
          max_followup_days
        ),
        
        # Shift the current follow-up time by the drug-specific lag period,
        # so that month 0 corresponds to treatment initiation
        plot_time_days = lag_days + postlag_time_days,
        
        month = plot_time_days / days_per_month,
        
        event = if_else(
          time_raw <= max_followup_days &
            .data[[outcome_var]] == 1,
          1L,
          0L
        )
        
      ) %>%
      filter(
        is.finite(month),
        month >= 0
      )
    
    if (nrow(df) < 5) {
      warning(
        paste0(
          "Fewer than 5 eligible medication users for ",
          drug_name, ", cluster ", cl
        )
      )
      next
    }
    
    fit <- survfit(
      Surv(month, event) ~ 1,
      data = df
    )
    
    surv_df <- broom::tidy(fit) %>%
      transmute(
        month = time,
        cuminc = 1 - estimate,
        subtype = factor(
          cluster_labels[as.character(cl)],
          levels = cluster_order
        ),
        drug = drug_name,
        drug_label = drug_titles[[drug_name]]
      )
    
    # Add two zero-incidence points so that the curve remains flat
    # from treatment initiation to the end of the lag period
    start_df <- tibble(
      month = c(0, lag_months),
      cuminc = c(0, 0),
      subtype = factor(
        cluster_labels[as.character(cl)],
        levels = cluster_order
      ),
      drug = drug_name,
      drug_label = drug_titles[[drug_name]]
    )
    
    all_surv[[as.character(cl)]] <- bind_rows(
      start_df,
      surv_df
    )
  }
  
  bind_rows(all_surv)
}

#========================================================
# 2. Extract data for all drugs
#========================================================

surv_plot_data <- bind_rows(
  lapply(
    drug_vars_plot,
    function(drug_name) {
      get_single_drug_surv_data(
        data = ukbdata,
        drug_name = drug_name,
        max_followup_days = max_followup_days
      )
    }
  )
) %>%
  mutate(
    subtype = factor(
      subtype,
      levels = cluster_order
    ),
    drug_label = factor(
      drug_label,
      levels = unname(drug_titles[drug_vars_plot])
    )
  )

# Check whether any drug panel is missing
missing_drugs <- setdiff(
  unname(drug_titles[drug_vars_plot]),
  unique(as.character(surv_plot_data$drug_label))
)

if (length(missing_drugs) > 0) {
  warning(
    paste(
      "No eligible survival data for:",
      paste(missing_drugs, collapse = ", ")
    )
  )
}

#========================================================
# Export plotting data for panel A
#========================================================

# panelA_plot_data <- surv_plot_data %>%
#   arrange(
#     drug_label,
#     subtype,
#     month
#   )
# 
# write.csv(
#   panelA_plot_data,
#   file = "PanelA_cumulative_incidence_plot_data.csv",
#   row.names = FALSE,
#   fileEncoding = "UTF-8"
# )

#========================================================
# 3. Determine a common y-axis limit
#========================================================

max_cuminc <- max(
  surv_plot_data$cuminc,
  na.rm = TRUE
)

# Add a small amount of space above the largest curve
y_upper <- min(
  1,
  ceiling(max_cuminc * 20) / 20 + 0.02
)

# Avoid overly compressed panels when cumulative incidence is low
y_upper <- max(y_upper, 0.10)

#========================================================
# 4. Nature-style cumulative incidence curves
#    Each panel retains its own x- and y-axes
#    Only left and bottom axis lines are shown
#========================================================

library(ggplot2)
library(scales)
library(grid)
library(ggh4x)

#--------------------------------------------------------
# 4.1 Define common axis ranges
#--------------------------------------------------------

# Common x-axis upper limit
x_upper <- ceiling(max_plot_months)

# Use regular y-axis intervals
y_break_step <- ifelse(
  y_upper <= 0.20,
  0.05,
  0.10
)

y_axis_upper <- ceiling(
  y_upper / y_break_step
) * y_break_step

y_breaks <- seq(
  from = 0,
  to   = y_axis_upper,
  by   = y_break_step
)

#--------------------------------------------------------
# 4.2 Plot
#--------------------------------------------------------

p_cuminc <- ggplot(
  data = surv_plot_data,
  mapping = aes(
    x = month,
    y = cuminc,
    color = subtype,
    group = subtype
  )
) +
  
  # Cumulative incidence curves
  geom_step(
    linewidth = 1.15,
    direction = "hv",
    lineend = "butt",
    na.rm = TRUE
  ) +
  
  # Repeat x- and y-axes in every panel
  ggh4x::facet_wrap2(
    ~ drug_label,
    ncol = 3,
    drop = FALSE,
    axes = "all",
    remove_labels = "none"
  ) +
  
  # Subtype colors
  scale_color_manual(
    values = my_colors,
    limits = cluster_order,
    breaks = cluster_order,
    labels = cluster_order,
    drop = FALSE
  ) +
  
  # X axis
  scale_x_continuous(
    limits = c(0, x_upper),
    breaks = seq(
      from = 0,
      to   = x_upper,
      by   = 2
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  # Y axis
  scale_y_continuous(
    limits = c(0, y_axis_upper),
    breaks = y_breaks,
    labels = percent_format(
      accuracy = 1
    ),
    expand = expansion(
      mult = c(0, 0)
    )
  ) +
  
  labs(
    x = "Time since treatment initiation (months)",
    y = "Cumulative incidence of MACE",
    color = NULL
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "Arial"
  ) +
  
  theme(
    #----------------------------------------------------
    # Facet titles
    #----------------------------------------------------
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 13,
      colour = "black",
      margin = margin(
        t = 2,
        r = 2,
        b = 7,
        l = 2
      )
    ),
    
    #----------------------------------------------------
    # Panel appearance
    #----------------------------------------------------
    panel.border = element_blank(),
    panel.grid = element_blank(),
    
    panel.spacing.x = unit(
      1.2,
      "lines"
    ),
    
    panel.spacing.y = unit(
      1.4,
      "lines"
    ),
    
    #----------------------------------------------------
    # Retain only left and bottom axes in each panel
    #----------------------------------------------------
    axis.line.x = element_line(
      colour = "black",
      linewidth = 0.85,
      lineend = "square"
    ),
    
    axis.line.y = element_line(
      colour = "black",
      linewidth = 0.85,
      lineend = "square"
    ),
    
    #----------------------------------------------------
    # Axis ticks
    #----------------------------------------------------
    axis.ticks.x = element_line(
      colour = "black",
      linewidth = 0.70
    ),
    
    axis.ticks.y = element_line(
      colour = "black",
      linewidth = 0.70
    ),
    
    axis.ticks.length = unit(
      2.3,
      "mm"
    ),
    
    #----------------------------------------------------
    # Tick labels shown in every panel
    #----------------------------------------------------
    axis.text.x = element_text(
      size = 10.5,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 5
      )
    ),
    
    axis.text.y = element_text(
      size = 10.5,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 5
      )
    ),
    
    #----------------------------------------------------
    # Overall shared axis titles
    #----------------------------------------------------
    axis.title.x = element_text(
      size = 13,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 11
      )
    ),
    
    axis.title.y = element_text(
      size = 13,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 11
      )
    ),
    
    #----------------------------------------------------
    # Legend
    #----------------------------------------------------
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 8,
      l = 0
    ),
    
    legend.key = element_blank(),
    
    legend.key.width = unit(
      1.4,
      "cm"
    ),
    
    legend.key.height = unit(
      0.42,
      "cm"
    ),
    
    legend.text = element_text(
      size = 12,
      face = "bold",
      colour = "black"
    ),
    
    #----------------------------------------------------
    # Overall layout
    #----------------------------------------------------
    plot.margin = margin(
      t = 6,
      r = 9,
      b = 8,
      l = 8
    ),
    
    # Slightly square-shaped panels
    aspect.ratio = 0.88
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.65
      )
    )
  )

# Display plot
print(p_cuminc)

ggsave(
  filename = "Cumulative_incidence_by_drug_and_CAD_subtype.pdf",
  plot = p_cuminc,
  width = 10,
  height = 6.8,
  units = "in",
  device = cairo_pdf
)

#========================================================
# PSM + Cox regression
#========================================================
library(survival)
library(dplyr)
library(MatchIt)
library(broom)

drug_vars <- c(
  "metformin",
  "SGLT2",
  "ACEIARBs",
  "beta_blocker",
  "CCB",
  "statin"
)

drug_superclass_map <- list(
  
  glucose_lowering = c(
    "metformin",
    "sulfonylureas",
    "SGLT2"
  ),
  
  antihypertensive = c(
    "ACEIARBs",
    "beta_blocker",
    "CCB"
  ),
  
  lipid_lowering = c(
    "statin"
  )
)

run_drug_cox_models_psm <- function(
    data,
    drug_vars
){
  
  res_list <- list()
  
  for (drug in drug_vars) {
    
    cat("\n=================================\n")
    cat("PSM + Cox:", drug, "\n")
    
    #====================================================
    # outcome
    #====================================================
    
    event_col <- paste0(
      drug,
      "_MACE"
    )
    
    time_col <- paste0(
      drug,
      "_MACE_TIME"
    )
    
    #====================================================
    # exposure-related variables
    #====================================================
    
    prior_var <- paste0(
      "prior_HF_",
      drug
    )
    
    concomitant_var <- paste0(
      drug,
      "_concomitant"
    )
    
    #====================================================
    # superclass
    #====================================================
    
    superclass_name <- NULL
    
    for (cls in names(drug_superclass_map)) {
      
      if (drug %in% drug_superclass_map[[cls]]) {
        
        superclass_name <- cls
      }
    }
    
    if (is.null(superclass_name))
      next
    
    #====================================================
    # outcome-specific PS covariates
    #====================================================
    
    #----------------------------
    # Antihypertensive medications
    #----------------------------
    
    if(superclass_name == "antihypertensive"){
      
      ps_covars <- c(
        
        "mace_baseline",
        
        prior_var,
        
        "age_at_cad",
        "sex",
        
        "SBP",
        "BMI",
        "eGFR",
        "FBG",
        "TG",
        
        concomitant_var,
        
        "prior_CCB",
        "prior_beta_blocker",
        "prior_ACEIARBs"
      )
    }
    
    #----------------------------
    # Glucose-lowering medications
    #----------------------------
    
    if(superclass_name == "glucose_lowering"){
      
      ps_covars <- c(
        
        "mace_baseline",
        
        prior_var,
        
        "age_at_cad",
        "sex",
        
        "FBG",
        "BMI",
        "eGFR",
        "SBP",
        "TG",
        
        concomitant_var,
        
        "prior_metformin",
        "prior_sulfonylureas",
        "prior_SGLT2"
      )
    }
    
    if(superclass_name == "lipid_lowering"){
      
      ps_covars <- c(
        
        "mace_baseline",
        
        prior_var,
        
        "age_at_cad",
        "sex",
        
        "FBG",
        "BMI",
        "eGFR",
        "SBP",
        "TG",
        
        concomitant_var,
        
        "prior_statin"
      )
    }
    
    ps_covars <- unique(ps_covars)
    
    ps_covars <- ps_covars[
      ps_covars %in% names(data)
    ]
    
    needed_vars <- unique(c(
      drug,
      event_col,
      time_col,
      ps_covars
    ))
    
    if (!all(needed_vars %in% names(data)))
      next
    
    dat_sub <- data %>%
      
      select(all_of(needed_vars)) %>%
      
      filter(
        .data[[drug]] %in% c(0,1)
      ) %>%
      
      filter(
        !is.na(.data[[time_col]]),
        !is.na(.data[[event_col]])
      ) %>%
      
      drop_na(all_of(ps_covars))
    
    if (nrow(dat_sub) < 50)
      next
    
    n_treat <- sum(
      dat_sub[[drug]] == 1
    )
    
    n_ctrl <- sum(
      dat_sub[[drug]] == 0
    )
    
    if (n_treat < 10 || n_ctrl < 10)
      next
    
    # PSM
    formula_ps <- as.formula(
      paste(
        drug,
        "~",
        paste(ps_covars, collapse = " + ")
      )
    )
    
    cat("\nPS model:\n")
    print(formula_ps)
    
    m.out <- MatchIt::matchit(
      formula = formula_ps,
      data = dat_sub,
      method = "nearest",
      ratio = 2,
      caliper = 0.3
    )
    
    matched_data <- MatchIt::match.data(
      m.out
    )
    
    n_treat_matched <- sum(
      
      matched_data[[drug]] == 1,
      
      na.rm = TRUE
      
    )
    
    n_ctrl_matched <- sum(
      
      matched_data[[drug]] == 0,
      
      na.rm = TRUE
      
    )
    
    if (nrow(matched_data) < 20)
      next
    
    n_event_treat <- sum(
      matched_data[[event_col]] == 1 &
        matched_data[[drug]] == 1,
      na.rm = TRUE
    )
    
    n_event_ctrl <- sum(
      matched_data[[event_col]] == 1 &
        matched_data[[drug]] == 0,
      na.rm = TRUE
    )
    
    formula_cox <- as.formula(
      paste0(
        "Surv(",
        time_col,
        ", ",
        event_col,
        ") ~ ",
        drug
      )
    )
    
    fit <- coxph(
      formula_cox,
      data = matched_data,
      robust = TRUE
    )
    
    res <- broom::tidy(
      fit,
      conf.int = TRUE
    ) %>%
      filter(term == drug) %>%
      mutate(
        drug = drug,
        estimate = exp(estimate),
        conf.low = exp(conf.low),
        conf.high = exp(conf.high),
        n_total = nrow(dat_sub),
        n_matched = nrow(matched_data),
        n_treat = n_treat,
        n_ctrl = n_ctrl,
        n_event_total = sum(
          matched_data[[event_col]] == 1
        ),
        n_event_treat = n_event_treat,
        n_event_ctrl = n_event_ctrl,
        n_treat_matched = n_treat_matched,
        n_ctrl_matched  = n_ctrl_matched
      ) %>%
      
      select(
        drug,
        estimate,
        conf.low,
        conf.high,
        p.value,
        n_total,
        n_matched,
        n_treat,
        n_ctrl,
        n_treat_matched,
        n_ctrl_matched,
        n_event_total,
        n_event_treat,
        n_event_ctrl
      )
    
    res_list[[drug]] <- res
  }
  
  bind_rows(res_list)
}

#========================================================
# overall
#========================================================

fg_results_psm <- run_drug_cox_models_psm(
  data = ukbdata,
  drug_vars = drug_vars
)

fg_results_psm

#========================================================
# cluster stratified
#========================================================

run_drug_cox_models_psm_by_cluster <- function(
    data,
    drug_vars
){
  
  res_list <- list()
  
  clusters <- sort(
    unique(data$clusterpamage)
  )
  
  for (cl in clusters) {
    
    cat(
      "\n===== Cluster:",
      cl,
      "=====\n"
    )
    
    dat_cl <- data %>%
      filter(clusterpamage == cl)
    
    res_cl <- run_drug_cox_models_psm(
      data = dat_cl,
      drug_vars = drug_vars
    ) %>%
      
      mutate(
        cluster = paste0(
          "Cluster ",
          cl
        )
      )
    
    res_list[[as.character(cl)]] <- res_cl
  }
  
  bind_rows(res_list)
}


fg_results_cluster <- run_drug_cox_models_psm_by_cluster(
  data = ukbdata,
  drug_vars = drug_vars
)

fg_results_cluster


write.csv(
  fg_results_cluster,
  "Figure 4-B-drug_MACE_results-main.csv",
  row.names = FALSE
)


# RHR calculate
library(dplyr)

rhr_results <- fg_results_cluster %>%
  
  mutate(
    cluster = gsub("Cluster ","",cluster),
    cluster = factor(cluster,
                     levels=c("1","2","3","4"),
                     labels=c("MH","HT","HG","RI")),
    
    logHR = log(estimate),
    
    se = (log(conf.high)-log(conf.low))/(2*1.96)
  ) %>%
  
  group_by(drug) %>%
  
  group_modify(~{
    
    mh <- filter(.x,cluster=="MH")
    
    if(nrow(mh)==0) return(NULL)
    
    mh_loghr <- mh$logHR
    mh_se    <- mh$se
    
    .x %>%
      
      mutate(
        
        logRHR = logHR-mh_loghr,
        
        se_rhr = sqrt(se^2+mh_se^2),
        
        RHR = exp(logRHR),
        
        lower = exp(logRHR-1.96*se_rhr),
        
        upper = exp(logRHR+1.96*se_rhr),
        
        z = logRHR/se_rhr,
        
        p = 2*pnorm(-abs(z))
        
      )
    
  }) %>%
  
  ungroup()

rhr_results

#========================================================
# Pairwise RHR heatmap
# Based on fg_results_cluster
#
# RHR = HR(comparison subtype) / HR(reference subtype)
#
# Full 4 × 4 matrix:
# - diagonal cells are blank and white
# - missing estimates are displayed as -Inf
# - y-axis order from top to bottom: MH, HT, HG, RI
#========================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(openxlsx)
library(grid)

#========================================================
# 0. Settings
#========================================================

cluster_order <- c("MH", "HT", "HG", "RI")

# ggplot places the first y-axis factor level at the bottom,
# so reverse the levels to obtain MH → HT → HG → RI from top to bottom
comparison_order <- rev(cluster_order)

drug_order <- c(
  "CCB",
  "ACEIARBs",
  "beta_blocker",
  "metformin",
  "SGLT2",
  "statin"
)

drug_labels <- c(
  CCB = "CCB",
  ACEIARBs = "ACEI/ARBs",
  beta_blocker = "\u03B2-blocker",
  metformin = "Metformin",
  SGLT2 = "SGLT2 inhibitor",
  statin = "Statin"
)

# Whether HR estimates with zero events in either group
# should be treated as unavailable
exclude_zero_event_estimates <- TRUE

# Heatmap color limits based on log(RHR)
fill_min <- log(0.25)
fill_max <- log(4)

#========================================================
# 1. Prepare subtype-specific HR results
#========================================================

hr_dat <- fg_results_cluster %>%
  mutate(
    cluster = gsub("Cluster\\s*", "", cluster),
    
    cluster = factor(
      cluster,
      levels = c("1", "2", "3", "4"),
      labels = cluster_order
    ),
    
    drug = factor(
      drug,
      levels = drug_order
    ),
    
    # Flag zero-event estimates
    zero_event_flag =
      is.na(n_event_treat) |
      is.na(n_event_ctrl) |
      n_event_treat == 0 |
      n_event_ctrl == 0,
    
    # Calculate log(HR) and its standard error only when valid
    logHR = case_when(
      is.na(estimate) ~ NA_real_,
      !is.finite(estimate) ~ NA_real_,
      estimate <= 0 ~ NA_real_,
      TRUE ~ log(estimate)
    ),
    
    se = case_when(
      is.na(conf.low) | is.na(conf.high) ~ NA_real_,
      !is.finite(conf.low) | !is.finite(conf.high) ~ NA_real_,
      conf.low <= 0 | conf.high <= 0 ~ NA_real_,
      TRUE ~ (log(conf.high) - log(conf.low)) / (2 * 1.96)
    )
  )

# Optionally remove HR estimates based on zero-event groups
if (exclude_zero_event_estimates) {
  
  hr_dat <- hr_dat %>%
    mutate(
      estimate = if_else(
        zero_event_flag,
        NA_real_,
        estimate
      ),
      
      conf.low = if_else(
        zero_event_flag,
        NA_real_,
        conf.low
      ),
      
      conf.high = if_else(
        zero_event_flag,
        NA_real_,
        conf.high
      ),
      
      logHR = if_else(
        zero_event_flag,
        NA_real_,
        logHR
      ),
      
      se = if_else(
        zero_event_flag,
        NA_real_,
        se
      )
    )
}

#========================================================
# 2. Calculate the complete pairwise RHR matrix
#========================================================

pairwise_rhr <- hr_dat %>%
  select(
    drug,
    cluster,
    estimate,
    conf.low,
    conf.high,
    logHR,
    se,
    n_total,
    n_treat,
    n_ctrl,
    n_event_total,
    n_event_treat,
    n_event_ctrl,
    zero_event_flag
  ) %>%
  group_by(drug) %>%
  group_modify(~{
    
    dat <- .x
    
    # Ensure one record per subtype
    dat <- dat %>%
      mutate(
        cluster_character = as.character(cluster)
      ) %>%
      arrange(
        match(cluster_character, cluster_order)
      )
    
    # Generate all 4 × 4 combinations
    pair_grid <- expand_grid(
      reference = cluster_order,
      comparison = cluster_order
    )
    
    out <- pair_grid %>%
      rowwise() %>%
      mutate(
        ref_logHR = dat$logHR[
          match(
            reference,
            dat$cluster_character
          )
        ],
        
        ref_se = dat$se[
          match(
            reference,
            dat$cluster_character
          )
        ],
        
        comp_logHR = dat$logHR[
          match(
            comparison,
            dat$cluster_character
          )
        ],
        
        comp_se = dat$se[
          match(
            comparison,
            dat$cluster_character
          )
        ],
        
        diagonal = reference == comparison,
        
        # Only calculate off-diagonal comparisons
        logRHR = case_when(
          diagonal ~ NA_real_,
          is.na(ref_logHR) | is.na(comp_logHR) ~ NA_real_,
          TRUE ~ comp_logHR - ref_logHR
        ),
        
        # This assumes subtype-specific log-HR estimates are independent
        se_logRHR = case_when(
          diagonal ~ NA_real_,
          is.na(ref_se) | is.na(comp_se) ~ NA_real_,
          TRUE ~ sqrt(comp_se^2 + ref_se^2)
        ),
        
        RHR = case_when(
          diagonal ~ NA_real_,
          is.na(logRHR) ~ NA_real_,
          TRUE ~ exp(logRHR)
        ),
        
        lower = case_when(
          diagonal ~ NA_real_,
          is.na(logRHR) | is.na(se_logRHR) ~ NA_real_,
          TRUE ~ exp(logRHR - 1.96 * se_logRHR)
        ),
        
        upper = case_when(
          diagonal ~ NA_real_,
          is.na(logRHR) | is.na(se_logRHR) ~ NA_real_,
          TRUE ~ exp(logRHR + 1.96 * se_logRHR)
        ),
        
        z = case_when(
          diagonal ~ NA_real_,
          is.na(logRHR) | is.na(se_logRHR) ~ NA_real_,
          se_logRHR <= 0 ~ NA_real_,
          TRUE ~ logRHR / se_logRHR
        ),
        
        p.value = case_when(
          diagonal ~ NA_real_,
          is.na(z) ~ NA_real_,
          TRUE ~ 2 * pnorm(-abs(z))
        )
      ) %>%
      ungroup() %>%
      mutate(
        reference = factor(
          reference,
          levels = cluster_order
        ),
        
        # Reversed factor levels give top-to-bottom MH, HT, HG, RI
        comparison = factor(
          comparison,
          levels = comparison_order
        )
      )
    
    out
  }) %>%
  ungroup()

#========================================================
# 3. Format labels and fill values
#========================================================

pairwise_rhr_plot <- pairwise_rhr %>%
  mutate(
    drug_label = recode(
      as.character(drug),
      !!!drug_labels
    ),
    
    drug_label = factor(
      drug_label,
      levels = unname(drug_labels[drug_order])
    ),
    
    # Identify unavailable off-diagonal comparisons
    unavailable = !diagonal & (
      is.na(RHR) |
        is.na(lower) |
        is.na(upper)
    ),
    
    # Significance labels
    sig = case_when(
      diagonal ~ "",
      unavailable ~ "",
      is.na(p.value) ~ "",
      p.value < 0.001 ~ "***",
      p.value < 0.01 ~ "**",
      p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    
    # Cell labels
    cell_label = case_when(
      # Diagonal is blank
      diagonal ~ "",
      
      # Missing estimates are shown as -Inf
      unavailable ~ "-Inf",
      
      # Available estimates
      TRUE ~ sprintf(
        "%.2f\n(%.2f\u2013%.2f)%s",
        RHR,
        lower,
        upper,
        sig
      )
    ),
    
    # Fill values:
    # - diagonal: NA, subsequently displayed as white
    # - unavailable estimates: color lower limit
    # - valid estimates: truncated log(RHR)
    fill_value = case_when(
      diagonal ~ NA_real_,
      
      unavailable ~ fill_min,
      
      logRHR < fill_min ~ fill_min,
      
      logRHR > fill_max ~ fill_max,
      
      TRUE ~ logRHR
    )
  )

#========================================================
# Add separate labels and text color
#========================================================

pairwise_rhr_plot <- pairwise_rhr_plot %>%
  mutate(
    # First line: RHR estimate
    hr_label = case_when(
      diagonal ~ "",
      unavailable ~ "NE",
      TRUE ~ sprintf("%.2f", RHR)
    ),
    
    # Second line: 95% CI and significance symbols
    ci_label = case_when(
      diagonal ~ "",
      unavailable ~ "",
      TRUE ~ sprintf(
        "(%.2f\u2013%.2f)%s",
        lower,
        upper,
        sig
      )
    ),
    
    # Text color according to background intensity
    # Unavailable cells retain the dark-blue fill and use white text
    text_colour = case_when(
      diagonal ~ "black",
      unavailable ~ "white",
      fill_value <= log(0.50) ~ "white",
      fill_value >= log(2.00) ~ "white",
      TRUE ~ "black"
    )
  )

#========================================================
# Export plotting data for panel B
#========================================================

panelB_plot_data <- pairwise_rhr_plot %>%
  arrange(
    drug_label,
    reference,
    comparison
  )

write.csv(
  panelB_plot_data,
  file = "PanelB_pairwise_RHR_heatmap_data.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

#========================================================
# 4. pairwise RHR heatmap
#    - RHR estimate shown in larger font
#    - 95% CI removed
#    - Significance symbols retained
#    - NE used for unavailable estimates
#    - Unavailable-cell fill colors remain unchanged
#    - White labels used on dark blue/red cells
#    - Each panel retains x- and y-axis labels
#    - Drug titles positioned above the top x-axis
#========================================================

library(ggplot2)
library(grid)
library(scales)
library(ggh4x)
library(dplyr)

#--------------------------------------------------------
# Create label: RHR estimate + significance symbols
#--------------------------------------------------------

pairwise_rhr_plot <- pairwise_rhr_plot %>%
  mutate(
    hr_sig_label = case_when(
      unavailable ~ "NE",
      diagonal ~ "",
      TRUE ~ paste0(
        hr_label,
        ifelse(is.na(sig), "", sig)
      )
    )
  )

#========================================================
# Heatmap
#========================================================

p_rhr_heatmap <- ggplot(
  data = pairwise_rhr_plot,
  mapping = aes(
    x = reference,
    y = comparison
  )
) +
  
  #------------------------------------------------------
# Heatmap cells
#------------------------------------------------------
geom_tile(
  aes(fill = fill_value),
  colour = "grey70",
  linewidth = 0.4,
  width = 1,
  height = 1
) +
  
  #------------------------------------------------------
# RHR estimate + significance symbols
#------------------------------------------------------
geom_text(
  aes(
    label = hr_sig_label,
    colour = text_colour
  ),
  size = 3.6,
  fontface = "bold",
  family = "Arial",
  show.legend = FALSE
) +
  
  #------------------------------------------------------
# Repeat x- and y-axis labels in every panel
#------------------------------------------------------
ggh4x::facet_wrap2(
  ~ drug_label,
  nrow = 2,
  ncol = 3,
  drop = FALSE,
  axes = "all",
  remove_labels = "none",
  strip.position = "top"
) +
  
  #------------------------------------------------------
# Reference subtype labels on top
#------------------------------------------------------
scale_x_discrete(
  drop = FALSE,
  position = "bottom"
) +
  
  scale_y_discrete(
    drop = FALSE
  ) +
  
  #------------------------------------------------------
# Diverging color scale
#------------------------------------------------------
scale_fill_gradient2(
  low = "#4A90C2",
  mid = "#F7F7F7",
  high = "#D6604D",
  midpoint = 0,
  limits = c(
    fill_min,
    fill_max
  ),
  breaks = log(
    c(
      0.25,
      0.5,
      1,
      2,
      4
    )
  ),
  labels = c(
    "0.25",
    "0.5",
    "1",
    "2",
    "4"
  ),
  name = "RHR",
  na.value = "white",
  oob = scales::squish
) +
  
  # Interpret "black" and "white" directly
  scale_colour_identity() +
  
  #------------------------------------------------------
# Square cells
#------------------------------------------------------
coord_fixed(
  ratio = 1,
  clip = "off"
) +
  
  labs(
    x = "Reference subtype",
    y = "Comparison subtype",
    fill = "RHR"
  ) +
  
  theme_minimal(
    base_size = 11,
    base_family = "Arial"
  ) +
  
  theme(
    #----------------------------------------------------
    # Drug titles
    #----------------------------------------------------
    strip.background = element_blank(),
    strip.placement = "outside",
    
    strip.text.x = element_text(
      face = "bold",
      size = 12.5,
      colour = "black",
      hjust = 0.5,
      margin = margin(
        t = 2,
        r = 2,
        b = 8,
        l = 2
      )
    ),
    
    #----------------------------------------------------
    # Panel appearance
    #----------------------------------------------------
    panel.grid = element_blank(),
    panel.border = element_blank(),
    
    panel.background = element_rect(
      fill = "white",
      colour = NA
    ),
    
    #----------------------------------------------------
    # Top x-axis labels
    #----------------------------------------------------
    axis.text.x.bottom = element_text(
      size = 10.5,
      face = "bold",
      colour = "black",
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(
        t = 4
      )
    ),
    
    #----------------------------------------------------
    # Left y-axis labels
    #----------------------------------------------------
    axis.text.y.left = element_text(
      size = 10.5,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 5
      )
    ),
    
    # Hide duplicate lower/right labels
    axis.text.y.right = element_blank(),
    
    #----------------------------------------------------
    # Axis ticks
    #----------------------------------------------------
    axis.ticks.x.bottom = element_line(
      colour = "black",
      linewidth = 0.45
    ),
  
    
    axis.ticks.y.left = element_line(
      colour = "black",
      linewidth = 0.45
    ),
    
    axis.ticks.x.top = element_blank(),
    axis.ticks.y.right = element_blank(),
    
    axis.ticks.length = unit(
      1.6,
      "mm"
    ),
    
    #----------------------------------------------------
    # Shared axis titles
    #----------------------------------------------------
    axis.title.x = element_text(
      size = 12.5,
      face = "bold",
      colour = "black",
      margin = margin(
        t = 8
      )
    ),
    
    axis.title.y = element_text(
      size = 12.5,
      face = "bold",
      colour = "black",
      margin = margin(
        r = 9
      )
    ),
    
    #----------------------------------------------------
    # Legend
    #----------------------------------------------------
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 11,
      colour = "black"
    ),
    
    legend.text = element_text(
      size = 9.5,
      colour = "black"
    ),
    
    legend.margin = margin(
      t = 0,
      r = 0,
      b = 0,
      l = 8
    ),
    
    #----------------------------------------------------
    # Panel spacing and margins
    #----------------------------------------------------
    panel.spacing.x = unit(
      1.1,
      "lines"
    ),
    
    panel.spacing.y = unit(
      1.5,
      "lines"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 10,
      l = 8
    )
  ) +
  
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(
        4.6,
        "cm"
      ),
      barwidth = unit(
        0.45,
        "cm"
      ),
      ticks = TRUE,
      frame.colour = "black",
      frame.linewidth = 0.35
    )
  )

print(p_rhr_heatmap)



#========================================================
# 5. Export PDF
#========================================================

ggsave(
  filename = "pairwise_RHR_heatmap.pdf",
  plot = p_rhr_heatmap,
  width = 10,
  height = 6.5,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

#========================================================
# 6. Export Excel table
#========================================================

pairwise_rhr_table <- pairwise_rhr %>%
  mutate(
    drug = as.character(drug),
    reference = as.character(reference),
    comparison = as.character(comparison),
    
    `RHR (95% CI)` = if_else(
      is.na(RHR),
      NA_character_,
      sprintf("%.2f (%.2f, %.2f)", RHR, lower, upper)
    ),
    
    `P value` = case_when(
      is.na(p.value) ~ NA_character_,
      p.value < 0.001 ~ "<0.001",
      TRUE ~ sprintf("%.3f", p.value)
    )
  ) %>%
  select(
    drug,
    reference,
    comparison,
    RHR,
    lower,
    upper,
    p.value,
    `RHR (95% CI)`,
    `P value`
  )

write.xlsx(
  list(
    "Pairwise RHR" = pairwise_rhr_table
  ),
  file = "pairwise_RHR_table.xlsx",
  overwrite = TRUE
)

pairwise_rhr_table
