###############################################################################
# CardioMet: Prognostic analysis across study cohorts
#
# Study:
# Unsupervised cardiometabolic phenotyping reveals prognostic and 
# treatment-associated heterogeneity in coronary artery disease
#
# Description:
# This script provides the common prognostic-analysis workflow used to
# evaluate cardiovascular outcomes across the four CardioMet cardiometabolic
# subtypes in the BIPASS discovery cohort and the IMPROVE, Mount Sinai Data
# Warehouse (MSDW), and UK Biobank (UKB) validation cohorts.
#
# The primary outcome is major adverse cardiovascular events (MACE).
# Non-cardiovascular death is treated as a competing event when applicable.
# Cohort-specific outcome definitions, follow-up variables, and data
# harmonization are completed before application of this common analytical
# workflow.
#
# The workflow includes:
#   1. Construction of the competing-risk outcome and follow-up time
#   2. Estimation of subtype-specific cumulative incidence of MACE using the
#      Aalen–Johansen estimator
#   3. Calculation of MACE event counts and crude event proportions by subtype
#   4. Overall comparison of cumulative incidence functions using Gray's test
#   5. Pairwise Gray's tests between cardiometabolic subtypes with correction
#      for multiple comparisons
#   6. Fine–Gray competing-risk regression to estimate subtype-specific
#      subdistribution hazard ratios for MACE, adjusted for age and sex
#   7. Estimation of follow-up duration using the reverse Kaplan–Meier method
#
# Cardiometabolic subtypes:
#   - MH: metabolically healthy
#   - HT: hypertensive
#   - HG: hyperglycemic
#   - RI: renal impairment
#
# Outcome:
#   - MACE: major adverse cardiovascular events
#
# Competing event:
#   - Non-cardiovascular death
#
# Main analytical methods:
#   - Aalen–Johansen cumulative incidence estimator
#   - Gray's test for equality of cumulative incidence functions
#   - Pairwise Gray's tests with multiple-comparison adjustment
#   - Fine–Gray subdistribution hazard regression
#   - Reverse Kaplan–Meier estimation of follow-up duration
#
# Input:
#   - A cohort-specific harmonized analytical dataset containing:
#       * participant identifier
#       * CardioMet subtype assignment
#       * MACE indicator
#       * time to MACE or censoring
#       * non-cardiovascular death indicator and event time
#       * age
#       * sex
#
# Output:
#   - Subtype-specific 1-year cumulative incidence estimates for MACE
#   - MACE event counts and crude event proportions by subtype
#   - Overall and pairwise Gray's test results
#   - Fine–Gray subdistribution hazard ratios with 95% confidence intervals
#   - Follow-up duration estimates
#
# Notes:
#   - The same prognostic-analysis framework is applied across the BIPASS,
#     IMPROVE, MSDW, and UKB cohorts after cohort-specific data preparation
#     and outcome harmonization.
#   - Follow-up duration and outcome ascertainment may differ across cohorts
#     according to study design and data availability.
#   - MACE is coded as the event of interest (status = 1),
#     non-cardiovascular death as the competing event (status = 2), and
#     event-free or censored observations as status = 0.
#   - The MH subtype is used as the reference category in Fine–Gray regression.
#   - Pairwise Gray's tests are adjusted using the Bonferroni method.
#     Benjamini–Hochberg adjusted P values are additionally calculated.
#
# Author:
# Xinhui Liu
#
# Institution:
# Department of Emergency Medicine,
# Qilu Hospital, Shandong University,
# Jinan, Shandong, China
#
# Correspondence:
# Jiali Wang, MD, PhD
# Yuguo Chen, MD, PhD
#
# Repository:
# https://github.com/lxinhui/CardioMet
#
###############################################################################


###############################################################################
# 1. Construct competing-risk outcome and estimate 1-year cumulative incidence
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(cmprsk)
library(tibble)

mainout  <- c("COMP")
maintime <- c("COMP_INHOSP")

# Merge cardiovascular outcomes and competing events
data_outcome <- data_mean_imputed[, c(
  "id",
  mainout,
  maintime,
  "NONCVD_DEATH_TIME",
  "NONCVD_DEATH"
)]

data_cluster <- merge(
  data_cluster,
  data_outcome,
  all.x = TRUE,
  by = "id"
)

# Restrict analysis to 1 year
data_cluster$time_1yr <- 365

table(data_cluster$COMP)


# Construct follow-up time for competing-risk analysis
data_cluster$time_COMP_CIF <- pmin(
  ifelse(
    data_cluster$COMP == 1,
    data_cluster$COMP_INHOSP,
    Inf
  ),
  ifelse(
    data_cluster$NONCVD_DEATH == 1,
    data_cluster$NONCVD_DEATH_TIME,
    Inf
  ),
  data_cluster$time_1yr,
  na.rm = TRUE
)


# Status:
# 0 = event-free/censored
# 1 = MACE
# 2 = non-cardiovascular death
data_cluster$status_COMP_CIF <- with(
  data_cluster,
  ifelse(
    time_COMP_CIF == COMP_INHOSP & COMP == 1,
    1,
    ifelse(
      time_COMP_CIF == NONCVD_DEATH_TIME &
        NONCVD_DEATH == 1,
      2,
      0
    )
  )
)


# Competing-event mapping
comp_event_map <- list(
  COMP = list(
    comp = "status_COMP_CIF",
    comptime = "time_COMP_CIF"
  )
)

outcome_labels <- c(
  COMP = "MACE"
)


# Function to calculate cumulative incidence
compute_cif <- function(
    event_var,
    data,
    time_point_years = 1
) {
  
  comp_var <- comp_event_map[[event_var]]$comp
  comptime_var <- comp_event_map[[event_var]]$comptime
  
  df <- data %>%
    mutate(
      ftime = .data[[comptime_var]],
      fstatus = .data[[comp_var]]
    ) %>%
    select(
      id,
      cluster,
      ftime,
      fstatus
    ) %>%
    filter(
      !is.na(ftime),
      !is.na(fstatus),
      !is.na(cluster)
    )
  
  cif_fit <- tryCatch(
    cmprsk::cuminc(
      ftime = df$ftime,
      fstatus = df$fstatus,
      group = df$cluster
    ),
    error = function(e) {
      message(
        "cuminc failed for ",
        event_var,
        ": ",
        e$message
      )
      return(NULL)
    }
  )
  
  if (is.null(cif_fit)) {
    return(NULL)
  }
  
  cif_long <- purrr::map_dfr(
    names(cif_fit),
    function(name) {
      
      if (grepl("Tests", name)) {
        return(NULL)
      }
      
      est <- cif_fit[[name]]
      parts <- strsplit(name, " ")[[1]]
      
      tibble(
        cluster = parts[1],
        event = parts[2],
        time = est$time,
        CIF = est$est * 100
      )
    }
  )
  
  time_point <- time_point_years * 365
  
  cif_1y <- cif_long %>%
    filter(event == "1") %>%
    group_by(cluster) %>%
    slice_min(
      abs(time - time_point),
      n = 1,
      with_ties = FALSE
    ) %>%
    summarise(
      Rate = mean(CIF, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Outcome = event_var
    )
  
  return(cif_1y)
}


# Estimate 1-year cumulative incidence
event_summary_list <- map(
  mainout,
  ~ compute_cif(
    event_var = .x,
    data = data_cluster,
    time_point_years = 1
  )
)

event_summary_long <- bind_rows(
  event_summary_list
)


cluster_map <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)


event_summary_long <- event_summary_long %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = names(outcome_labels)
    ),
    
    Outcome_label =
      outcome_labels[as.character(Outcome)],
    
    cluster = factor(
      cluster_map[as.character(cluster)],
      levels = c(
        "MH",
        "HT",
        "HG",
        "RI"
      )
    )
  )


cat(
  "\n1-Year Cumulative Incidence by Subtype:\n"
)

print(
  event_summary_long,
  n = 20
)


###############################################################################
# 2. Calculate MACE events and sample size by subtype
###############################################################################

case_total_by_cluster <- data_cluster %>%
  filter(
    !is.na(cluster)
  ) %>%
  mutate(
    
    cluster_label = factor(
      cluster_map[as.character(cluster)],
      levels = c(
        "MH",
        "HT",
        "HG",
        "RI"
      )
    ),
    
    mace_1yr = ifelse(
      status_COMP_CIF == 1,
      1,
      0
    ),
    
    competing_death_1yr = ifelse(
      status_COMP_CIF == 2,
      1,
      0
    ),
    
    censored_or_event_free_1yr = ifelse(
      status_COMP_CIF == 0,
      1,
      0
    )
    
  ) %>%
  group_by(
    cluster_label
  ) %>%
  summarise(
    
    total_n = n(),
    
    mace_cases = sum(
      mace_1yr,
      na.rm = TRUE
    ),
    
    competing_deaths = sum(
      competing_death_1yr,
      na.rm = TRUE
    ),
    
    censored_or_event_free = sum(
      censored_or_event_free_1yr,
      na.rm = TRUE
    ),
    
    crude_mace_percent =
      mace_cases / total_n * 100,
    
    case_total = paste0(
      mace_cases,
      "/",
      total_n
    ),
    
    case_total_percent = paste0(
      mace_cases,
      "/",
      total_n,
      " (",
      sprintf(
        "%.1f",
        crude_mace_percent
      ),
      "%)"
    ),
    
    .groups = "drop"
  )


print(case_total_by_cluster)


###############################################################################
# 3. Overall and pairwise Gray's tests
###############################################################################

outcome_vars <- c("COMP")
outcome_labels <- c("MACE")

comp_event_map <- list(
  COMP = list(
    comp = "status_COMP_CIF",
    comptime = "time_COMP_CIF"
  )
)

pval_list <- list()


for (i in seq_along(outcome_vars)) {
  
  message(
    "Processing outcome ",
    i,
    ": ",
    outcome_labels[i]
  )
  
  main_event <- outcome_vars[i]
  
  comp_event <-
    comp_event_map[[main_event]]$comp
  
  comp_time <-
    comp_event_map[[main_event]]$comptime
  
  
  if (
    !all(
      c(
        comp_event,
        comp_time
      ) %in% names(data_cluster)
    )
  ) {
    
    message(
      "Skipped: missing required variables for ",
      outcome_labels[i]
    )
    
    next
  }
  
  
  tmp_df <- data_cluster %>%
    mutate(
      ftime = .data[[comp_time]],
      fstatus = .data[[comp_event]],
      group = factor(cluster)
    ) %>%
    filter(
      !is.na(ftime),
      !is.na(fstatus),
      !is.na(group)
    )
  
  
  if (
    nrow(tmp_df) < 50 ||
    length(unique(tmp_df$group)) < 2
  ) {
    next
  }
  
  
  # Overall Gray's test
  fit_overall <- tryCatch(
    
    cuminc(
      ftime = tmp_df$ftime,
      fstatus = tmp_df$fstatus,
      group = tmp_df$group
    ),
    
    error = function(e) NULL
  )
  
  
  if (
    is.null(fit_overall) ||
    is.null(fit_overall$Tests)
  ) {
    next
  }
  
  
  gray_test <- as.data.frame(
    fit_overall$Tests
  )
  
  overall_p <- gray_test$pv[1]
  
  
  # Pairwise Gray's tests
  clusters <- sort(
    unique(tmp_df$group)
  )
  
  pairwise_results <- list()
  
  
  for (
    a in 1:(length(clusters) - 1)
  ) {
    
    for (
      b in (a + 1):length(clusters)
    ) {
      
      sub_df <- tmp_df %>%
        filter(
          group %in%
            c(
              clusters[a],
              clusters[b]
            )
        )
      
      
      fit_sub <- tryCatch(
        
        cuminc(
          ftime = sub_df$ftime,
          fstatus = sub_df$fstatus,
          group = sub_df$group
        ),
        
        error = function(e) NULL
      )
      
      
      if (
        is.null(fit_sub) ||
        is.null(fit_sub$Tests)
      ) {
        next
      }
      
      
      sub_test <- as.data.frame(
        fit_sub$Tests
      )
      
      pval <- sub_test$pv[1]
      
      
      pairwise_results[
        [
          paste0(
            clusters[a],
            "_vs_",
            clusters[b]
          )
        ]
      ] <- data.frame(
        
        Group1 = clusters[a],
        Group2 = clusters[b],
        p_raw = pval
      )
    }
  }
  
  
  if (
    length(pairwise_results) == 0
  ) {
    next
  }
  
  
  # Multiple-comparison correction
  p_df <- bind_rows(
    pairwise_results
  ) %>%
    mutate(
      
      p_bonferroni = p.adjust(
        p_raw,
        method = "bonferroni"
      ),
      
      p_bh = p.adjust(
        p_raw,
        method = "BH"
      ),
      
      Outcome =
        outcome_labels[i],
      
      p_overall =
        overall_p
    )
  
  
  pval_list[
    [outcome_labels[i]]
  ] <- p_df
}


pval_table <- bind_rows(
  pval_list
)

cat(
  "\nOverall and Pairwise Gray's Test Results:\n"
)

print(pval_table)


###############################################################################
# 4. Fine–Gray competing-risk regression
###############################################################################

library(broom)

mainout <- c("COMP")

comp_event_map <- list(
  COMP = list(
    comp = "status_COMP_CIF",
    comptime = "time_COMP_CIF"
  )
)


# Recode cluster labels
data_cluster$cluster <- factor(
  data_cluster$cluster
)

cluster_map <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)

data_cluster$cluster <- factor(
  cluster_map[
    as.character(
      data_cluster$cluster
    )
  ],
  levels = c(
    "MH",
    "HT",
    "HG",
    "RI"
  )
)

data_cluster$cluster <- relevel(
  data_cluster$cluster,
  ref = "MH"
)


# Add sex variable if required
data_add <- datafix[, c(
  "id",
  "sex"
)]

data_cluster <- merge(
  data_cluster,
  data_add,
  all.x = TRUE,
  by = "id"
)


# Fine–Gray regression function
run_finegray <- function(
    comp_event,
    comp_time
) {
  
  df_sub <- data_cluster %>%
    mutate(
      time =
        .data[[comp_time]],
      
      status =
        .data[[comp_event]]
    )
  
  
  model_matrix <- model.matrix(
    ~ cluster + sex.y + age,
    data = df_sub
  )[, -1, drop = FALSE]
  
  
  model <- crr(
    ftime = df_sub$time,
    fstatus = df_sub$status,
    cov1 = model_matrix
  )
  
  
  res <- broom::tidy(
    model,
    conf.int = TRUE,
    exponentiate = TRUE
  )
  
  
  res %>%
    filter(
      grepl(
        "cluster",
        term
      )
    ) %>%
    mutate(
      
      term = gsub(
        "^cluster",
        "",
        term
      ),
      
      HR_CI = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      )
    ) %>%
    select(
      term,
      estimate,
      conf.low,
      conf.high,
      HR_CI,
      p.value
    )
}


results_list <- list()


for (
  outcome in mainout
) {
  
  comp_event <-
    comp_event_map[[outcome]]$comp
  
  comp_time <-
    comp_event_map[[outcome]]$comptime
  
  
  res_i <- run_finegray(
    comp_event = comp_event,
    comp_time = comp_time
  )
  
  
  res_i$outcome <- outcome
  
  results_list[[outcome]] <-
    res_i
}


final_result <- bind_rows(
  results_list
) %>%
  mutate(
    outcome = case_when(
      outcome == "COMP" ~ "MACE",
      TRUE ~ outcome
    )
  )


cat(
  "\nFine-Gray Competing-Risk Regression Results:\n"
)

print(final_result)


###############################################################################
# 5. Estimate median follow-up duration
###############################################################################

library(survival)

fit_followup <- survfit(
  
  Surv(
    time_COMP_CIF,
    status_COMP_CIF == 0
  ) ~ 1,
  
  data = data_cluster
)


median_followup <- fit_followup$table[
  "median"
]


cat(
  "\nMedian Follow-up Duration:\n"
)

print(
  median_followup
)
