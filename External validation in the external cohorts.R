###############################################################################
# CardioMet: External validation pipeline
#
# Study:
# Unsupervised cardiometabolic phenotyping reveals prognostic and 
# treatment-associated heterogeneity in coronary artery disease
#
# Description:
# This script provides the common analytical pipeline used for external
# validation of the CardioMet phenotyping framework in independent cohorts
# of participants with coronary artery disease. The same analytical workflow
# was applied to the IMPROVE, Mount Sinai Data Warehouse (MSDW), and UK
# Biobank (UKB) validation cohorts after cohort-specific data extraction,
# preprocessing, and harmonization.
#
# The workflow includes:
#   1. Loading a harmonized external-validation dataset
#   2. Exclusion of clinically implausible values
#   3. Age- and sex-adjustment and standardization of clustering variables
#   4. Independent K-means clustering using the prespecified six-variable
#      feature set and four-cluster solution
#   5. Assignment and characterization of the four cardiometabolic subtypes:
#      metabolically healthy (MH), hypertensive (HT), hyperglycemic (HG),
#      and renal impairment (RI)
#   6. Visualization of subtype distributions and metabolic profiles
#   7. Application of the multinomial ridge regression classifier developed
#      in the BIPASS discovery cohort
#   8. Evaluation of individual-level subtype assignment using one-versus-rest
#      receiver operating characteristic (ROC) analyses and area under the
#      curve (AUC) estimates
#
# Required input variables:
#   - Participant identifier
#   - Age
#   - Sex
#   - Body mass index (BMI)
#   - Systolic blood pressure (SBP)
#   - Diastolic blood pressure (DBP)
#   - Triglycerides (TG)
#   - Fasting blood glucose (FBG)
#   - Estimated glomerular filtration rate (eGFR)
#
# Input:
#   - A cohort-specific harmonized analytical dataset containing the required
#     variables listed above
#   - ridgeregression-bipass.rds
#       Multinomial ridge regression classifier, age- and sex-adjustment
#       models, and standardization parameters derived from the BIPASS
#       discovery cohort
#
# Output:
#   - Independently derived cluster assignments in the validation cohort
#   - CardioMet subtype assignments predicted by the BIPASS-derived classifier
#   - Descriptive summaries of clustering variables by subtype
#   - Subtype distribution and metabolic-profile visualizations
#   - One-versus-rest ROC curves and subtype-specific AUC estimates with
#     95% confidence intervals
#
# Notes:
#   - Cohort-specific raw-data extraction, variable construction, and data
#     harmonization are performed before this common validation pipeline and
#     are not included in this script.
#   - K-means clustering is fitted independently in each external cohort using
#     the same prespecified six variables and four-cluster specification as in
#     the BIPASS discovery cohort.
#   - Because K-means cluster numbering is arbitrary, independently derived
#     clusters are mapped to the four prespecified cardiometabolic subtype
#     labels according to their characteristic metabolic profiles.
#   - The BIPASS-derived ridge classifier is not refitted in external cohorts.
#     Age- and sex-adjustment models, standardization parameters, regression
#     coefficients, and the selected lambda value are inherited from BIPASS.
#   - Classifier performance is evaluated against subtype labels derived from
#     independently fitted K-means clustering in each external cohort.
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
load(file = "data_mean_imputed.Rdata")

cluster_vars_numeric <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI", "age")

# Convert sex coding in the new dataset from 0 = female / 1 = male
# to 1 = male / 2 = female to match the training dataset
data_mean_imputed$sex <- ifelse(
  data_mean_imputed$sex == 1, 1,
  ifelse(data_mean_imputed$sex == 0, 2, NA)
)

ukbdata <- data_mean_imputed

data_cluster <- data_mean_imputed[, c("eid",c(cluster_vars_numeric,"sex"))]
data_cluster <- na.omit(data_cluster)

id_col <- data_cluster[c("eid")]
data_values <- data_cluster[, c(cluster_vars_numeric,"sex")]

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
  age  = c(low = 15,   high = 95)
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
cluster_vars <- setdiff(colnames(data_values_clean), c("age", "sex"))
data_scaled <- data.frame(row.names = rownames(data_values_clean))

for (v in cluster_vars) {
  
  fit <- lm(
    formula = as.formula(paste(v, "~ age + sex")),
    data = data_values_clean
  )

  res <- resid(fit)
  res_scaled <- as.numeric(scale(res))
  data_scaled[[v]] <- res_scaled
}

str(data_scaled)

data_cluster_scaled <- cbind(id = id_col_clean, as.data.frame(data_scaled))

ukbdata <- subset(ukbdata,eid %in% data_cluster_scaled$eid)
data_cluster <- subset(data_cluster,eid %in% data_cluster_scaled$eid)#####

summary(data_cluster_scaled)

####---------------3. Distributional characteristics of candidate clustering variables-------------------####
library(dplyr)
library(gtsummary)

vars1 <- c(
  "age", "sex", "BMI", "SBP", "DBP",
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

####---------------4. Clustering analysis-------------------####
set.seed(20250814)
km_res <- kmeans(data_scaled, centers = 4, nstart = 25)

data_cluster$cluster <- factor(km_res$cluster)
ukbdata$clusterpamage <- factor(km_res$cluster)

recode_cluster <- function(x) {
  x_num <- as.numeric(as.character(x))
  x_new <- ifelse(
    x_num == 2, 3,
    ifelse(
      x_num == 3, 4,
      ifelse(
        x_num == 4, 2,
        x_num
      )
    )
  )
  return(x_new)
}

data_cluster$cluster <- recode_cluster(data_cluster$cluster)

ukbdata$clusterpamage <- recode_cluster(ukbdata$clusterpamage)

####---------------5. Pie chart-------------------####
library(ggplot2)
library(dplyr)

cluster_levels <- c(
  "MH",
  "HT",
  "HG",
  "RI"
)

cluster_data <- data.frame(
  Cluster = factor(cluster_levels, levels = cluster_levels),
  Count = as.vector(table(data_cluster$cluster))
)

cluster_data <- cluster_data %>%
  mutate(Percent = Count / sum(Count) * 100,
         Label = paste0(round(Percent), "%"))

custom_colors <- c(
  "#4E79A7",  # Hypertension-driven metabolic subtype
  "#B07AA1",  # Aging-related renal-impairment subtype
  "#E15759",  # Glucose-dysregulation subtype
  "#F28E2B"  # Metabolically healthy middle-aged group
)

library(ggplot2)
library(dplyr)

cluster_data <- cluster_data %>%
  arrange(Cluster) %>%
  mutate(
    Cluster = factor(Cluster, levels = unique(Cluster))
  )

cluster_data <- cluster_data %>%
  mutate(
    fraction = Count / sum(Count),
    ymax = cumsum(fraction),
    ymin = lag(ymax, default = 0),
    ymid = (ymax + ymin) / 2,
    label_y = (1 - ymid) * sum(Count),
    
    label_y_adj = case_when(
      Cluster == "MH" ~ label_y + 250,
      Cluster == "HG" ~ label_y + 120,
      Cluster == "HT" ~ label_y - 120,
      Cluster == "RI" ~ label_y - 250,
      TRUE ~ label_y
    ),
    
    text_x = case_when(
      Cluster %in% c("HG", "RI") ~ 1.75,
      TRUE ~ 1.75
    ),
    
    line_xend = case_when(
      Cluster %in% c("HG", "RI") ~ 1.75,
      TRUE ~ 1.75
    )
  )

p_cluster_pie <- ggplot(cluster_data, aes(x = 1, y = Count, fill = Cluster)) +
  
  geom_col(width = 1, color = "white", linewidth = 0.8) +
  
  coord_polar(theta = "y", start = 0, direction = -1, clip = "off") +
  
  geom_segment(
    aes(
      x = 1,
      xend = line_xend,
      y = label_y_adj,
      yend = label_y_adj
    ),
    color = "black",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      x = text_x,
      y = label_y_adj,
      label = Label
    ),
    color = "black",
    size = 4,
    hjust = 0
  ) +
  
  scale_fill_manual(
    values = custom_colors,
    breaks = levels(cluster_data$Cluster)
  ) +
  
  labs(
    title = "Validation cohort \n (UKB)",
    fill = NULL
  ) +
  
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    plot.margin  = margin(10, 80, 10, 10)
  )

p_cluster_pie

ggsave(
  filename = "Figure2_piechart_UKB.pdf",
  plot     = p_cluster_pie,
  device   = "pdf",
  width    = 6,
  height   = 3,
  units    = "in"
)

####---------------5. Box plot-------------------####
library(ggplot2)
library(dplyr)
library(patchwork)
library(rlang) 

vars <- c(
  "BMI",
  "SBP",
  "DBP",
  "FBG",
  "TG",
  "eGFR"
)

my_colors <- c(
  "1" = "#4E79A7",
  "2" = "#B07AA1",
  "3" = "#E15759",
  "4" = "#F28E2B"
)

data_cluster$cluster <- factor(
  data_cluster$cluster,
  levels = c("1","2","3","4")
)

ylim_list <- list(
  age        = c(20, 80),
  BMI        = c(15, 45),   
  SBP        = c(80, 200),  
  LDL        = c(0, 6),      
  eGFR       = c(10, 130),  
  FBG        = c(3, 15)     
)

var_titles <- c(
  BMI  = "BMI",
  SBP  = "SBP",
  DBP  = "DBP",
  FBG  = "FBG",
  TG   = "TG",
  eGFR = "eGFR"
)

var_units <- c(
  BMI  = "kg/m²",
  SBP  = "mmHg",
  DBP  = "mmHg",
  FBG  = "mmol/L",
  TG   = "mmol/L",
  eGFR = "mL/min/1.73 m²"
)

get_whisker_limits <- function(data, var, group) {
  
  df <- data %>%
    dplyr::select(!!sym(var), !!sym(group)) %>%
    dplyr::rename(value = 1, grp = 2)
  
  whiskers <- df %>%
    dplyr::group_by(grp) %>%
    dplyr::summarise(
      Q1  = quantile(value, 0.25, na.rm = TRUE),
      Q3  = quantile(value, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      lower = Q1 - 1.5 * IQR,
      upper = Q3 + 1.5 * IQR,
      .groups = "drop"
    )
  
  c(
    min(whiskers$lower, na.rm = TRUE),
    max(whiskers$upper, na.rm = TRUE)
  )
}

plot_list <- lapply(vars, function(var) {
  
  ylim_use <- get_whisker_limits(
    data_cluster,
    var,
    "cluster"
  )
  
  ggplot(
    data_cluster,
    aes(
      x = cluster,
      y = !!sym(var),
      fill = cluster
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      linewidth = 0.5
    ) +
    
    coord_cartesian(
      ylim = ylim_use
    ) +
    
    scale_fill_manual(
      values = my_colors
    ) +
    
    scale_x_discrete(
      labels = c(
        "1" = "MH",
        "2" = "HG",
        "3" = "HT",
        "4" = "RI"
      )
    ) +
    
    labs(
      title = var_titles[[var]],  
      x = NULL,
      y = var_units[[var]]     
    ) +
    
    theme_minimal(base_size = 11) +
    
    theme(
      
      panel.grid = element_blank(),
      
      axis.line = element_line(
        color = "black",
        linewidth = 0.6
      ),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.5
      ),
      
      axis.text = element_text(
        color = "black"
      ),
      
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      
      axis.title.y = element_text(
        size = 10
      ),
      
      legend.position = "none",
      
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
})

final_plot <- wrap_plots(
  plot_list,
  ncol = 6
)

final_plot

ggsave(
  filename = "Figure2_boxplot_UKB.pdf",
  plot     = final_plot,
  device   = cairo_pdf,  
  width    = 10,         
  height   = 4,        
  units    = "in"
)


# median output
library(dplyr)

vars <- c(
  "BMI",
  "SBP",
  "DBP",
  "FBG",
  "TG",
  "eGFR"
)

median_by_cluster <- data_cluster %>%
  group_by(cluster) %>%
  summarise(
    across(
      all_of(vars),
      ~ median(.x, na.rm = TRUE),
      .names = "{.col}_median"
    ),
    .groups = "drop"
  )

median_by_cluster

###############------------- 6. Radar chart  --------------###############

library(dplyr)
library(scales)
library(ggplot2)
library(ggradar)
library(patchwork)

#========================================================
# 1. Calculate cluster means
#========================================================

data_cluster_means <- data_cluster %>%
  group_by(cluster) %>%
  summarise(
    across(
      c(FBG, SBP, DBP, eGFR, TG, BMI),
      ~ mean(as.numeric(.), na.rm = TRUE)
    ),
    .groups = "drop"
  )

#========================================================
# 2. Rescale variables to 0–1
#    For eGFR, reverse the direction:
#    lower eGFR = higher renal impairment score
#========================================================

df_scaled <- data_cluster_means %>%
  mutate(across(-cluster, ~ as.numeric(.))) %>%
  mutate(
    FBG_scaled = rescale(FBG),
    SBP_scaled = rescale(SBP),
    DBP_scaled = rescale(DBP),
    TG_scaled  = rescale(TG),
    BMI_scaled = rescale(BMI),
    
    # Reverse eGFR direction
    # higher value = worse kidney function
    Reduced_eGFR = 1 - rescale(eGFR)
  ) %>%
  select(
    cluster,
    BMI = BMI_scaled,
    SBP = SBP_scaled,
    DBP = DBP_scaled,
    FBG = FBG_scaled,
    TG  = TG_scaled,
    `Reduced eGFR` = Reduced_eGFR
  ) %>%
  mutate(
    cluster = factor(cluster, levels = c("1", "2", "3", "4"))
  ) %>%
  arrange(cluster)

df_scaled <- as.data.frame(df_scaled)

#========================================================
# 3. Colors
#========================================================

my_colors <- c(
  "1" = "#4E79A7",
  "2" = "#B07AA1",
  "3" = "#E15759",
  "4" = "#F28E2B"
)

#========================================================
# 4. Radar plot
#========================================================

p_radar_all <- ggradar(
  df_scaled %>% rename(group = cluster),
  
  grid.min = 0,
  grid.mid = 0.5,
  grid.max = 1,
  
  group.line.width = 1.2,
  group.point.size = 3,
  
  background.circle.colour = "white",
  gridline.min.colour = "gray",
  gridline.mid.colour = "gray",
  gridline.max.colour = "gray",
  gridline.min.linetype = "solid",
  gridline.mid.linetype = "solid",
  gridline.max.linetype = "solid",
  
  grid.label.size = 0,
  
  group.colours = my_colors,
  fill = TRUE,
  fill.alpha = 0.25
) +
  labs(title = "") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none"
  )

print(p_radar_all)

ggsave(
  filename = "Figure2_radar_UKB.pdf",
  plot     = p_radar_all,
  width    = 2.7,
  height   = 2.7,
  units    = "in",
  dpi      = 300
)

####---------------7. Development of ridge regression model-------------------####
cluster_vars_numeric <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI")

library(glmnet)
library(pROC)
library(dplyr)
library(ggplot2)

model_obj <- readRDS(file = "ridgeregression-bipass.rds")

adjust_models <- model_obj$adjust_models
center <- model_obj$center
scale  <- model_obj$scale
vars   <- model_obj$vars

ridge_model <- model_obj$ridge_model
lambda      <- model_obj$lambda

compute_residual_matrix <- function(newdata, adjust_models, center, scale, vars) {
  
  newdata <- as.data.frame(newdata)
  
  n <- nrow(newdata)
  
  resid_matrix <- matrix(NA, nrow = n, ncol = length(vars))
  colnames(resid_matrix) <- vars
  
  for (v in vars) {
    
    fit <- adjust_models[[v]]
    coef_vec <- coef(fit)
    
    X <- model.matrix(~ sex + age, data = newdata)
    common_cols <- intersect(colnames(X), names(coef_vec))
    
    if (length(common_cols) == 0) {
      pred <- rep(coef_vec[1], n)
    } else {
      pred <- as.vector(
        X[, common_cols, drop = FALSE] %*%
          coef_vec[common_cols]
      )
    }
    resid_matrix[, v] <- newdata[[v]] - pred
  }
  
  resid_scaled <- scale(
    resid_matrix,
    center = center,
    scale  = scale
  )
  
  as.matrix(resid_scaled)
}

X_new <- compute_residual_matrix(
  newdata = data_cluster,
  adjust_models = adjust_models,
  center = center,
  scale = scale,
  vars = vars
)

prob_array <- predict(
  ridge_model,
  newx = X_new,
  s = lambda,
  type = "response"
)

prob_matrix <- prob_array[, , 1]
predicted_cluster <- colnames(prob_matrix)[max.col(prob_matrix)]
data_cluster$pred_cluster <- predicted_cluster
cluster_levels <- colnames(prob_matrix)
auc_results <- data.frame()

for (cl in cluster_levels) {
  
  true_label <- ifelse(data_cluster$cluster == cl, 1, 0)
  
  roc_obj <- roc(true_label, prob_matrix[, cl])
  
  auc_results <- rbind(
    auc_results,
    data.frame(
      cluster = cl,
      auc = auc(roc_obj),
      ci_low = ci.auc(roc_obj)[1],
      ci_high = ci.auc(roc_obj)[3]
    )
  )
}

auc_results

library(pROC)
library(dplyr)
library(ggplot2)

cluster_name_map <- c(
  "1" = "MH",
  "2" = "HT",
  "3" = "HG",
  "4" = "RI"
)

roc_list <- list()

for (cl in cluster_levels) {
  
  true_label <- ifelse(data_cluster$cluster == cl, 1, 0)
  
  roc_obj <- roc(true_label, prob_matrix[, cl], quiet = TRUE)
  
  roc_df <- data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    cluster = cluster_name_map[as.character(cl)]
  )
  
  roc_list[[cl]] <- roc_df
}

roc_plot_df <- bind_rows(roc_list)

auc_results <- auc_results %>%
  mutate(
    cluster = cluster_name_map[as.character(cluster)]
  )

mean_auc <- round(mean(auc_results$auc, na.rm = TRUE), 2)

legend_title <- paste0("AUC: ", sprintf("%.2f", mean_auc))

roc_plot_df$cluster <- factor(
  roc_plot_df$cluster,
  levels = c("MH", "HT", "HG", "RI")
)

my_colors <- c(
  "MH" = "#4E79A7",
  "HT" = "#B07AA1",
  "HG" = "#E15759",
  "RI" = "#F28E2B"
)

legend_title <- sprintf("AUC = %.2f", mean_auc)

p <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = cluster)) +
  
  geom_line(size = 1.3) +
  
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    color = "gray"
  ) +
  
  scale_color_manual(
    values = my_colors,
    name = legend_title
  ) +
  
  labs(
    x = "1 − Specificity",
    y = "Sensitivity"
  ) +
  
  coord_equal() +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = c(0.75, 0.25),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    axis.title = element_text(face = "bold")
  )

p

ggsave(
  filename = "Figure2-ROC_curve-UKB.pdf",
  plot = p,
  device = cairo_pdf,   # 🔥 推荐：嵌入字体更稳定（期刊友好）
  width = 5,            # 英寸（常用单栏宽度）
  height = 5,           # 保持正方形（你用了 coord_equal）
  units = "in",
  dpi = 600
)

# 输出auc具体数值
auc_results <- data.frame()

for (cl in cluster_levels) {
  
  true_label <- ifelse(data_cluster$cluster == cl, 1, 0)
  
  roc_obj <- roc(true_label, prob_matrix[, cl], quiet = TRUE)
  
  # AUC CI
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- ci.auc(roc_obj)
  
  # 最佳cutoff（Youden index）
  best_coords <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c(
      "threshold",
      "sensitivity",
      "specificity",
      "accuracy"
    )
  )
  
  auc_results <- rbind(
    auc_results,
    data.frame(
      
      cluster = cluster_name_map[as.character(cl)],
      
      auc = round(auc_val, 3),
      
      ci_low  = round(as.numeric(auc_ci[1]), 3),
      ci_high = round(as.numeric(auc_ci[3]), 3),
      
      cutoff = round(best_coords["threshold"], 3),
      sensitivity = round(best_coords["sensitivity"], 3),
      specificity = round(best_coords["specificity"], 3),
      accuracy = round(best_coords["accuracy"], 3)
    )
  )
}

auc_results

