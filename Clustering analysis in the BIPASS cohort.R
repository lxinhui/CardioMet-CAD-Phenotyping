###############################################################################
# CardioMet: Clustering analysis in the BIPASS cohort
#
# Study:
# Unsupervised cardiometabolic phenotyping reveals prognostic and 
# treatment-associated heterogeneity in coronary artery disease
#
# Description:
# This script performs the complete clustering workflow in the BIPASS
# discovery cohort, including:
#   1. Data preparation and quality control
#   2. Exclusion of clinically implausible outliers
#   3. Age- and sex-adjustment of candidate metabolic variables
#   4. Selection of the optimal number of clusters
#   5. K-means clustering and cluster stability assessment
#   6. Visualization of cardiometabolic subtypes
#      (pie chart, UMAP, boxplots, radar chart, heatmap)
#   7. Development of the multinomial ridge regression classifier for
#      subtype assignment in external cohorts
#
# Input:
#   - Cleaned BIPASS dataset
#
# Output:
#   - Cluster assignments
#   - Figures for the manuscript
#   - Trained multinomial ridge regression model (.rds)
#
# Author: Xinhui Liu
# Institution: Department of Emergency Medicine,
#              Qilu Hospital, Shandong University
#
# Correspondence:
# Jiali Wang, MD, PhD
# Yuguo Chen, MD, PhD
#
# Repository:
# https://github.com/lxinhui/CardioMet
#
###############################################################################

###############------------- 1. data prepare --------------###############
###############------------- 2. Clustering analysis  --------------###############
library(cluster)
library(factoextra)
library(fpc)
library(tidyverse)

vars <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI", "age")
cluster_vars_numeric <- vars

# Exclusion of extreme outliers
data_cluster <- data_mean_imputed[, c("id",c(cluster_vars_numeric,"sex"))]
data_cluster <- na.omit(data_cluster)

id_col <- data_cluster["id"]   
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


datafix <- subset(datafix,id %in% data_cluster_scaled$id)
ukbdata <- subset(ukbdata,id %in% data_cluster_scaled$id)
datafixfix <- subset(datafixfix,id %in% data_cluster_scaled$id)#####
data_cluster <- subset(data_cluster,id %in% data_cluster_scaled$id)#####

###############------------- 3. Distributional characteristics of candidate clustering variables  --------------###############
library(dplyr)
library(gtsummary)

vars <- c(
  "age", "sex", "BMI", "SBP", "DBP",
  "HDL", "LDL", "TC", "TG", "FBG", "eGFR"
)

table1 <- ukbdata %>%
  select(all_of(vars)) %>%
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

###############------------- 4. Selection of the optimal number of clusters  --------------###############
library(cluster)
library(fpc)

data_use <- data_scaled
sil_result <- data.frame(k = 3:8, silhouette = NA, jaccard = NA)
d <- dist(data_use)

for (k in 3:8) {

  cat("Running k =", k, "...\n")

  sil_tmp <- numeric(30)

  for (i in 1:30) {
    km_fit <- kmeans(data_use, centers = k, nstart = 25)
    ss <- silhouette(km_fit$cluster, d)
    sil_tmp[i] <- mean(ss[,3])
  }

  mean_sil <- mean(sil_tmp)

  ## Jaccard stability
  cb <- clusterboot(
    data_use,
    B = 50,
    bootmethod = "boot",
    clustermethod = kmeansCBI,
    krange = k,
    seed = 20260127
  )

  mean_jaccard <- mean(cb$bootmean)

  sil_result[sil_result$k == k, "silhouette"] <- mean_sil
  sil_result[sil_result$k == k, "jaccard"] <- mean_jaccard
}

sil_result


###############------------- 5. Clustering analysis  --------------###############
set.seed(20250814)

km_res <- kmeans(data_scaled, centers = 4, nstart = 25)
# saveRDS(km_res, file = "20251017-kmeans_model-fbg-bipass.rds")

data_cluster$cluster <- factor(km_res$cluster)
datafixfix$clusterpamage <- factor(km_res$cluster)
ukbdata$clusterpamage <- factor(km_res$cluster)

# Cluster-specific Jaccard similarity indices
library(fpc)

set.seed(20250814)

cb <- clusterboot(
  data_scaled,
  B = 200,
  bootmethod = "boot",
  clustermethod = kmeansCBI,
  krange = 4,
  seed = 20260127
)

cb$bootmean

# Cluster labels were reassigned according to the plotting order for ease of interpretation.
recode_cluster <- function(x) {
  dplyr::recode(x,
                `1` = 4,
                `2` = 2,
                `3` = 1,
                `4` = 3)
}

data_cluster$cluster <- recode_cluster(data_cluster$cluster)
datafixfix$clusterpamage <- recode_cluster(datafixfix$clusterpamage)
ukbdata$clusterpamage    <- recode_cluster(ukbdata$clusterpamage)


###############------------- 6. Pie chart  --------------###############
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
  mutate(
    Percent = Count / sum(Count) * 100,
    Label = paste0(round(Percent), "%")
  )

custom_colors <- c(
  "#4E79A7", 
  "#B07AA1", 
  "#E15759", 
  "#F28E2B" 
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
    title = "Discovery cohort \n (BIPass)",
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
  filename = "Figure2_piechart_BIPass.pdf",
  plot     = p_cluster_pie,
  device   = "pdf",
  width    = 6,
  height   = 3,
  units    = "in"
)


###############------------- 7. UMAP  --------------###############
library(umap)
library(ggplot2)
library(dplyr)
library(ggrepel)

vars <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI", "age")
cluster_vars_numeric <- vars
vars <- vars[-which(vars == "age")]

data_umap <- data_cluster %>%
  select(all_of(vars), cluster) %>%
  na.omit()

data_scaled <- scale(data_umap[, vars])

set.seed(20250814)  

umap_res <- umap(data_scaled, 
                 n_neighbors = 20,     
                 min_dist = 0.1,       
                 metric = "euclidean") 

plot_df <- data.frame(
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  
  cluster = factor(
    data_umap$cluster,
    levels = c(1, 2, 3, 4),
    labels = c("MH","HT","HG","RI"
    )
  )
)

my_colors <- c(
  "MH" = "#4E79A7",
  "HT" = "#B07AA1",
  "HG" = "#E15759",
  "RI" = "#F28E2B"
)

p_umap <- ggplot(
  plot_df,
  aes(
    x = UMAP1,
    y = UMAP2,
    color = cluster
  )
) +
  
  geom_point(
    alpha = 0.8,
    size = 0.8
  ) +
  
  stat_ellipse(
    aes(fill = cluster),
    geom = "polygon",
    alpha = 0.15,
    color = NA
  ) +
  
  scale_color_manual(
    values = my_colors
  ) +
  
  scale_fill_manual(
    values = my_colors
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    
    panel.background = element_rect(
      fill = "white",
      color = "black",
      linewidth = 0.8
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
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
      size = 14
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 11
    )
  ) +
  
  labs(
    title = "",
    
    x = "Dimension 1",
    y = "Dimension 2",
    
    color = "Cluster",
    fill  = "Cluster"
  )

print(p_umap)

ggsave(
  filename = "Figure 2-A_UMAP_clusters.pdf",
  plot     = p_umap,
  width    = 6,
  height   = 4.5,
  units    = "in"
)


###############------------- 8. Box plot  --------------###############
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
  filename = "Figure2-A_boxplot_bipass.pdf",
  plot     = final_plot,
  device   = cairo_pdf,  
  width    = 10,        
  height   = 3.2,         
  units    = "in"
)


###############------------- 9. Radar chart  --------------###############
library(dplyr)
library(scales)
library(ggplot2)
library(ggradar)
library(patchwork)

data_cluster_means <- data_cluster %>%
  group_by(cluster) %>%
  summarise(across(c(FBG,SBP,DBP,eGFR,TG,BMI),
                   ~ mean(as.numeric(.), na.rm = TRUE))) %>%
  ungroup()

df_scaled <- data_cluster_means %>%
  mutate(across(-cluster, ~ as.numeric(.))) %>%
  mutate(across(-cluster, ~ rescale(.))) %>%
  mutate(cluster = as.factor(cluster))

df_scaled <- as.data.frame(df_scaled)

my_colors <- c(
  "1" = "#4E79A7",
  "2" = "#B07AA1",
  "3" = "#E15759",
  "4" = "#F28E2B"
)

cluster_list <- split(df_scaled, df_scaled$cluster)

df_scaled_plot <- df_scaled %>%
  mutate(
    cluster = factor(cluster, levels = c("1", "2", "3", "4"))
  ) %>%
  arrange(cluster)

p_radar_all <- ggradar(
  df_scaled_plot %>% rename(group = cluster),
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
  filename = "Figure2_radar_bipass.pdf",
  plot     = p_radar_all,
  width    = 2.7,
  height   = 2.7,
  units    = "in",
  dpi      = 300
)


###############------------- 9. Radar chart  --------------###############

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
  filename = "Figure2_radar_bipass.pdf",
  plot     = p_radar_all,
  width    = 2.7,
  height   = 2.7,
  units    = "in",
  dpi      = 300
)

###############------------- 10. Heatmap of pairwise comparisons of clustering variables across CAD subtypes  --------------###############
library(dplyr)
library(tidyr)
library(ggplot2)
library(rstatix)

data_cluster$cluster <- factor(data_cluster$cluster)

cluster_pairs <- combn(levels(data_cluster$cluster), 2, simplify = FALSE)

# Cohen’s d + p value
results <- lapply(vars, function(var) {
  lapply(cluster_pairs, function(pair) {
    g1 <- data_cluster %>% filter(cluster == pair[1]) %>% pull(!!sym(var))
    g2 <- data_cluster %>% filter(cluster == pair[2]) %>% pull(!!sym(var))
    
    mean_diff <- median(g1, na.rm = TRUE) - median(g2, na.rm = TRUE)
    pooled_sd <- sqrt((var(g1, na.rm = TRUE) + var(g2, na.rm = TRUE)) / 2)
    d_value <- mean_diff / pooled_sd
    p_val <- wilcox.test(g1, g2)$p.value
    
    data.frame(
      Variable = var,
      Group1 = pair[1],
      Group2 = pair[2],
      d = d_value,
      p = p_val
    )
  }) %>% bind_rows()
}) %>% bind_rows()

cluster_labels <- c("1" = "MH", "2" = "HT", "3" = "HG", "4" = "RI")

results <- results %>%
  mutate(
    Group1_label = cluster_labels[as.character(Group1)],
    Group2_label = cluster_labels[as.character(Group2)],
    Comparison = paste0(Group1_label, " vs ", Group2_label),
    Sig = ifelse(p < 0.0084, "*", "")
  )

results <- results %>%
  mutate(
    d_cap = pmax(pmin(d, 0.99), -0.99)
  )

breaks_d <- c(-0.8, -0.5, -0.2, 0, 0.2, 0.5, 0.8)
labels_d <- c(
  "Large ≤ -0.80",
  "Medium -0.79~-0.50",
  "Small -0.49~-0.20",
  "Negligible -0.19~0.19",
  "Small 0.20~0.49",
  "Medium 0.50~0.79",
  "Large ≥ 0.80"
)

my_colors <- c(
  "#B6D4E3", "#DDEAF2", "#F7F7F7", "#FBE3D1", "#F4A582", "#CA0020"
)

distance <- ggplot(results, aes(x = Variable, y = Comparison, fill = d_cap)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Sig), color = "black", size = 4, fontface = "bold") +
  scale_fill_gradientn(
    colors  = my_colors,
    breaks  = breaks_d,
    labels  = labels_d,
    limits  = c(-0.99, 0.99),
    name    = "Cohen’s d"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "",
    subtitle = ""
  ) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid    = element_blank(),
    axis.text.x   = element_text(angle = 45, hjust = 1, face = "bold", color = "black"),
    axis.text.y   = element_text(face = "bold", color = "black"),
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    legend.text     = element_text(size = 9)
  )
distance

ggsave(
  filename = "Cluster_pairwise_effectsize_heatmap.pdf",
  plot     = distance,
  width    = 7,     # 根据变量个数可调整
  height   = 3,      # 根据 cluster 数量可调整
  units    = "in"
)



###############------------- 11. Development of ridge regression model  --------------###############
library(glmnet)

df <- data_cluster
vars <- c("FBG", "SBP", "DBP", "eGFR", "TG", "BMI")

adjust_models <- list()
residual_matrix <- matrix(NA, nrow = nrow(df), ncol = length(vars))
colnames(residual_matrix) <- vars

train_center <- numeric(length(vars))
train_scale  <- numeric(length(vars))

names(train_center) <- vars
names(train_scale)  <- vars

for (v in vars) {
  formula_str <- paste0(v, " ~ sex + age")
  fit <- lm(as.formula(formula_str), data = df)
  adjust_models[[v]] <- fit
  resid_v <- resid(fit)
  residual_matrix[, v] <- resid_v
  train_center[v] <- mean(resid_v, na.rm = TRUE)
  train_scale[v]  <- sd(resid_v, na.rm = TRUE)
}

X_train <- scale(residual_matrix,
                 center = train_center,
                 scale  = train_scale)

X_train <- as.matrix(X_train)

y_train <- factor(df$cluster) 

w <- rep(1, nrow(df))
w[y_train == "3"] <- 5   
w[y_train == "4"] <- 5   

cv_ridge <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "multinomial",
  alpha = 0,
  weights = w,
  standardize = FALSE,
  nfolds = 5
)

best_lambda <- cv_ridge$lambda.min
best_lambda

ridge_cluster_model <- list(
  
  ridge_model   = cv_ridge,
  lambda        = best_lambda,
  
  adjust_models = adjust_models,  
  
  center        = train_center,
  scale         = train_scale,
  
  vars          = vars,
  
  sex_levels    = levels(df$sex),  
  
  age_range     = range(df$age, na.rm = TRUE)
)

saveRDS(
  ridge_cluster_model,
  file = "ridgeregression-bipass.rds"
)
