library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(FactoMineR)
library(effsize)
library(reticulate)
library(mclust)
library(reshape2)
library(corrplot)
library(RColorBrewer)
library(circlize)

options(digits=3)
options(max.print=1000)

process_file <- function(file) {
  df <- read.csv(file)
  
  check_pvals <- function(pvals, stats, variable_name) {
    if (any(pvals < 0.05)) {
      cat(sprintf("Significant Results for %s:\n", variable_name))
      for (i in seq_along(pvals)) {
        cat(sprintf("Stat: %f, P-value: %f \n", stats[i], pvals[i]))
      }
    }
  }
  
  test_result_sex <- prop.test(table(df$cluster_labels_area, df$sex))
  test_result_Dx <- prop.test(table(df$Dx, df$cluster_labels_area))
  test_result_age <- wilcox.test(age_at_scan ~ cluster_labels_area, data = df)
  
  pvals_sex <- test_result_sex$p.value
  pvals_Dx <- test_result_Dx$p.value
  pvals_age <- test_result_age$p.value
  
  check_pvals(pvals_sex, test_result_sex$statistic, "Sex")
  check_pvals(pvals_Dx, test_result_Dx$statistic, "Dx")
  check_pvals(pvals_age, test_result_age$statistic, "Age")
}

process_file('df1.csv')
process_file('df2.csv')

compute_cohen_d <- function(df, cluster_label_col, start_col, end_col) {
  d_values <- sapply(start_col:end_col, function(i) {
    group1 <- df[df[[cluster_label_col]] == "cluster_1", i]
    group2 <- df[df[[cluster_label_col]] == "cluster_2", i]
    cohen.d(group1, group2, pooled = TRUE)$estimate
  })
  return(d_values)
}

df1 <- read.csv('df1.csv')
df2 <- read.csv('df2.csv')

df1_area_effect_size <- compute_cohen_d(df1, "cluster_labels_area", 22, 97)
df2_area_effect_size <- compute_cohen_d(df2, "cluster_labels_area", 22, 97)

write.csv(df1_area_effect_size, "df1_area_effect_size.csv", row.names = FALSE)
write.csv(df2_area_effect_size, "df2_area_effect_size.csv", row.names = FALSE)

age_stats <- df1 %>%
  group_by(Dx) %>%
  summarize(
    Age_Median = median(age_at_scan, na.rm = TRUE),
    Age_IQR = IQR(age_at_scan, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Lower_Bound = Age_Median - Age_IQR/2,
    Upper_Bound = Age_Median + Age_IQR/2
  )
print(age_stats)

cluster_labels <- c("cluster_labels_area", "cluster_labels_subcortical", "cluster_labels_thickness", "cluster_labels_volume")
aris_hbn <- outer(cluster_labels, cluster_labels, Vectorize(function(x, y) adjustedRandIndex(df2[[x]], df2[[y]])))
rownames(aris_hbn) <- cluster_labels
colnames(aris_hbn) <- cluster_labels

corrplot(aris_hbn, type="upper", order="hclust", col=rev(brewer.pal(n = 4, name = "RdYlBu")),
         tl.cex=0.6, tl.col="black", tl.srt=45, is.corr=FALSE)
