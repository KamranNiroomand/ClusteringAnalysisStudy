library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(factoextra)
library(FactoMineR)
options(digits=3)
options(max.print=1000)

PCA_dims <- function(data, variance_threshold = 0.75) {
  pca_result <- PCA(data, ncp = dim(data)[2], graph = FALSE)
  explained_variance <- pca_result$eig[1, ] / sum(pca_result$eig)
  num_components <- sum(cumsum(explained_variance) < variance_threshold) + 1
  full_pcas = pca_result$ind$coord
  thresholded_pcas = pca_result$ind$cord[,1:num_components]
  full_loads = pca_result$var$cor
  return(list(full_pcas = full_pcas, thresholded_pcas = thresholded_pcas, full_loads = full_loads))
}

cosine_similarity_matrix <- function(data) {
  normalized_data <- scale(data, center = FALSE, scale = TRUE)
  similarity_matrix <- normalized_data %*% t(normalized_data)
  similarity_matrix[similarity_matrix > 1] <- 1
  similarity_matrix[similarity_matrix < -1] <- -1
  return(similarity_matrix)
}

data <- read_csv("data.csv")
data$sex <- as.numeric(replace(data$sex, data$sex == 'Male', 0))
data$sex <- as.numeric(replace(data$sex, data$sex == 'Female', 1))

data_POND <- filter(data, Dataset == "POND")
data_HBN <- filter(data, Dataset == "HBN")

area_POND <- dplyr::select(data_POND,'lobeArea40mm_Left_Gyrus_Rectus':'lobeArea40mm_Right_Insula')
vol_POND <- dplyr::select(data_POND,'lobeVolume_Left_Gyrus_Rectus':'lobeVolume_Right_Insula')
thickness_POND <- dplyr::select(data_POND,'lobeThickness_Left_Gyrus_Rectus':'lobeThickness_Right_Insula')
subvol_POND <- dplyr::select(data_POND,'leftNucleusAccumbens_VentralStriatum':'rightThalamus')

area_HBN <- dplyr::select(data_HBN,'lobeArea40mm_Left_Gyrus_Rectus':'lobeArea40mm_Right_Insula')
vol_HBN <- dplyr::select(data_HBN,'lobeVolume_Left_Gyrus_Rectus':'lobeVolume_Right_Insula')
thickness_HBN <- dplyr::select(data_HBN,'lobeThickness_Left_Gyrus_Rectus':'lobeThickness_Right_Insula')
subvol_HBN <- dplyr::select(data_HBN,'leftNucleusAccumbens_VentralStriatum':'rightThalamus')

area_POND_scaled <- scale(area_POND, center=TRUE, scale=TRUE)
vol_POND_scaled <- scale(vol_POND, center=TRUE, scale=TRUE)
thickness_POND_scaled <- scale(thickness_POND, center=TRUE, scale=TRUE)
subvol_POND_scaled <- scale(subvol_POND, center=TRUE, scale=TRUE)

area_HBN_scaled <- scale(area_HBN, center=TRUE, scale=TRUE)
vol_HBN_scaled <- scale(vol_HBN, center=TRUE, scale=TRUE)
thickness_HBN_scaled <- scale(thickness_HBN, center=TRUE, scale=TRUE)
subvol_HBN_scaled <- scale(subvol_HBN, center=TRUE, scale=TRUE)

pca_area_pond <- PCA_dims(area_POND_scaled, variance_threshold = 0.75)$full_pcas
pca_volume_pond <- PCA_dims(vol_POND_scaled, variance_threshold = 0.75)$full_pcas
pca_thickness_pond <- PCA_dims(thickness_POND_scaled, variance_threshold = 0.75)$full_pcas
pca_subcort_pond <- PCA_dims(subvol_POND_scaled, variance_threshold = 0.75)$full_pcas

pca_area_hbn <- PCA_dims(area_HBN_scaled, variance_threshold = 0.75)$full_pcas
pca_volume_hbn <- PCA_dims(vol_HBN_scaled, variance_threshold = 0.75)$full_pcas
pca_thickness_hbn <- PCA_dims(thickness_HBN_scaled, variance_threshold = 0.75)$full_pcas
pca_subcort_hbn <- PCA_dims(subvol_HBN_scaled, variance_threshold = 0.75)$full_pcas

write.csv(as.data.frame(pca_area_pond), "pca_area_pond.csv", row.names = FALSE)
write.csv(as.data.frame(pca_volume_pond), "pca_volume_pond.csv", row.names = FALSE)
write.csv(as.data.frame(pca_thickness_pond), "pca_thickness_pond.csv", row.names = FALSE)
write.csv(as.data.frame(pca_subcort_pond), "pca_subcort_pond.csv", row.names = FALSE)
write.csv(as.data.frame(pca_area_hbn), "pca_area_hbn.csv", row.names = FALSE)
write.csv(as.data.frame(pca_volume_hbn), "pca_volume_hbn.csv", row.names = FALSE)
write.csv(as.data.frame(pca_thickness_hbn), "pca_thickness_hbn.csv", row.names = FALSE)
write.csv(as.data.frame(pca_subcort_hbn), "pca_subcort_hbn.csv", row.names = FALSE)

Dimension <- 3
Dimnesionx <- 3
load_data <- PCA_dims(subvol_POND_scaled, variance_threshold = 0.75)$full_loads
load_datax <- PCA_dims(subvol_HBN_scaled, variance_threshold = 0.75)$full_loads
cor_test <- cor.test(load_data[,Dimension], load_datax[,Dimnesionx])
cor_value <- cor_test$estimate
p_value <- cor_test$p.value

plot_data <- data.frame(POND = load_data[,Dimension], HBN = load_datax[,Dimnesionx])

p <- ggplot(plot_data, aes(x=POND, y=HBN)) +
  geom_point(shape = 16, size = 3, color = 'darkgreen') +
  geom_smooth(method='lm', se=FALSE, color='black', linetype="dashed") + 
  labs(
    title = "Correlation Plot Between Same Dimension Loading in POND and HBN",
    x = paste("POND Loading for Component", Dimension),
    y = paste("HBN Loading for Component", Dimnesionx),
    caption = paste("Correlation is", round(cor_value, digits=3), "& p-value is", format.pval(p_value, digits=3))
  ) +
  theme_minimal(base_family = "Roboto") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    plot.caption = element_text(hjust = 0.5, size = 10),
    legend.position = "none"
  )

print(p)
