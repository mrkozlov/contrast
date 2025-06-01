Contrast Analysis for Species Distribution Modeling
This R script (contrast.R) performs a contrast analysis for species distribution modeling by generating pseudo-absences using Surface Range Envelope (SRE) and calculating statistical differences between presence and pseudo-absence points across environmental covariates. It leverages the sf, terra, stats, and writexl packages to process spatial data and export results.

Features
Pseudo-absence generation: Uses SRE to generate pseudo-absences at a 1:1 ratio with presence points, based on selected environmental layers.
Spatial processing: Crops raster layers to a study area polygon, ensuring consistent CRS (WGS84).
Statistical analysis: Computes means, standard deviations, standard errors, t-statistics, and p-values for environmental covariates, comparing presence (nf=1) and pseudo-absence (nf=0) points.
Output: Saves pseudo-absences as CSV and shapefile, and exports statistical results to an Excel file (or CSV if Excel fails).
Requirements
R (version 4.0 or higher recommended)
Required R packages:
R

install.packages(c("sf", "terra", "stats", "writexl"))
Input data:
train_balanced: An sf object with presence points (or a shapefile train_balanced.shp).
15_11Layers.shp: A polygon shapefile defining the study area.
Environmental rasters in t30/ directory (e.g., wc2.1_30s_bio_8.tif, wc2.1_30s_bio_14.tif, wc2.1_30s_bio_3.tif, etc.).
AFTER that split data by code:
# Подключение библиотек
library(FNN)      # Для KNN
library(ggplot2)  # Для визуализации
library(svglite)  # Для сохранения графиков
library(dplyr)    # Для манипуляций с данными
library(terra)    # Для работы с пространственными данными
library(sf)

# Simulated data setup (replace with actual combined_data)
# For demonstration, create a data frame with 970 observations and 21 variables
set.seed(42)
n <- 970
combined_data <- data.frame(
  wc2.1_30s_bio_1 = rnorm(n, mean = 10, sd = 5),
  wc2.1_30s_bio_10 = rnorm(n, mean = 20, sd = 5),
  wc2.1_30s_bio_11 = rnorm(n, mean = 5, sd = 5),
  wc2.1_30s_bio_12 = rnorm(n, mean = 600, sd = 100),
  wc2.1_30s_bio_13 = rnorm(n, mean = 100, sd = 30),
  wc2.1_30s_bio_14 = rnorm(n, mean = 30, sd = 15),
  wc2.1_30s_bio_15 = rnorm(n, mean = 30, sd = 15),
  wc2.1_30s_bio_16 = rnorm(n, mean = 200, sd = 50),
  wc2.1_30s_bio_17 = rnorm(n, mean = 100, sd = 30),
  wc2.1_30s_bio_18 = rnorm(n, mean = 150, sd = 50),
  wc2.1_30s_bio_19 = rnorm(n, mean = 200, sd = 50),
  wc2.1_30s_bio_2 = rnorm(n, mean = 8, sd = 2),
  wc2.1_30s_bio_3 = rnorm(n, mean = 30, sd = 10),
  wc2.1_30s_bio_4 = rnorm(n, mean = 500, sd = 200),
  wc2.1_30s_bio_5 = rnorm(n, mean = 25, sd = 5),
  wc2.1_30s_bio_6 = rnorm(n, mean = 0, sd = 5),
  wc2.1_30s_bio_7 = rnorm(n, mean = 25, sd = 5),
  wc2.1_30s_bio_8 = rnorm(n, mean = 15, sd = 5),
  wc2.1_30s_bio_9 = rnorm(n, mean = 10, sd = 5),
  wc2.1_30s_elev = rnorm(n, mean = 200, sd = 100),
  nf = sample(c(1, 100000000000), n, replace = TRUE, prob = c(0.9, 0.1))
)

# Simulate x and y coordinates (replace with actual coordinates if available)
combined_data$x <- runif(n, min = -180, max = 180)  # Simulated longitude
combined_data$y <- runif(n, min = -90, max = 90)    # Simulated latitude

# Parameters for splitting
params <- list(
  test_size = 0.1,            # 10% for test set (~97 points)
  max_attempts = 90000,       # Maximum attempts for optimization
  n_clusters = 12,            # Number of clusters
  w_cross_dist = 7.0,         # Weight for cross-group distance
  w_coverage = 9.0,           # Weight for cluster coverage
  w_balance = 60.0,           # Weight for balance
  w_dist_diff = 0.1,          # Weight for distance difference
  max_points_per_cluster = 2  # Maximum points per cluster in test set
)

# Function to compute spatial metrics
check_spatial_metrics <- function(data) {
  if (nrow(data) < 2) {
    return(list(mean_dist = 0, min_dist = 0, max_dist = 0))
  }
  
  coords <- as.matrix(data[, c("x", "y")])
  dist_matrix <- as.matrix(dist(coords))
  
  return(list(
    mean_dist = mean(dist_matrix[upper.tri(dist_matrix)]),
    min_dist = min(dist_matrix[upper.tri(dist_matrix)]),
    max_dist = max(dist_matrix[upper.tri(dist_matrix)])
  ))
}

# Function to compute cross-group distance
check_cross_distance <- function(train_data, test_data) {
  if (nrow(train_data) == 0 || nrow(test_data) == 0) {
    return(0)
  }
  
  knn <- get.knnx(train_data[, c("x", "y")], 
                  test_data[, c("x", "y")], 
                  k = 1)
  return(mean(knn$nn.dist))
}

# Function for balanced cluster-based splitting
balanced_cluster_split <- function(pts, params) {
  # Check for x and y columns
  if (!all(c("x", "y") %in% colnames(pts))) {
    stop("Данные должны содержать колонки 'x' и 'y'")
  }
  
  # Check numeric type
  if (!is.numeric(pts$x) || !is.numeric(pts$y)) {
    stop("Столбцы 'x' и 'y' должны быть числовыми")
  }
  
  # Clustering
  set.seed(42)
  coords <- as.matrix(pts[, c("x", "y")])
  clusters <- kmeans(coords, centers = params$n_clusters)$cluster
  pts$cluster <- as.factor(clusters)
  
  # Print cluster distribution
  cat("\nРаспределение всех точек по кластерам:\n")
  print(table(pts$cluster))
  
  # Target test points per cluster
  target_test <- round(table(pts$cluster) * params$test_size)
  target_test[target_test == 0] <- 1  # Minimum 1 point per cluster
  
  # Available clusters
  cluster_counts <- table(pts$cluster)
  available_clusters <- names(cluster_counts[cluster_counts > 0])
  
  best_split <- NULL
  best_score <- -Inf
  
  for (attempt in 1:params$max_attempts) {
    # Sample points from each cluster
    test_indices <- unlist(lapply(names(target_test), function(cl) {
      cluster_pts <- which(pts$cluster == cl)
      sample(cluster_pts, min(target_test[cl], length(cluster_pts)), replace = FALSE)
    }))
    
    # Ensure cluster coverage
    temp_test_data <- pts[test_indices, ]
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(pts$cluster == cl & !(1:nrow(pts) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    
    # Limit points per cluster
    temp_test_data <- pts[test_indices, ]
    cluster_counts_test <- table(temp_test_data$cluster)
    for (cl in names(cluster_counts_test)) {
      if (cluster_counts_test[cl] > params$max_points_per_cluster) {
        cluster_indices <- which(temp_test_data$cluster == cl)
        excess <- cluster_counts_test[cl] - params$max_points_per_cluster
        remove_indices <- sample(cluster_indices, excess)
        test_indices <- test_indices[!test_indices %in% temp_test_data$index[remove_indices]]
      }
    }
    
    # Add index for tracking
    pts$index <- 1:nrow(pts)
    temp_test_data <- pts[test_indices, ]
    
    # Ensure coverage after limiting
    cluster_coverage <- table(temp_test_data$cluster)
    if (!all(available_clusters %in% names(cluster_coverage))) {
      missing_clusters <- available_clusters[!available_clusters %in% names(cluster_coverage)]
      for (cl in missing_clusters) {
        cluster_pts <- which(pts$cluster == cl & !(1:nrow(pts) %in% test_indices))
        if (length(cluster_pts) > 0) {
          test_indices <- c(test_indices, sample(cluster_pts, 1))
        }
      }
    }
    
    # Create train and test sets
    train_data <- pts[-test_indices, ]
    test_data <- pts[test_indices, ]
    
    # Remove temporary index
    train_data$index <- NULL
    test_data$index <- NULL
    pts$index <- NULL
    
    # Compute metrics
    train_metrics <- check_spatial_metrics(train_data)
    test_metrics <- check_spatial_metrics(test_data)
    cross_dist <- check_cross_distance(train_data, test_data)
    coverage <- mean(table(test_data$cluster) > 0)
    balance <- 1 - sd(table(test_data$cluster)) / mean(table(test_data$cluster))
    dist_diff <- abs(train_metrics$mean_dist - test_metrics$mean_dist)
    
    # Combined score
    score <- (cross_dist * params$w_cross_dist + 
                coverage * params$w_coverage + 
                balance * params$w_balance - 
                dist_diff * params$w_dist_diff)
    
    if (score > best_score) {
      best_split <- list(train = train_data, test = test_data)
      best_score <- score
    }
  }
  
  return(best_split)
}

# Extent calculation function
make_extent_from_pts <- function(pts, digits = 1) {
  if (!all(c("x", "y") %in% names(pts))) {
    stop("Датафрейм должен содержать колонки 'x' и 'y'")
  }
  extent(
    round(min(pts$x), digits),
    round(max(pts$x), digits),
    round(min(pts$y), digits),
    round(max(pts$y), digits)
  )
}

# Apply splitting
pts <- combined_data
set.seed(42)
split_result <- balanced_cluster_split(pts, params)
train_data <- split_result$train
test_data <- split_result$test

# Compute extent
e <- make_extent_from_pts(pts, digits = 1)
cat("\nЭкстент точек:\n")
print(e)

# Compute metrics
train_metrics <- check_spatial_metrics(train_data)
test_metrics <- check_spatial_metrics(test_data)
cross_dist <- check_cross_distance(train_data, test_data)
coverage <- length(unique(test_data$cluster)) / params$n_clusters

# Visualization
knn <- get.knnx(train_data[, c("x", "y")], test_data[, c("x", "y")], k = 1)
segment_data <- data.frame(
  x = test_data$x,
  y = test_data$y,
  xend = train_data$x[knn$nn.index],
  yend = train_data$y[knn$nn.index]
)

final_plot <- ggplot() +
  geom_point(data = pts, aes(x, y), color = "grey90", size = 2) +
  geom_point(data = train_data, aes(x, y, color = "Training"), size = 3, alpha = 0.7) +
  geom_point(data = test_data, aes(x, y, color = "Testing"), shape = 17, size = 4) +
  geom_segment(
    data = segment_data,
    aes(x = x, y = y, xend = xend, yend = yend),
    color = "darkgrey", linetype = "dashed", alpha = 0.5
  ) +
  scale_color_manual(
    name = "Dataset",
    values = c("Training" = "#3575b5", "Testing" = "#e74c3c")
  ) +
  labs(
    title = "Balanced Cluster-Based Data Split",
    subtitle = sprintf(
      "Training: %d points | Testing: %d points | Cross-dist: %.2f | Cluster coverage: %.0f%%",
      nrow(train_data), nrow(test_data), cross_dist, coverage * 100
    ),
    x = "X coordinate",
    y = "Y coordinate"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Save plot
ggsave("balanced_cluster_split_plot.png", final_plot, width = 8, height = 6, dpi = 300)
ggsave("balanced_cluster_split_plot.tiff", final_plot, width = 8, height = 6, dpi = 300, compression = "lzw")
ggsave("balanced_cluster_split_plot.pdf", final_plot, width = 8, height = 6, dpi = 300)

# Save data
write.table(test_data, file = "testAria3.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(train_data, file = "trainAria3.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Generate report
report_text <- paste(
  "=== SPLITTING RESULTS ===",
  "\nParameters:",
  sprintf("- Test size: %.1f", params$test_size),
  sprintf("- Max attempts: %d", params$max_attempts),
  sprintf("- Number of clusters: %d", params$n_clusters),
  sprintf("- Cross-distance weight: %.1f", params$w_cross_dist),
  sprintf("- Coverage weight: %.1f", params$w_coverage),
  sprintf("- Balance weight: %.1f", params$w_balance),
  sprintf("- Distance difference weight: %.1f", params$w_dist_diff),
  sprintf("- Max points per cluster: %d", params$max_points_per_cluster),
  "\nDataset sizes:",
  sprintf("- Training set: %d points (%.1f%%)", nrow(train_data), nrow(train_data)/nrow(pts)*100),
  sprintf("- Testing set: %d points (%.1f%%)", nrow(test_data), nrow(test_data)/nrow(pts)*100),
  "\nSpatial metrics:",
  sprintf("- Mean distance in train: %.3f", train_metrics$mean_dist),
  sprintf("- Mean distance in test: %.3f", test_metrics$mean_dist),
  sprintf("- Cross-group distance: %.3f", cross_dist),
  sprintf("- Cluster coverage: %.1f%%", coverage * 100),
  "\nTesting points per cluster:",
  capture.output(print(table(test_data$cluster))),
  sep = "\n"
)

# Save report
writeLines(report_text, "splitting_report.txt")

# Output plot and report
print(final_plot)
cat(report_text)