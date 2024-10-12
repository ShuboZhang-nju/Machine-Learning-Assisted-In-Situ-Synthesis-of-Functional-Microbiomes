library(gbm)
library(dplyr)
library(ggplot2)
library(caret)
library(rBayesianOptimization)
library(pdp)  # For ICE analysis
library(plotly)

# Load data 
GBRT <- read.csv("GBRT.csv")

# Ensure all variables are numeric
GBRT[] <- lapply(GBRT, function(x) as.numeric(as.character(x)))

# Split data into training and testing sets
set.seed(1)
train_index <- sample(nrow(GBRT), 0.75 * nrow(GBRT))
GBRT_train <- GBRT[train_index, ]
GBRT_test <- GBRT[-train_index, ]

# Define function for hyperparameter tuning
optimize_gbm <- function(interaction.depth, n.trees, shrinkage, bag.fraction) {
  tryCatch({
    # Print current parameters
    cat("interaction.depth:", interaction.depth, "\n")
    cat("n.trees:", n.trees, "\n")
    cat("shrinkage:", shrinkage, "\n")
    cat("bag.fraction:", bag.fraction, "\n")
    
    # Train GBRT model
    model <- gbm(CSI ~ TN + NH4 + pH + DO + T + SRT + HRT + TP + TOC,
                 data = GBRT_train, distribution = "gaussian",
                 n.trees = as.integer(n.trees),
                 interaction.depth = as.integer(interaction.depth),
                 shrinkage = shrinkage,
                 bag.fraction = bag.fraction,
                 train.fraction = 0.75, cv.folds = 5, verbose = FALSE)
    
    # Use cross-validation to select the best iteration
    best_iter <- gbm.perf(model, method = "cv", plot.it = FALSE)
    pred <- predict(model, GBRT_test, n.trees = best_iter)
    rmse <- RMSE(pred, GBRT_test$CSI)
    list(Score = -rmse, Pred = pred)
  }, error = function(e) {
    print(e)
    return(list(Score = NA, Pred = NA))
  })
}

# Set hyperparameter bounds and perform Bayesian optimization (optional)
bounds <- list(
  interaction.depth = c(3L, 15L),
  n.trees = c(1000L, 5000L),
  shrinkage = c(0.001, 0.1),
  bag.fraction = c(0.5, 0.9)
)

set.seed(123)
opt_results <- BayesianOptimization(
  FUN = optimize_gbm,
  bounds = bounds,
  init_points = 10,
  n_iter = 20,
  acq = "ucb",
  kappa = 2.5,
  eps = 0.0,
  verbose = TRUE
)

# Extract the best parameters and convert to appropriate types
best_params <- as.list(opt_results$Best_Par)
n.trees <- as.integer(best_params["n.trees"])
interaction.depth <- as.integer(best_params["interaction.depth"])
shrinkage <- as.numeric(best_params["shrinkage"])
bag.fraction <- as.numeric(best_params["bag.fraction"])

# Train the final model with the best parameters
best_model <- gbm(CSI ~ TN + NH4 + pH + DO + T + SRT + HRT + TP + TOC,
                  data = GBRT_train, distribution = "gaussian",
                  n.trees = n.trees,
                  interaction.depth = interaction.depth,
                  shrinkage = shrinkage,
                  bag.fraction = bag.fraction,
                  train.fraction = 0.75, cv.folds = 5, verbose = FALSE)

# Select the best iteration
best_iter <- gbm.perf(best_model, method = "cv", plot.it = FALSE)

# Predict on the test set
GBRT_test$Ynpk_pred_1 <- predict(best_model, GBRT_test, n.trees = best_iter)

# Plot scatter plot of predicted vs. observed values
plot(GBRT_test$Ynpk_pred_1, GBRT_test$CSI,
     xlab = "Predicted CSI", ylab = "Observed CSI",
     main = "Predicted vs Observed",
     xlim = c(min(GBRT_test$CSI), max(GBRT_test$CSI)),
     ylim = c(min(GBRT_test$CSI), max(GBRT_test$CSI)))
abline(a = 0, b = 1, col = "red")
fitline <- lm(Degree ~ Ynpk_pred_1, data = GBRT_test)
abline(fitline, lty = 2)
summary(fitline)

# Calculate model evaluation metrics
r2_test <- summary(lm(CSI ~ Ynpk_pred_1, data = GBRT_test))$r.squared
rmse_test <- RMSE(GBRT_test$Ynpk_pred_1, GBRT_test$CSI)
mae_test <- MAE(GBRT_test$Ynpk_pred_1, GBRT_test$CSI)
mape_test <- mean(abs((GBRT_test$CSI - GBRT_test$Ynpk_pred_1) / GBRT_test$Degree)) * 100

model_evaluation <- data.frame(
  Dataset = c("Test"),
  R2 = c(r2_test),
  RMSE = c(rmse_test),
  MAE = c(mae_test),
  MAPE = c(mape_test)
)

print(model_evaluation)

# Perform ICE analysis and generate plots
features <- c("NH4", "pH", "DO", "T", "SRT", "HRT", "TOC", "TP", "TN")
for (feature in features) {
  try({
    cat("Processing feature:", feature, "\n")
    ice_obj <- partial(best_model, pred.var = feature, ice = TRUE, n.trees = best_iter, grid.resolution = 100)
    ice_plot <- autoplot(ice_obj, rug = TRUE, train = GBRT_train) + 
      ggtitle(paste("ICE Plot for", feature))
    print(ice_plot)
  }, silent = TRUE)
}

# Generate 3D interactive plot and heatmap based on ICE analysis
features <- c("DO", "SRT")
pdp_result <- partial(best_model, pred.var = features, n.trees = best_iter, grid.resolution = 100)

# Extract unique values
XX <- unique(pdp_result[[features[1]]])
YY <- unique(pdp_result[[features[2]]])
yhat_values <- matrix(pdp_result$yhat, nrow = length(YY), byrow = TRUE)  # Note: arranged by row

# Downsample extracted values
downsample <- function(data, step=10) {
  sapply(seq(1, length(data), by=step), function(i) mean(data[i:min(i+step-1, length(data))]))
}
XX_downsampled <- downsample(XX, step = 10)
YY_downsampled <- downsample(YY, step = 10)

# Create downsampled ZZ matrix
ZZ_downsampled <- matrix(0, nrow = length(YY_downsampled), ncol = length(XX_downsampled))

for (i in 1:length(YY_downsampled)) {
  for (j in 1:length(XX_downsampled)) {
    ZZ_downsampled[i, j] <- mean(yhat_values[
      ((i - 1) * 10 + 1):min(i * 10, length(YY)),
      ((j - 1) * 10 + 1):min(j * 10, length(XX))
    ])
  }
}

# Colorscale definition
colorscale <- list(
  list(0, '#619DB8'),
  list(0.2, '#aecdd7'),
  list(0.4, '#e3eeef'),
  list(0.6, '#fae7d9'),
  list(0.8, '#f0b79a'),
  list(1, '#c85d4d')
)

# Plot 3D surface plot
fig_surface <- plot_ly(
  x = ~XX_downsampled,
  y = ~YY_downsampled,
  z = ~ZZ_downsampled,
  type = 'surface',
  colorscale = colorscale
)

# Adjust 3D surface plot layout and font
fig_surface <- fig_surface %>% layout(
  scene = list(
    xaxis = list(
      title = list(text = 'DO', font = list(family = "Arial", size = 18)),
      tickfont = list(family = "Arial", size = 15)
    ),
    yaxis = list(
      title = list(text = 'SRT', font = list(family = "Arial", size = 18)),
      tickfont = list(family = "Arial", size = 15)
    ),
    zaxis = list(
      title = list(text = 'Partial Dependence', font = list(family = "Arial", size = 18)),
      tickfont = list(family = "Arial", size = 15)
    ),
    aspectratio = list(x = 1, y = 1, z = 0.7)
  ),
  coloraxis_colorbar = list(
    title = list(text = 'Partial Dependence', font = list(family = "Arial", size = 15)),
    tickvals = seq(0, 1, by = 0.25),
    ticktext = c('Low', '', '', '', 'High'),
    x = 1.1,
    lenmode = 'fraction',
    len = 0.1,
    thickness = 2
  ),
  margin = list(t = 0),
  showlegend = FALSE
)

# Plot heatmap
fig_heatmap <- plot_ly(
  x = ~XX,
  y = ~YY,
  z = ~yhat_values,
  type = 'heatmap',
  colorscale = colorscale
)

# Adjust heatmap layout and font
fig_heatmap <- fig_heatmap %>% layout(
  xaxis = list(title = list(text = 'DO', font = list(family = "Arial", size = 18))),
  yaxis = list(title = list(text = 'SRT', font = list(family = "Arial", size = 18))),
  coloraxis_colorbar = list(
    title = list(text = 'Partial Dependence', font = list(family = "Arial", size = 15)),
    tickvals = seq(0, 1, by = 0.25),
    ticktext = seq(0, 1, by = 0.2), 
    lenmode = 'fraction',
    len = 0.1,  # Shorten colorbar length
    thickness = 2  # Narrow colorbar width
  )
)
