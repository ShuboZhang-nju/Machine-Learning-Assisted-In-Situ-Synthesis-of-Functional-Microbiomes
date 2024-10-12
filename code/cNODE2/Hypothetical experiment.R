library(data.table)
library(parallel)
library(doParallel)
library(pbapply)


# load data
data <- fread('OTU.txt', header = TRUE)
setDF(data)
rownames(data) <- data[, 1]
data <- data[, -1]
N <- nrow(data) 
M <- ncol(data) 
steady_state_absolute <- as.matrix(data)

absent_collection <- matrix(nrow = N, ncol = 0)
species_id <- c()
sample_id <- c()

cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

clusterExport(cl, varlist = c("N", "M", "steady_state_absolute"))

results <- pblapply(1:N, function(j1) {
  local_absent_collection <- matrix(nrow = N, ncol = 0)
  local_species_id <- c()
  local_sample_id <- c()
  
  for (j2 in 1:M) {
    if (steady_state_absolute[j1, j2] > 0) {
      y_0 <- steady_state_absolute[, j2]
      y_0_binary <- y_0 > 0
      if (sum(y_0_binary) > 1) {
        y_0[j1] = 0  
        local_species_id <- c(local_species_id, j1)
        local_sample_id <- c(local_sample_id, j2)
        y_0[y_0 > 0] = 1
        local_absent_collection <- cbind(local_absent_collection, y_0 / sum(y_0))
      }
    }
  }
  print(paste("Completed task for species:", j1))  
  list(absent_collection = local_absent_collection, species_id = local_species_id, sample_id = local_sample_id)
}, cl = cl)


stopCluster(cl)

final_absent_collection <- matrix(nrow = N, ncol = 0)
final_species_id <- c()
final_sample_id <- c()

for (result in results) {
  final_absent_collection <- cbind(final_absent_collection, result$absent_collection)
  final_species_id <- c(final_species_id, result$species_id)
  final_sample_id <- c(final_sample_id, result$sample_id)
}

print(paste("Final absent collection dimensions:", dim(final_absent_collection)))
print(paste("Final species ID length:", length(final_species_id)))
print(paste("Final sample ID length:", length(final_sample_id)))

fwrite(as.data.table(final_absent_collection), file = 'Ztest.csv', row.names = FALSE, col.names = FALSE, sep = ",")
fwrite(as.data.table(final_species_id), file = 'Species_id.csv', row.names = FALSE, col.names = FALSE, sep = ",")
fwrite(as.data.table(final_sample_id), file = 'Sample_id.csv', row.names = FALSE, col.names = FALSE, sep = ",")