library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(ggExtra)
library(data.table)

# load data
species_id = fread('Species_id2.csv', header = FALSE, sep = ",")
sample_id = fread('Sample_id2.csv', header = FALSE, sep = ",")
Ptrain = fread('OTU.csv', header = FALSE, sep = ",")
ptst = fread('Ztest.csv', header = FALSE, sep = ",")
qtst = fread('qtst.csv', header = FALSE, sep = ",")
qtrn = fread('qtrn.csv', header = FALSE, sep = ",")

print(dim(ptst))
print(dim(qtst))
print(dim(sample_id))
print(dim(species_id))

keystone_predicted = c()
keystone_true = c()

for (i in 1:nrow(qtst)) {
  # predicted null composition
  q_i = qtrn[sample_id$V1[i], ]
  q_i_null = q_i
  q_i_null[species_id$V1[i]] = 0
  q_i_null = q_i_null / sum(q_i_null)
  
  # true null composition
  p_i = Ptrain[, sample_id$V1[i]]
  p_i_null = p_i
  p_i_null[species_id$V1[i]] = 0
  p_i_null = p_i_null / sum(p_i_null)
  
  # The bray-curtis distance was calculated
  BC_true = sum(abs(p_i_null - ptst[, i])) / sum(abs(p_i_null + ptst[, i]))
  BC_pred = sum(abs(q_i_null - qtst[i, ])) / sum(abs(q_i_null + qtst[i, ]))
  
  # Calculate keystoneness
  keystone_predicted = c(keystone_predicted, BC_pred * as.numeric(1 - Ptrain[species_id$V1[i], sample_id$V1[i]]))
  keystone_true = c(keystone_true, BC_true * as.numeric(1 - Ptrain[species_id$V1[i], sample_id$V1[i]]))
}

keystoness = data.frame(str_pred = keystone_predicted, str_true = keystone_true)

write.table(keystoness, file = "keystoness.csv", sep = ",", row.names = FALSE, col.names = TRUE)
