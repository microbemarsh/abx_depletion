# microbe-metabolome

library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(ggsci)
library(microViz)
library(reshape2)
library(tidyr)
library(pheatmap)
library(Hmisc)
library(qvalue)

# Load data from nf-core ampliseq
raw_physeq <- readRDS("dada2_phyloseq.rds")

# Now we read in the phylo tree computed by qiime2
tree <- read.tree("tree.nwk")

# And add it to our phyloseq object
physeq <- merge_phyloseq(raw_physeq, phy_tree(tree))

good_data <- subset_taxa(physeq, Kingdom == "Bacteria")

clean_data <- phyloseq_validate(good_data) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")

View(otu_table(clean_data))

# This is super ugly but worked at one point in time :/

write.csv(clean_data@otu_table, file = "microbe-metabo.csv", row.names = TRUE)

microb_metabo = read.csv("microbe-metabo.csv", row.names = 1) %>%
  select(matches("^Villapol_5D"))

write.csv(microb_metabo, file = "microbiome4metbo.csv", row.names = TRUE)

# Load your microbial counts and metabolomic data
# Replace 'microbial_counts.csv' and 'metabolomic_data.csv' with your actual file paths
microbial_counts <- read.csv("microbiome4metbo.csv", row.names = 1)  # Assuming sample IDs are row names
metabolomic_data <- read.csv("scfa_notreat.csv", row.names = 1, check.names = FALSE)  # Assuming sample IDs are row names
metabolomic_data <- metabolomic_data[, sapply(metabolomic_data, is.numeric)]

############## heatmap #############

# Load the OTU table (replace 'otu_table.csv' with your file name)
otu_table <- read.csv("good_micrbbe4scfa.csv", row.names=1, check.names=FALSE)

# Calculate the total counts in each sample
sample_sums <- colSums(otu_table)

# Normalize the OTU table by dividing each count by the total count in its respective sample
normalized_otu_table <- t(t(otu_table) / sample_sums)

# Optionally, you may want to transform the data (e.g., log transformation) for further analysis
normalized_otu_table <- log2(normalized_otu_table + 1)  # Log2 transformation

# Write the normalized OTU table to a new file (replace 'normalized_otu_table.csv' with your desired output file name)
write.csv(normalized_otu_table, file="good_norm_microbe4_metabo.csv")


########################## USE THIS FOR SCFA / TAXA CORRELATION #############################


# read in good .csv
data = read.csv("mynorm_ABX_scfa_n_taxa.csv", header = TRUE, check.names = FALSE)

data_long <- data %>%
  pivot_longer(cols = c(Acetate:Caproate), names_to = "SCFA", values_to = "SCFA_value") %>%
  inner_join(pivot_longer(data, cols = c(Bacteroides:Lactobacillus), names_to = "Taxa", values_to = "Taxa_value"), by = "Sample", relationship = "many-to-many")

# Calculate correlation between SCFAs and taxa
correlation <- data_long %>%
  group_by(SCFA, Taxa) %>%
  summarise(correlation = cor(SCFA_value, Taxa_value)) %>%
  pivot_wider(names_from = Taxa, values_from = correlation)

# Create heatmap
heatmap_data <- correlation %>%
  select(-SCFA) %>%
  tibble::column_to_rownames(var = "SCFA")

# Plot heatmap
pheatmap(heatmap_data, 
         main = "Correlation between SCFAs and Taxa", 
         fontsize = 8,
         show_colnames = TRUE,
         show_rownames = TRUE)

# Compute correlation between SCFAs and taxa

data_no_name = data |>
  tibble::column_to_rownames("Sample")

correlation <- cor(data_no_name, use = "pairwise.complete.obs")

# Compute p-values for each correlation coefficient
p_values <- sapply(2:ncol(data_no_name), function(i) {
  cor.test(data_no_name[, 1], data_no_name[, i], method = "pearson", use = "pairwise.complete.obs")$p.value
})

# Replace non-significant p-values with NA
p_values[!p.adjust(p_values, method = "BH") < 0.05] <- NA

# Print correlation matrix with p-values
print(correlation)
print(p_values)

write.csv(p_values, file = "mynorm_abx_htmp_pvals.csv")

# Write correlation matrix with p-values to CSV file
correlation_with_pvalues <- correlation
correlation_with_pvalues$p_values <- p_values
write.csv(correlation_with_pvalues, "mynorm_abx_microbes_corr_matrix_with_pvalues.csv")

######################## Trying it ######################

# Remove the 'Sample' column for correlation analysis
data <- data[, -1]

# Calculate correlation matrix
cor_matrix <- cor(data, use = "complete.obs", method = "pearson")

# Calculate p-values matrix
p_matrix <- matrix(nrow = ncol(data), ncol = ncol(data), dimnames = list(colnames(data), colnames(data)))
for (i in 1:ncol(data)) {
  for (j in i:ncol(data)) {
    test <- cor.test(data[[i]], data[[j]], method = "pearson")
    p_matrix[i, j] <- test$p.value
    p_matrix[j, i] <- test$p.value  # since the matrix is symmetric
  }
}

# Adjust p-values using the Benjamini-Hochberg method
adj_p_matrix <- apply(p_matrix, 1, function(p) p.adjust(p, method = "BH"))

# Save results to CSV files
write.csv(cor_matrix, "nynorm_abx_correlation_coefficients.csv", row.names = TRUE)
write.csv(adj_p_matrix, "nynorm_abx_adjusted_p_values.csv", row.names = TRUE)

# Calculate correlations and p-values
cor_results <- rcorr(as.matrix(data[, -1]), type = "pearson")  # Assuming you drop the 'Sample' column
cor_matrix <- cor_results$r  # Correlation coefficients
p_matrix <- cor_results$P  # p-values

# Adjust the p-values for multiple testing
adj_p_matrix <- p.adjust(p_matrix, method = "BH")

# Set a significance level
alpha <- 0.05

# Create a matrix to identify significant correlations
significant <- adj_p_matrix < alpha

# Combine the data to display significant correlations and their coefficients only if significant
significant_correlations <- cor_matrix
significant_correlations[!significant] <- NA  # Replace non-significant correlations with NA

# Save significant correlations to CSV
write.csv(significant_correlations, "mynorm_abx_significant_correlations.csv", row.names = TRUE)

# Optionally, print significant correlations
print(significant_correlations)

########################################

# Load the data
data_no_name <- data[, -1]

# Extract SCFA and taxa columns
scfa_columns <- data_no_name[, 1:10]
taxa_columns <- data_no_name[, 11:ncol(data_no_name)]

# Compute correlation matrix
correlation <- cor(scfa_columns, taxa_columns, use = "pairwise.complete.obs")

# Compute p-values for each correlation coefficient
p_values <- sapply(1:ncol(scfa_columns), function(i) {
  cor.test(scfa_columns[, i], taxa_columns[, i], method = "pearson", use = "pairwise.complete.obs")$p.value
})

# Adjust p-values for multiple testing using Holm's correction
adjusted_p_values <- p.adjust(p_values, method = "holm")

# Identify significant correlations (p < 0.05 after correction)
significant_correlations <- which(adjusted_p_values < 0.05)

# Print significant correlations
print("\nSignificant Correlations:")
if (length(significant_correlations) > 0) {
  for (i in significant_correlations) {
    scfa_name <- colnames(scfa_columns)[i]
    taxa_name <- colnames(taxa_columns)[i]
    cat("SCFA:", scfa_name, "and Taxa:", taxa_name, "have a significant correlation (p =", p_values[i], ", adjusted p =", adjusted_p_values[i], ")\n")
  }
} else {
  print("No significant correlations found.")
}


###############################

# Load the data
data_no_name <- data[, -1]

# Compute correlation matrix
correlation <- cor(data_no_name, use = "pairwise.complete.obs")

# Compute p-values for each correlation coefficient
p_values <- sapply(2:ncol(data_no_name), function(i) {
  cor.test(data_no_name[, 1], data_no_name[, i], method = "pearson", use = "pairwise.complete.obs")$p.value
})

# Adjust p-values for multiple testing using Holm's correction
adjusted_p_values <- p.adjust(p_values, method = "holm")

# Identify significant correlations (p < 0.05 after correction)
significant_correlations <- which(adjusted_p_values < 0.05)

# Print significant correlations
print("\nSignificant Correlations:")
if (length(significant_correlations) > 0) {
  for (i in significant_correlations) {
    col_name <- colnames(data_no_name)[i + 1]  # Skip the first column (Sample)
    scfa_name <- colnames(data_no_name)[1]
    cat(scfa_name, "and", col_name, "have a significant correlation (p =", p_values[i], ", adjusted p =", adjusted_p_values[i], ")\n")
  }
} else {
  print("No significant correlations found.")
}

##################################################

# read in good .csv
data_w_name = read.csv("mynorm_ABX_scfa_n_taxa.csv", header = TRUE, check.names = FALSE)

data = data_w_name[, -1]

set.seed(123)
num_samples <- 18
num_variables <- 20
data <- matrix(rnorm(num_samples * num_variables), nrow = num_samples)

# Perform statistical tests (e.g., correlation)
# Replace this with your actual statistical test
p_values <- matrix(NA, nrow = num_variables, ncol = num_variables)
for (i in 1:num_variables) {
  for (j in 1:num_variables) {
    if (i != j) {
      # Example: Pearson correlation test
      p_values[i, j] <- cor.test(data[, i], data[, j])$p.value
    }
  }
}

# Adjust p-values for multiple comparisons
# Example: Benjamini-Hochberg procedure
adjusted_p_values <- p.adjust(as.vector(p_values), method = "BH")
adjusted_p_values <- matrix(adjusted_p_values, nrow = num_variables, ncol = num_variables)


# Write results to CSV
write.csv(result_df, "mynorm_abx_adjusted_p_values_results.csv", row.names = FALSE)
