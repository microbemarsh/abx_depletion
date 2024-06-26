#####################################################################################
########################__ORDINATION_&_RELATIVE_ABUNDANCE__##########################
#####################################################################################

library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(ape)
library(ggsci)
library(microViz)

# Load data from nf-core ampliseq
raw_physeq <- readRDS("/path/to/dada2_phyloseq.rds")

# Now we read in the phylo tree computed by qiime2
tree <- read.tree("/path/to/tree.nwk")

# And add it to our phyloseq object
physeq <- merge_phyloseq(raw_physeq, phy_tree(tree))

#####################################################################################
############################____ALPHA__DIVERSITY____#################################

# Define the character strings and columns for filtering
filter_column1 <- "Timepoint"
filter_string1 <- "7dpi1"

filter_column2 <- "Treatment"
filter_string2 <- "VEH"

# Filter the phyloseq object
filtered_physeq <- subset_samples(depl, 
                                    get(filter_column1) == filter_string1 & 
                                      get(filter_column2) == filter_string2)

# Print the filtered sample data to check you did it right
print(sample_data(filtered_physeq))

#Prune it
unfilt_data <- prune_species(speciesSums(filtered_physeq) > 0, filtered_physeq)

# Calculate richness measures
richness_data <- estimate_richness(unfilt_data)

# Now make the plots using this data in GraphPad Prism
write.csv(richness_data, file = "depl_7dpi_VEH_alpha_diversity.csv", row.names = FALSE)

#####################################################################################
#########################_____BETA__DIVERSITY______##################################

# Remove OTUs that appear less than 5 times in all samples 
da = physeq
wh0 = genefilter_sample(da, filterfun_sample(function(x) x > 5), A=0.5*nsamples(da))
da1 = prune_taxa(wh0, da)

# Transform to even sampling depth
da1 = transform_sample_counts(da1, function(x) 1E6 * x/sum(x))

# Keep 5 most abundant taxa in all samples
phylum.sum = tapply(taxa_sums(da1), tax_table(da1)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
da1 = prune_taxa((tax_table(da1)[, "Phylum"] %in% top5phyla), da1)

##############################################################################################
####################_________ORDINATION_TIMEPOINTS_&_TREATMENTS__________#####################

# This section was repeated for each timepoint with the independent variable being the ABX

# Subsetting data
da1_baseline <- subset_samples(da1, Timepoint == "5dpi2")

#Compute weighted unifrac distances
wunifrac_dist <- distance(da1_baseline, method = "wunifrac")

# Plot PCoA
pcoa_res <- ordinate(da1_baseline, method = "PCoA", distance = "wunifrac")
p1 <- plot_ordination(da1_baseline, pcoa_res, color = "Treatment") +
  geom_point(size = 3) +
  ggtitle("Gut microbiome of 2X TBI 5dpi") +
  theme(panel.background = element_rect(colour = "black", size=1)) +
  stat_ellipse(type = "t", linetype = 2)

p1

# Perform dist matrix
distance_matrix <- phyloseq::distance(da1_baseline, method = "wunifrac")

# Perform PERMANOVA
grouping_variable <- "Treatment"
permanova_result1 <- adonis2(distance_matrix ~ get(grouping_variable), as(sample_data(da1_baseline), "data.frame"))

# Save PERMANOVA results
write.csv(permanova_result1, file = "5dpi_permanova.csv", row.names = FALSE)

######################################################################################
#############################___CUTSOM__COLOR__PALETTES__#############################

good_data <- subset_taxa(physeq, Kingdom == "Bacteria")

clean_data <- phyloseq_validate(good_data) %>%
  tax_fix() %>%
  tax_agg(rank = "Genus")

genus_palette <- tax_palette(
  data = clean_data, rank = "Genus", n = 25, pal = "greenArmytage",
  add = c(Other = "purple")
)

phylum_palette <- tax_palette(
  data = clean_data, rank = "Phylum", n = 6, pal = "greenArmytage",
  add = c(Other = "purple")
)

##############################______PHYLUM______###################################

ba_veh_phylum <- clean_data %>%
  ps_filter(Timepoint == "7dpi1", Treatment == "ABX", .keep_all_taxa = TRUE)

ba_veh_phylum %>%
  comp_barplot("Phylum", n_taxa = 5, merge_other = TRUE, label = NULL, palette = phylum_palette) +
  facet_wrap(vars(Treatment), scales = "free") + # scales = "free" is IMPORTANT!
  ggtitle(
    "7dpi ABX Group",
  ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))

##################################_____GENUS_____##################################

ba_veh_genus <- clean_data %>%
  ps_filter(Timepoint == "BA", Treatment == "ABX", .keep_all_taxa = TRUE)

ba_veh_genus %>%
  comp_barplot("Genus", n_taxa = 10, merge_other = TRUE, label = NULL, palette = genus_palette) +
  facet_wrap(vars(Treatment), scales = "free") + # scales = "free" is IMPORTANT!
  ggtitle(
    "BA ABX Group",
  ) +
  theme(axis.ticks.y = element_blank(), strip.text = element_text(face = "bold"))
