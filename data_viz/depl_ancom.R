#####################################################################
####################__ANCOMBC2_depletion_study__#####################
#####################################################################

library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(ape)
library(microViz)
library(DT)

# Load in data from nf-core ampliseq
raw_physeq <- readRDS("/path/to/dada2_phyloseq.rds")

# Now we read in the phylo tree computed by qiime2
tree <- read.tree("/path/to/tree.nwk")

# And add it to our phyloseq object
physeq <- merge_phyloseq(raw_physeq, phy_tree(tree))

View(sample_data(physeq))

# Change this up depending on whatever metadata you want to compare
gooder = subset_samples(physeq = physeq, Timepoint == "BA" & Treatment == "VEH" 
                        | Timepoint == "5dpi2" & Treatment == "VEH")

# Read in phyloseq for ancombc
tse <- mia::makeTreeSummarizedExperimentFromPhyloseq(gooder)

# set seed for ancombc
set.seed(123)

# run ancombc
output3 <- ancombc2(data = tse, assay_name = "counts", tax_level = "Genus",
                   fix_formula = "Timepoint", rand_formula = NULL,
                   p_adj_method = "holm", pseudo_sens = TRUE,
                   prv_cut = 0.1, lib_cut = 10000, s0_perc = 0.05,
                   group = "Timepoint", struc_zero = TRUE, neg_lb = TRUE,
                   alpha = 0.05, n_cl = 2, verbose = TRUE,
                   global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                   iter_control = list(tol = 1e-2, max_iter = 20, 
                                       verbose = TRUE),
                   em_control = list(tol = 1e-5, max_iter = 100),
                   lme_control = lme4::lmerControl(),
                   mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                   trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(-1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE),
                                                        matrix(c(1, 0, 1, -1),
                                                               nrow = 2, 
                                                               byrow = TRUE)),
                                        node = list(2, 2, 1),
                                        solver = "ECOS",
                                        B = 100))
# How many structural zeroes?
tab_zero <- output3$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")

#Print output to dataframe
res_prim3 <- output3$res

#Save ancombc stats for each subset
write.csv(res_prim3, file = "depl_ancom_timepoint_BA25dpi2_onlyVEH.csv", row.names = FALSE)

#I then used GraphPad prism to make pretty plots using this data
