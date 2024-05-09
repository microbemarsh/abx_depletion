# SCFA with MetaboAnalystR

library(MetaboAnalystR)
library(dplyr)
library(ellipse)
library(iheatmapr)
library(pls)

# Load in data
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "bcm_istdnorm_scfaonly.csv", "rowu", "disc")

# Check to make sure it was loaded properly
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)

# Missing values
mSet <- RemoveMissingPercent(mSet, percent=0.5)
mSet <- ImputeMissingVar(mSet, method="exclude")
mSet<-IsSmallSmplSize(mSet)

#Normalization to ISTD (can also perform others)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ref=NULL, ratio=FALSE, ratioNum=20)

# Normalization Viz
mSet<-PlotNormSummary(mSet, "feature_norm", format="png", dpi=72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "sample_norm", format="pdf", width=NA)

# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
mSet$analSet$fc$fc.log

# Perform PCA analysis
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "pdf", dpi = 72, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "pdf", dpi = 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "pdf", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "pdf", 72, width=NA, 1,2,3, 40)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "pdf", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "pdf", dpi = 72, width=NA, 1, 2)

# PLS Analysis
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "pdf", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,1,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "pdf", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "pdf", 72, width=NA, 1, 2);
mSet<-PLSDA.CV(mSet, "5", 5,5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "pdf", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "pdf", 72, width=NA, "vip", "Comp. 1", 15, FALSE)
mSet<-PLSDA.Permut(mSet, 1000, "accu")
mSet<-PlotPLS.Permutation(mSet, "pls_perm_1_", "pdf", 72, width=NA)

# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "pdf", dpi=72, width=NA, 1,2,0.95,1,0)
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "pdf", 72, width=NA);
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "pdf", 72, width=NA, "vip", "tscore", 15,FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "pdf", dpi=72, width=NA)
mSet<-OPLSDA.Permut(mSet, 1000)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "pdf", dpi=72, width=NA)

######################## Now my normalization #############################################

# Load in data
mSet <- InitDataObjects("pktable", "stat", FALSE)
mSet <- Read.TextData(mSet, "bcm_istdnorm_scfaonly.csv", "rowu", "disc")

# Check to make sure it was loaded properly
mSet <- SanityCheckData(mSet)
mSet <- ReplaceMin(mSet)

# Missing values
mSet <- RemoveMissingPercent(mSet, percent=0.5)
mSet <- ImputeMissingVar(mSet, method="exclude")
mSet<-IsSmallSmplSize(mSet)

#Normalization to ISTD (can also perform others)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ref=NULL, ratio=FALSE, ratioNum=20)

# Normalization Viz
mSet<-PlotNormSummary(mSet, "feature_norm", format="png", dpi=72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "sample_norm", format="pdf", width=NA)

# Perform fold-change analysis on uploaded data, unpaired
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
mSet$analSet$fc$fc.log

# Perform PCA analysis
mSet<-PCA.Anal(mSet)
mSet<-PlotPCAPairSummary(mSet, "pca_pair_0_", format = "pdf", dpi = 72, width=NA, 5)
mSet<-PlotPCAScree(mSet, "pca_scree_0_", "pdf", dpi = 72, width=NA, 5)
mSet<-PlotPCA2DScore(mSet, "pca_score2d_0_", format = "pdf", dpi=72, width=NA, 1, 2, 0.95, 1, 0)
mSet<-PlotPCA3DScoreImg(mSet, "pca_score3d_0_", "pdf", 72, width=NA, 1,2,3, 40)
mSet<-PlotPCALoading(mSet, "pca_loading_0_", "pdf", 72, width=NA, 1,2);
mSet<-PlotPCABiplot(mSet, "pca_biplot_0_", format = "pdf", dpi = 72, width=NA, 1, 2)

# PLS Analysis
mSet<-PLSR.Anal(mSet, reg=TRUE)
mSet<-PlotPLSPairSummary(mSet, "pls_pair_0_", "pdf", 72, width=NA, 5)
mSet<-PlotPLS2DScore(mSet, "pls_score2d_0_", "pdf", 72, width=NA, 1,2,0.95,1,0)
mSet<-PlotPLS3DScoreImg(mSet, "pls_score3d_0_", "pdf", 72, width=NA, 1,2,3, 40)
mSet<-PlotPLSLoading(mSet, "pls_loading_0_", "pdf", 72, width=NA, 1, 2);
mSet<-PLSDA.CV(mSet, "5", 5,5, "Q2")
mSet<-PlotPLS.Classification(mSet, "pls_cv_0_", "pdf", 72, width=NA)
mSet<-PlotPLS.Imp(mSet, "pls_imp_0_", "pdf", 72, width=NA, "vip", "Comp. 1", 15, FALSE)
mSet<-PLSDA.Permut(mSet, 1000, "accu")
mSet<-PlotPLS.Permutation(mSet, "pls_perm_1_", "pdf", 72, width=NA)

# Perform oPLS-DA analysis
mSet<-OPLSR.Anal(mSet, reg=TRUE)
mSet<-PlotOPLS2DScore(mSet, "opls_score2d_0_", format = "pdf", dpi=72, width=NA, 1,2,0.95,1,0)
mSet<-PlotOPLS.Splot(mSet, "opls_splot_0_", "all", "pdf", 72, width=NA);
mSet<-PlotOPLS.Imp(mSet, "opls_imp_0_", "pdf", 72, width=NA, "vip", "tscore", 15,FALSE)
mSet<-PlotOPLS.MDL(mSet, "opls_mdl_0_", format = "pdf", dpi=72, width=NA)
mSet<-OPLSDA.Permut(mSet, 1000)
mSet<-PlotOPLS.Permutation(mSet, "opls_perm_2_", format = "pdf", dpi=72, width=NA)

# Heatmap
mSet<-PlotHeatMap2(mSet, "heatmap2_0_", "norm", "row", "png", 72, width=NA, 
                   "euclidean","ward.D","bwm","overview","mean",2000, F, F, F, F)

# Perform hierarchical clustering and plot dendogram
mSet<-PlotHCTree(mSet, "tree_0_", "png", dpi=72, width=NA, "euclidean", "ward.D")

mSet<-PlotCorrHeatMap(mSet, "corr_0_", "png", 72, width=NA, fix.col = TRUE, 
                      "col", "pearson", "bwm", "overview", F, F, 0.0)

mSet<-PlotHeatMap(mSet, "heatmap_0_", "png", 72, width=NA, 
                  "norm", "row", "euclidean", "ward.D","bwm", 8, "overview", T, T, NULL, T, F, T, T, T)

normalized_set <- mSet[["dataSet"]][["norm"]]
ordered_normalized_set <- normalized_set[order(row.names(normalized_set)), ]
write.csv(normalized_set, file = "depletion_scfa_norm_mSet.csv")
