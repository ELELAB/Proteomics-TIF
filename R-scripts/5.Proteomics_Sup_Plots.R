# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PROTEOMICS PLOTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Load packages
library(viridis)
library(heatmap.plus)
library(ggplot2)
library(biomaRt)
library(sva)
library(openxlsx)
library(reshape2)
library(topGO)
library(GOSim)
library(biomaRt)
library(gprofiler2)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR PLOTTING HEATMAPS
# Take arguments:
# my.df = matrix with expression/abundance values
# my.prots = vector of differentially abundant genes/proteins
# my.groups = vector of group IDs, must be as.factor()
# my.cvar = vector with covariate for batch correction
# my.col = vector with one color for each group level
# my.s and my.c = scale and center values for plotting TRUE/FALSE
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

HeatmapFunction <- function(my.df, my.prots, my.group, my.cvar, my.col, my.s, my.c) {
  
  my.df <- combat_corrections(my.df, my.group, my.cvar)
  my.df <- my.df[rownames(my.df) %in% my.prots,]
  print(dim(my.df))
  heat.cols <- viridis(n=100, option = "magma")
  my.col <- get_colors(my.group, my.col)
  my.spacer <- as.matrix(replicate(nrow(my.col), "white"))
  my.cols <- cbind(my.spacer, my.col)
  
  heatmap.plus(as.matrix(scale(my.df, scale = my.s, center = my.c)), col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labRow="", labCol="", ColSideColors=my.cols, margins = c(14,8), cexCol=1.2, cexRow = 1.3)
}




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load sets of differentially abundant proteins for plotting.

my.wd <- ""
setwd(my.wd)

H_L_DA <- read.delim("Subtypes/Corrected_for_Pool_permissive/HER2_LumA_DA_corrected_pool.txt", header = TRUE)
H_T_DA <- read.delim("Subtypes/Corrected_for_Pool_permissive/HER2_TNBC_DA_corrected_pool.txt", header = TRUE)
L_T_DA <- read.delim("Subtypes/Corrected_for_Pool_permissive/LumA_TNBC_DA_corrected_pool.txt", header = TRUE)

STP_DA <- rbind(H_L_DA, H_T_DA, L_T_DA)

ER_DA <- read.delim("covariates/PGR_ER_HER2/ERp_ERn_DA_corrected_pool.txt", header = TRUE)
PgR_DA <- read.delim("covariates/PGR_ER_HER2/PGRp_PGRn_DA_corrected_pool.txt", header = TRUE)
Her2_DA <- read.delim("covariates/PGR_ER_HER2/H23_H01_corrected_pool.txt", header = TRUE)

TILs_DA <- read.delim("covariates/TILs_and_grade/TILs_High_Low_DA_corrected_pool.txt", header = TRUE)
GR_DA <- read.delim("covariates/TILs_and_grade/GR_High_Low.txt", header = TRUE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf("HEAMP_subtype.pdf")
HeatmapFunction(tifdata_small, unique(STP_DA$Accession), diag, pool, c("#FAA275","#DDDC71", "#BF5454"), FALSE, TRUE)
dev.off()

pdf("HEAMP_ER.pdf")
HeatmapFunction(tifdata_small, ER_DA$Accession, ER, pool, c("#FFFCF7", "grey60"), FALSE, TRUE)
dev.off()

pdf("HEAMP_PgR.pdf")
HeatmapFunction(tifdata_small, PgR_DA$Accession, PGR, pool, c("#FFFCF7", "grey60"), FALSE, FALSE)
dev.off()

pdf("HEAMP_Her2.pdf")
HeatmapFunction(tifdata_small, Her2_DA$Accession, HER2, pool, c("#FFFCF7", brewer.pal(4,"Greys")[2:4]), FALSE, FALSE)
dev.off()

pdf("HEAMP_TILs.pdf")
HeatmapFunction(tifdata_small, TILs_DA$Accession, TILs, pool, c("#032A63","#0071AA"), FALSE, FALSE)
dev.off()

pdf("HEAMP_GR.pdf")
HeatmapFunction(tifdata_small, GR_DA$Accession, Gr, pool, c("#032A63","#0071AA"), TRUE, FALSE)
dev.off()



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH CORRECTION - MULTIDIMENSIONAL SCALING PLOTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

pdf("BeforeBC.pdf")
myMDSplot(tifdata_small, diag, "", colorcode_diag)
dev.off()

pdf("AfterBC.pdf")
myMDSplot(batch_corr_diag, diag, "", colorcode_diag)
dev.off()



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH CORRECTION - BOX PLOTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

s <- sample(seq(nrow(tifdata_small)), 20)

IDconvert <- uniprot_to_name(rownames(tifdata_small[s,]))[1:12,]
prots <- IDconvert$Accession
name <-IDconvert$name


Before <- data.frame(t(tifdata_small[rownames(tifdata_small) %in% prots,]))
colnames(Before) <- name
Before$diag <- diag
Before$group <- rep("Before Batch Correction", nrow(Before))
Before <- melt(Before)

After <- data.frame(t(batch_corr_diag[rownames(batch_corr_diag) %in% prots,]))
colnames(After) <- name
After$diag <- diag
After$group <- rep("After Batch Correction", nrow(After))
After <- melt(After)

dat <- rbind(Before, After)
dat$group <- factor(dat$group, levels = c("Before Batch Correction", "After Batch Correction"))
dat$variable <- as.factor(dat$variable)
p2 <- ggplot(dat, aes(x=diag, y=value, fill=group)) + scale_fill_manual(values = c("#F7A072", "#645E9D")) +
  theme_bw() +
  geom_boxplot() + theme(panel.grid.major.y = element_blank(), strip.background = element_rect(color="black", fill="white")) +
  facet_wrap(~variable)

pdf("BatchCorrectionBoxplot.pdf", height = 6, width = 10)
p2
dev.off()




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GENE ONTOLOGY ENRICHMENT PLOTS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION EXTRACTING GO-TERMS
# Take arguments:
# univ.genes  = background set of genes/proteins
# DA.set = sets if differentially expressed/abundant genes/proteins
# BG.BP and BG.MF = backgrounds sets "Biological Process" and "Molecular Function".
# n = number top terms to return
# comparison = groups comparison denoted by a string
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


GO_BP_MF <- function(univ.genes, DA.set, BG.BP, BG.MF, n, comparison){
  BP <- TOPGO("BP", univ.genes, DA.set$name, BG.BP, n)
  MF <- TOPGO("MF", univ.genes, DA.set$name, BG.MF, n)
  res <- rbind(BP, MF)
  res$Comparison <- rep(comparison, nrow(res))
  return(res)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Univers
univ_genes <- sort(unique(tif_full$name))[-1]

# Background molecular function and biological process
Gene_background_BP <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "symbol")
Gene_background_MF <- annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "symbol")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get Terms

H_L_Gos <- GO_BP_MF(univ_genes, H_L_DA, Gene_background_BP, Gene_background_MF, 20, "Her2 vs Luminal")
H_T_Gos <- GO_BP_MF(univ_genes, H_T_DA, Gene_background_BP, Gene_background_MF, 20, "Her2 vs TNBC")
L_T_Gos <- GO_BP_MF(univ_genes, L_T_DA, Gene_background_BP, Gene_background_MF, 20, "Lumnial vs TNBC")

ER_Gos <- GO_BP_MF(univ_genes, ER_DA, Gene_background_BP, Gene_background_MF, 20, "ER+ vs ER-")
#PgR_Gos <- GO_BP_MF(univ_genes, PgR_DA, Gene_background_BP, Gene_background_MF, 20, "PgR+ vs PgR-")
#Her2_Gos <- GO_BP_MF(univ_genes, Her2_DA, Gene_background_BP, Gene_background_MF, 20, "High Her2 vs Low Her2")

TILs_Gos <- GO_BP_MF(univ_genes, TILs_DA, Gene_background_BP, Gene_background_MF, 20, "High TILs vs Low TILs")
#GR_Gos <- GO_BP_MF(univ_genes, GR_DA, Gene_background_BP, Gene_background_MF, 20, "High Grade vs Low Grade")

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# bind terms together for plotting

Godat <- rbind(H_L_Gos, H_T_Gos, L_T_Gos, ER_Gos, TILs_Gos)
Godat$Fraction <- (Godat$Significant/Godat$Annotated)+1
#Godat$Type <- as.factor(Godat$Type)
Godat <- Godat[order(Godat$Fraction, decreasing = FALSE),]
Godat$Term <- factor(Godat$Term, levels = unique(as.character(Godat$Term)))
Godat$Comparison <- factor(Godat$Comparison, levels = c("Her2 vs Luminal", "Her2 vs TNBC", "Lumnial vs TNBC", "ER+ vs ER-", "High TILs vs Low TILs"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# plot terms
pdf("GOplot.pdf", height = 10, width = 14)
ggplot(Godat, aes(x=Type, 
           y=Term, 
           colour=FDR_fisher, 
           size=Fraction)) + theme_bw() + scale_color_viridis() +
  geom_point() + facet_wrap(~Comparison, scales = "free", ncol = 3) + theme(panel.grid.minor.x = element_blank(), strip.background = element_rect(color="black", fill="white")) + xlab("") + ylab("")
dev.off()

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# write results to tables
write.table(H_L_Gos, "H_L_Gos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(H_T_Gos, "H_T_Gos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(L_T_Gos, "L_T_Gos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(ER_Gos, "ER_Gos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(TILs_Gos, "TILs_Gos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Pathways Enrichment Analysis
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Use the gprofiler2 package to get pathways


# Her2 vs Luminal
H_L_pw <- gost(query = H_L_DA$name, 
               organism = "hsapiens", ordered_query = FALSE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "custom", custom_bg = univ_genes, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE)


H_L_pw <- H_L_pw$result[grep("KEGG|REAC|WP", H_L_pw$result$term_id),c(9,11,4,6,3,10)]


# Her2 vs TNBC
H_T_pw <- gost(query = H_T_DA$name, 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", 
                domain_scope = "custom", custom_bg = univ_genes, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)

H_T_pw <- H_T_pw$result[grep("KEGG|REAC|WP", H_T_pw$result$term_id),c(9,11,4,6,3,10)]


# Luminal vs TNBC
L_T_pw <- gost(query = L_T_DA$name, 
               organism = "hsapiens", ordered_query = FALSE, 
               multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
               measure_underrepresentation = FALSE, evcodes = FALSE, 
               user_threshold = 0.05, correction_method = "fdr", 
               domain_scope = "custom", custom_bg = univ_genes, 
               numeric_ns = "", sources = NULL, as_short_link = FALSE)

L_T_pw <- L_T_pw$result[grep("KEGG|REAC|WP", L_T_pw$result$term_id),c(9,11,4,6,3,10)]


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# write results to tables
write.table(H_L_pw, "H_L_pw.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(H_T_pw, "H_T_pw.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(L_T_pw, "L_T_pw.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
