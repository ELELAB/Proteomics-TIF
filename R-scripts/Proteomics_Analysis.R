#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# ------------------------------------------------------------------------------------
#                         ANALYSIS OF PROTEOMICS DATA
# ------------------------------------------------------------------------------------




# ------------------------------------------------------------------------------------
# GET UNIPROT ID - GENE ID
# ------------------------------------------------------------------------------------
tif_full_id <- as.character(sort(rownames(tifdata_small)))
tif_full_name <- unique(uniprot_to_name(tif_full_id)$name)
tif_full  <- unique(uniprot_to_name(tif_full_id))




#---------------------------------------------------------------------------
# GET OVERLAP BETWEEN PLASMA; SECRETED N-GLYCOSYLATED AND TIF.
# -------------------------------------------------------------------------------------

nglyco_name <- unique(uniprot_to_name(nglyco)$name)
secreted_name <- unique(uniprot_to_name(secreted)$name)
plasma_name <- unique(uniprot_to_name(plasma)$name)
exosomes_name <- sort(as.character(exosomes))


intsec.list <- list(tif_full_name, nglyco_name, plasma_name, secreted_name, exosomes_name)
names(intsec.list) <- c("TIF", "Nglyco", "Plasma", "Secreted", "Exosomes")


colov <- c("#08605F", "#177E89", "#598381","#8E936D", "#5D737E")
highly_expressed <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_TIF_Sec_Plasma_Nglyco", TRUE, FALSE)




# -------------------------------------------------------------------------------------
# CONVERT PROTIEN ID TO GENE SYSMBOL AND OVERLAP WITH EXOSOMES FROM WILLMS et al.
# -------------------------------------------------------------------------------------


intsec.list <- list(tif_full_name, exosomes)
names(intsec.list) <- c("TIF", "exosomes")

TIF_exo <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_TIF_Exosomes", TRUE, FALSE)




# -------------------------------------------------------------------------------------
# GET OVERLAP BETWEEN FULL TIF, POOLED NIF AND POOLED FIF.
# -------------------------------------------------------------------------------------


intsec.list <- list(tif_full_name, nif_pooled_id, fif_pooled_id)
names(intsec.list) <- c("TIF", "NIF", "FIF")

TIF_NIF_FIF <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_TIF_NIF_FIF", TRUE, FALSE)




# -------------------------------------------------------------------------------------
# GENE ONTOLOGY ENRICHMENT ANALTSIS OF OVERLAP BETWEEN TIF, NIF AND FIF: 
# -------------------------------------------------------------------------------------

proint_fifnif <- rbind(tif_full_fif, tif_full_nif)
proint_fifnif <- as.vector(unique(proint_fifnif[order(proint_fifnif$overlap) , ]))

univ_fifnif <- c(fif_pooled_id, nif_pooled_id, rownames(tifdata_small))
univ_fifnif <- unique(sort(univ_fifnif))

tif_nif_fif_GOobject <- GOobject("BP", univ_fifnif, proint_fifnif, GO_background) 

tif_nif_fif_GO <- TOPGO("BP", univ_fifnif, proint_fifnif, GO_background, 50) 




# GO_TERM GRAPHS

test.stat = new("elimCount", testStatistic = GOFisherTest, name = "elim test", cutOff = 0.01)
resultElim = getSigGroups(tif_nif_fif_GOobject, test.stat)

showSigOfNodes(tif_nif_fif_GOobject, score(resultElim), firstSigNodes = 5, useInfo ='all')
printGraph(tif_nif_fif_GOobject, resultElim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

barplot_func(tif_nif_fif_GO)

corrplot_func(tif_nif_fif_GO)





# -----------------------------------------------------------------------------------------------
# BATCH CORRECTION AND MDS PLOT
# -----------------------------------------------------------------------------------------------



# BATCH CORRECTION USING POOL, TUMOR PERCENTAGE AND TUMOR INFILTRATING LYMPHOCYTES.


# MODEL ON DIAGNOSIS
batch_corr_diag <- combat_corrections(tifdata_small, diag, pool)
batch_corr_diag_TILs <- combat_corrections(tifdata_small, diag, pool, TILs)


# MODEL ON RECEPTOR STATUS 
batch_corr_receptor<- combat_corrections(tifdata_small, receptor_status, pool, TILs)


# MULTIDIMENSIONAL SCALING, LABELED BY EITHER DIAGNOSIS OR RECEPTOR STATUS
myMDSplot(tifdata_small, diag, diag, colorcode_diag)
myMDSplot(batch_corr_diag_TILs, diag, diag, colorcode_diag)

myMDSplot(tifdata_small, receptor_status, receptor_status, colorcode_receptor)
myMDSplot(batch_corr_receptor, receptor_status, receptor_status, colorcode_receptor)



# -----------------------------------------------------------------------------------------------
# LABLE CHECK
# -----------------------------------------------------------------------------------------------

# COMPARE THE TWO TYPES OF MODELLING (DIAG; RECEPTOR STATUS):
lab_check_diag <- pvclust(batch_corr_diag_TILs, method.dist="euclidean", method.hclust="ward.D2", nboot=1000)
plot(lab_check_diag, hang=-1, labels=as.factor(diag))
pvrect(lab_check_diag, alpha=0.90)

lab_check_receptor_status <- pvclust(batch_corr_receptor, method.dist="euclidean", method.hclust="ward.D2", nboot=1000)
plot(lab_check_receptor_status, hang=-1, labels=as.factor(receptor_status))
pvrect(lab_check_receptor_status, alpha=0.90)



# -----------------------------------------------------------------------------------------------
# HIERARCHICAL CLUSTERING
# -----------------------------------------------------------------------------------------------

# heatmap colors in blue
heat.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)


# White spacer
my.spacer <- as.matrix(replicate(nrow(tifinfo), "white"))

# Color schemes for heatmap
my.TS <- get_colors(diag, colorcode_diag)
my.ER <- get_colors(ER, c("seashell", "grey60"))
my.PGR <- get_colors(PGR, c("seashell","grey60"))
my.GR <- get_colors(Gr, c("navyblue", "lightblue"))
my.TILS <- get_colors(TILs, c("orange", "yellow"))
#my.Age <- get_colors(AgeITV, c(brewer.pal(3,"Purples")))
#my.TP <- get_colors(tp, c(brewer.pal(4,"Reds")))
#my.AR <- get_colors(AR, c("seashell","grey60"))
#my.HER2 <- get_colors(as.factor(tifinfo$HER2), c(brewer.pal(4,"Greys")))
#my.pool <- get_colors(pool, c("olivedrab","bisque3", "azure3", "seagreen"))

my.cols <- cbind(my.TILS, my.spacer, my.GR, my.spacer, my.ER, my.spacer, my.PGR, my.spacer, my.TS)


mod_design <- model.matrix(~diag)
batch_corr <- ComBat(tifdata_small, pool, mod_design, par.prior=TRUE,prior.plots=FALSE)

pdf("proteinclusters.pdf")
heatmap.plus(as.matrix(scale(batch_corr , scale = FALSE)),col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labRow="", labCol="", ColSideColors=my.cols, margins = c(14,8), cexCol=1.2, cexRow = 1.3)
legend("bottomleft", legend = levels(diag), ncol=1, lty=c(1,1), lwd=c(3,3), cex=0.7, col=colorcode_diag, bty = "n")
legend("bottom", legend = c("High TILs", "Low TILs", "High Grade", "Low Grade", "ER+ & PGR+", "ER- & PGR-"), ncol=3, lty=c(1,1), lwd=c(3,3), cex=0.7, col=c("orange", "yellow", "navyblue", "lightblue", "grey60", "seashell"), bty = "n")
dev.off()






# ---------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING DIAGNOSIS.
# ---------------------------------------------------------------------------------------------


# CALCULATING DIFFERENCES BETWEEN LINEAR MODELS: LogFC > 1 and FDR < 0.5 OR LogFC < -1 and FDR < 0.5. 



# MAKE DESIGN MATRIX
diag_design <- model.matrix(~0+diag+pool)
DA_diag1 <- DA_all_contrasts(tifdata_small, diag_design, diag, "diag", 1, 0.05)

diag_design <- model.matrix(~0+diag+pool+TILs)
DA_diag2 <- DA_all_contrasts(tifdata_small, diag_design, diag, "diag", 1, 0.05)

#diag_design <- model.matrix(~0+diag+pool+TILs+Gr)
#DA_diag3 <- DA_all_contrasts(tifdata_small, diag_design, diag, "diag", 1, 0.05)


# ---------------------------------------------------------------------------------------------
# Writing out tables with DA proteins
# ---------------------------------------------------------------------------------------------

#write_out(rbind(DA_diag$`diagHER2-diagLumA`[[1]],DA_diag$`diagHER2-diagLumA`[[2]]), "HER2_LumA_DA")
#write_out(rbind(DA_diag$`diagHER2-diagTNBC`[[1]],DA_diag$`diagHER2-diagTNBC`[[2]]), "HER2_TNBC_DA")
#write_out(rbind(DA_diag$`diagLumA-diagTNBC`[[1]], DA_diag$`diagLumA-diagTNBC`[[2]]), "LumA_TNBC_DA")


# ---------------------------------------------------------------------------------------------
# MAKEING FULL DIAG DA DATASET
# ---------------------------------------------------------------------------------------------

#sig_d <- unique(sort(c(rownames(DA_diag$`diagHER2-diagLumA`[[1]]), rownames(DA_diag$`diagHER2-diagLumA`[[2]]), rownames(DA_diag$`diagHER2-diagTNBC`[[1]]), rownames(DA_diag$`diagHER2-diagTNBC`[[2]]), rownames(DA_diag$`diagLumA-diagTNBC`[[1]]), rownames(DA_diag$`diagLumA-diagTNBC`[[2]]))))
#sig_diag <- batch_corr_diag_TILs[rownames(batch_corr_diag_TILs) %in% sig_d, ]
#sig_diag <- tifdata_small[rownames(tifdata_small) %in% sig_d, ]
#sig_diag_scaled<- scale(sig_diag, center = TRUE, scale = FALSE)


# ---------------------------------------------------------------------------------------------
# H-CLUSTERING OF FULL DIAG DATASET
# ---------------------------------------------------------------------------------------------

#tree(sig_diag_scaled, 3, colorcode_diag, diag)


# hierarchical clustering
# heatmap colors in blue
#heat.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)


# HEATMAP.PLUS COLOR-MATRIX
#col_diag <- get_colors(diag, colorcode_diag)
#col_diag <- cbind(replicate(nrow(col_diag), "white"), col_diag)


# HEATMAP
#heatmap.plus(as.matrix(sig_diag_scaled), dendrogram="col", col= heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow="", ColSideColors=col_diag, margins = c(14,8))


# ---------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING HORMONE RECEPTORS
# ---------------------------------------------------------------------------------------------

# Make Design Matrix
ER_design <- model.matrix(~0+ER+pool)
DA_ER <- DA_all_contrasts(tifdata_small, ER_design, ER, "ER", 1, 0.05)

# Make Design Matrix
PGR_design <- model.matrix(~0+PGR+pool)
DA_PGR <- DA_all_contrasts(tifdata_small, PGR_design, PGR, "PGR", 1, 0.05)

# Make Design Matrix
AR_design <- model.matrix(~0+AR+pool)
DA_AR <- DA_all_contrasts(tifdata_small, AR_design, AR, "AR", 1, 0.05)

# Make Design Matrix
HER2_design <- model.matrix(~0+HER2+pool)
DA_HER2 <- DA_all_contrasts(tifdata_small, HER2_design, HER2, "HER2", 1, 0.05)

# ---------------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------------
#write_out(rbind(DA_ER$`ERERm-ERERp`[[1]],DA_ER$`ERERm-ERERp`[[2]]), "ER_DA_corrected")
#write_out(rbind(DA_PGR$`PGRPGRm-PGRPGRp`[[1]],DA_PGR$`PGRPGRm-PGRPGRp`[[2]]), "PGR_DA_corrected")
#write_out(rbind(DA_HER2$`HER2H0-HER2H2`[[1]], DA_HER2$`HER2H0-HER2H2`[[2]], DA_HER2$`HER2H1-HER2H2`[[2]]), "H01_H2_DA_corrected")
#write_out(rbind(DA_HER2$`HER2H0-HER2H3`[[1]], DA_HER2$`HER2H0-HER2H3`[[2]]), "H01_H3_DA_corrected_pool")



# ---------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING TILS AND TUMOUR GRADE
# ---------------------------------------------------------------------------------------------

#  Make Design Matix
TILs_design <- model.matrix(~0+TILs+pool)
DA_TILs <- DA_all_contrasts(tifdata_small, TILs_design, TILs, "TILs", 1, 0.05)

#  Make Design Matix
Gr_design <- model.matrix(~0+Gr+pool)
DA_Gr <- DA_all_contrasts(tifdata_small, Gr_design, Gr, "Gr", 1, 0.05)


# ---------------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------------
#write_out(rbind(DA_TILs$`TILshigh-TILslow`[[1]],DA_TILs$`TILshigh-TILslow`[[2]]), "TILS_pool_DA")



# ---------------------------------------------------------------------------------------------
# Overlap DA
# ---------------------------------------------------------------------------------------------


# HER2 vs Luminal A 

intsec.list <- list(rownames(rbind(DA_diag1$`diagHER2-diagLumA`[[1]],DA_diag1$`diagHER2-diagLumA`[[2]])), rownames(rbind(DA_diag2$`diagHER2-diagLumA`[[1]],DA_diag2$`diagHER2-diagLumA`[[2]])), rownames(rbind(DA_TILs$`TILshigh-TILslow`[[1]],DA_TILs$`TILshigh-TILslow`[[2]])), rownames(rbind(DA_ER$`ERERm-ERERp`[[1]], DA_ER$`ERERm-ERERp`[[2]])), rownames(rbind(DA_PGR$`PGRPGRm-PGRPGRp`[[1]], DA_PGR$`PGRPGRm-PGRPGRp`[[2]])),  rownames(rbind(DA_HER2$`HER2H0-HER2H3`[[1]], DA_HER2$`HER2H0-HER2H3`[[2]])))
names(intsec.list) <- c("HER2LumA", "HER2LumATILs", "TILs", "ER", "PGR", "HER2")

HER2_LumA_ov <- plot_upsetR(intsec.list, c("HER2LumA", "HER2LumATILs"), "Overlap_HER2_LumA", TRUE, FALSE)



# HER2 vs TNBC 

intsec.list <- list(rownames(rbind(DA_diag1$`diagHER2-diagTNBC`[[1]],DA_diag1$`diagHER2-diagTNBC`[[2]])), rownames(rbind(DA_diag2$`diagHER2-diagTNBC`[[1]],DA_diag2$`diagHER2-diagTNBC`[[2]])), rownames(rbind(DA_TILs$`TILshigh-TILslow`[[1]],DA_TILs$`TILshigh-TILslow`[[2]])),rownames(rbind(DA_ER$`ERERm-ERERp`[[1]], DA_ER$`ERERm-ERERp`[[2]])), rownames(rbind(DA_PGR$`PGRPGRm-PGRPGRp`[[1]], DA_PGR$`PGRPGRm-PGRPGRp`[[2]])),rownames(rbind(DA_HER2$`HER2H0-HER2H3`[[1]], DA_HER2$`HER2H0-HER2H3`[[2]])))
names(intsec.list) <- c("HER2TNBC", "HER2TNBCTILs", "TILs", "ER", "PGR", "HER2")

HER2_TNBC_ov <- plot_upsetR(intsec.list, "HER2TNBC", "Overlap_HER2_TNBC", TRUE, FALSE)




# Luminal A vs TNBC

intsec.list <- list(rownames(rbind(DA_diag1$`diagLumA-diagTNBC`[[1]],DA_diag1$`diagLumA-diagTNBC`[[2]])), rownames(rbind(DA_diag2$`diagLumA-diagTNBC`[[1]],DA_diag2$`diagLumA-diagTNBC`[[2]])), rownames(rbind(DA_TILs$`TILshigh-TILslow`[[1]],DA_TILs$`TILshigh-TILslow`[[2]])),rownames(rbind(DA_ER$`ERERm-ERERp`[[1]], DA_ER$`ERERm-ERERp`[[2]])), rownames(rbind(DA_PGR$`PGRPGRm-PGRPGRp`[[1]], DA_PGR$`PGRPGRm-PGRPGRp`[[2]])),rownames(rbind(DA_HER2$`HER2H0-HER2H3`[[1]], DA_HER2$`HER2H0-HER2H3`[[2]])))
names(intsec.list) <- c("LumATNBC", "LumATNBCTILs", "TILs", "ER", "PGR", "HER2")

LumA_TNBC_ov <- plot_upsetR(intsec.list, "HER2LumA", "Overlap_LumA_TNBC", TRUE, FALSE)





# ---------------------------------------------------------------------------------------
# GENE ONTOLOGY ENRICHMENT USIGN DIAGNOSIS.
# ---------------------------------------------------------------------------------------

univ_tif <- rownames(tifdata_small)
univ_tif <- unique(sort(univ_tif))


H_L_GO_up <- GO_visual(rownames(DA_diag$`diagHER2-diagLumA`[[1]]), univ_tif, GO_background, 30, "H_L_up", TRUE, TRUE)
H_L_GO_down <- GO_visual(rownames(DA_diag$`diagHER2-diagLumA`[[2]]), univ_tif, GO_background, 30, "H_L_down", TRUE, TRUE)
H_T_GO_up <- GO_visual(rownames(DA_diag$`diagHER2-diagTNBC`[[1]]), univ_tif, GO_background, 30, "H_T_up", TRUE, TRUE)
H_T_GO_down <- GO_visual(rownames(DA_diag$`diagHER2-diagTNBC`[[2]]), univ_tif, GO_background, 30, "H_T_down", TRUE, TRUE)
L_T_GO_up <- GO_visual(rownames(DA_diag$`diagLumA-diagTNBC`[[1]]), univ_tif, GO_background, 30, "L_T_up", TRUE, TRUE)
L_T_GO_down <- GO_visual(rownames(DA_diag$`diagLumA-diagTNBC`[[2]]), univ_tif, GO_background, 30, "L_T_down", TRUE, TRUE)




intsec.list <- list(rownames(DA_diag1$`diagHER2-diagLumA`[[2]]), rownames(DA_diag1$`diagHER2-diagTNBC`[[2]]), rownames(DA_diag1$`diagLumA-diagTNBC`[[2]]), rownames(DA_diag1$`diagHER2-diagLumA`[[1]]),rownames(DA_diag1$`diagHER2-diagTNBC`[[1]]),rownames(DA_diag1$`diagLumA-diagTNBC`[[1]]))
names(intsec.list) <- c("HER2LumADown", "HER2TNBCDown", "LumATNBCDown", "HER2LumAUp", "HER2TNBCUp", "LumATNBCUp")

Ov_Subtypes <- plot_upsetR(intsec.list, c("HER2TNBCDown", "HER2LumADown"), "Overlap_Subtypes", TRUE, FALSE)


# ---------------------------------------------------------------------------------------
# ODDS RATIO DA PROTEINS MICRO VESICLES VS. TIF  
# ---------------------------------------------------------------------------------------

exo_tif_ov <- intersect(tif_full_name, exosomes)
exo_tif_no_ov <- tif_full_name[!tif_full_name %in% exo_tif_ov]

sig_d_genes <- unique(sort(uniprot_to_name(sig_d)$name)) 

get_stats(sig_d_genes, exo_tif_no_ov, exo_tif_ov)




# ---------------------------------------------------------------------------------------
# PCA ANALYSIS USING DIAG
# ---------------------------------------------------------------------------------------


# ROTATIONS OF PROTEINS IN N-DIMENSIONAL SPACE
#rot_diag <- rotation_func(batch_corr_diag_TILs)


# BOXPLOT TO IDENTIFY MOST IMPORTANT PCs:
#boxplots(batch_corr_diag_TILs, diag, colorcode_diag)


# EXTRACTING THE ROTATION FROM MOST INTERESTING PCs
#PCA_diag <- sort(c(as.character(rownames(head(rot_diag[order(rot_diag$PC1),], n=50))), as.character(rownames(tail(rot_diag[order(rot_diag$PC1),], n=50))), as.character(rownames(head(rot_diag[order(rot_diag$PC2),], n=50))), as.character(rownames(tail(rot_diag[order(rot_diag$PC2),], n=50))), as.character(rownames(head(rot_diag[order(rot_diag$PC4),], n=50))), as.character(rownames(tail(rot_diag[order(rot_diag$PC4),], n=50)))))
#PCA_data_diag <- batch_corr_diag_TILs[rownames(batch_corr_diag_TILs) %in% PCA_diag, ]


# OVERLAP PCA RESULTS WITH DA RESULTS 
#venn <- venn.diagram(list(A=sig_d, B=PCA_diag), category.names = c("diag", "PCA_diag"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 2, cex=1.5, fill=rainbow(2))
#grid.draw(venn)


# OVERLAP PCA RESULTS WITH DA RESULTS - GET GENE INFO
#TH_up_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(T_H[[1]]$up), PCA_diag), ]
#TH_down_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(T_H[[2]]$down), PCA_diag), ]

#TL_up_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(T_L[[1]]$up), PCA_diag), ]
#TL_down_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(T_L[2]$down), PCA_diag), ]

#LH_up_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(L_H[[1]]$up), PCA_diag), ]
#LH_down_DA_PCA <- hugo[hugo$Accession %in% intersect(as.character(L_H[[2]]$down), PCA_diag), ]




# ---------------------------------------------------------------------------------------
# LASSO USING DIAGNOSIS
# ---------------------------------------------------------------------------------------


group_diag <- as.integer(diag)

LASSO_BC1 <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_diag, TRUE), LASSO_protein(40, batch_corr_diag, group_diag, TRUE), LASSO_protein(99, batch_corr_diag, group_diag, TRUE), LASSO_protein(640, batch_corr_diag, group_diag, TRUE), LASSO_protein(6565, batch_corr_diag, group_diag, TRUE)))
LASSO_BC2 <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag_TILs, group_diag, TRUE), LASSO_protein(40, batch_corr_diag_TILs, group_diag, TRUE), LASSO_protein(99, batch_corr_diag_TILs, group_diag, TRUE), LASSO_protein(640, batch_corr_diag_TILs, group_diag, TRUE), LASSO_protein(6565, batch_corr_diag_TILs, group_diag, TRUE)))
LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_diag, TRUE), LASSO_protein(40, tifdata_small, group_diag, TRUE), LASSO_protein(99, tifdata_small, group_diag, TRUE), LASSO_protein(640, tifdata_small, group_diag, TRUE), LASSO_protein(6565, tifdata_small, group_diag, TRUE)))
LASSO_diag <- intersect(LASSO_BC1, LASSO_BC2)
#write_out(LASSO_diag, "LASSO_subtypes_corrected_pool_TILs")

# ---------------------------------------------------------------------------------------
# LASSO USING ER, PGR and TILS
# ---------------------------------------------------------------------------------------

# LASSO ER
group_ER <- as.integer(ER)
LASSO_BC <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_ER, FALSE), LASSO_protein(40, batch_corr_diag, group_ER,  FALSE), LASSO_protein(99, batch_corr_diag, group_ER,   FALSE), LASSO_protein(640, batch_corr_diag, group_ER,  FALSE), LASSO_protein(6565, batch_corr_diag, group_ER,  FALSE)))
LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_ER,  FALSE), LASSO_protein(40, tifdata_small, group_ER,  FALSE), LASSO_protein(99, tifdata_small, group_ER,  FALSE), LASSO_protein(640, tifdata_small, group_ER,  FALSE), LASSO_protein(6565, tifdata_small, group_ER,  FALSE)))
LASSO_ER <- unique(sort(c(LASSO_BC, LASSO_NonBC)))
#write_out(LASSO_ER, "LASSO_ER_corrected_pool")

# ---------------------------------------------------------------------------------------

# LASSO PGR
group_PGR <- as.integer(PGR)
LASSO_BC <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_PGR, FALSE), LASSO_protein(40, batch_corr_diag, group_PGR, FALSE), LASSO_protein(99, batch_corr_diag, group_PGR, FALSE), LASSO_protein(640, batch_corr_diag, group_PGR,  FALSE), LASSO_protein(6565, batch_corr_diag, group_PGR,  FALSE)))
LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_PGR, FALSE), LASSO_protein(40, tifdata_small, group_PGR, FALSE), LASSO_protein(99, tifdata_small, group_PGR, FALSE), LASSO_protein(640, tifdata_small, group_PGR,  FALSE), LASSO_protein(6565, tifdata_small, group_PGR, FALSE)))
LASSO_PGR <- unique(sort(c(LASSO_BC, LASSO_NonBC)))
#write_out(LASSO_PGR, "LASSO_PGR_corrected_pool")

# ---------------------------------------------------------------------------------------

# LASSO HER2
HER2 <- as.factor(ifelse(HER2 %in% c("H2", "H3"), "H23", as.character(HER2)))
group_HER2 <- as.integer(HER2)
LASSO_BC <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_HER2, TRUE), LASSO_protein(40, batch_corr_diag, group_HER2, TRUE), LASSO_protein(99, batch_corr_diag, group_HER2, TRUE), LASSO_protein(640, batch_corr_diag, group_HER2, TRUE), LASSO_protein(6565, batch_corr_diag, group_HER2, TRUE)))
LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_HER2, TRUE), LASSO_protein(40, tifdata_small, group_HER2, TRUE), LASSO_protein(99, tifdata_small, group_HER2, TRUE), LASSO_protein(640, tifdata_small, group_HER2, TRUE), LASSO_protein(6565, tifdata_small, group_HER2, TRUE)))
LASSO_HER2 <- unique(sort(c(LASSO_BC, LASSO_NonBC)))
#write_out(LASSO_HER2, "LASSO_HER2_corrected_pool")

# ---------------------------------------------------------------------------------------

# LASSO TILS
group_TILs <- as.integer(TILs)
LASSO_BC <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_TILs, TRUE), LASSO_protein(40, batch_corr_diag, group_TILs, TRUE), LASSO_protein(99, batch_corr_diag, group_TILs, TRUE), LASSO_protein(640, batch_corr_diag, group_TILs, TRUE), LASSO_protein(6565, batch_corr_diag, group_TILs, TRUE)))
LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_TILs, TRUE), LASSO_protein(40, tifdata_small, group_TILs, TRUE), LASSO_protein(99, tifdata_small, group_TILs, TRUE), LASSO_protein(640, tifdata_small, group_TILs, TRUE), LASSO_protein(6565, tifdata_small, group_TILs, TRUE)))
LASSO_TILs <- unique(sort(c(LASSO_BC, LASSO_NonBC)))
#write_out(LASSO_TILs, "LASSO_TILs_corrected_pool")

# ---------------------------------------------------------------------------------------

# LASSO Grade
#group_Gr <- as.integer(Gr)
#LASSO_BC <- Reduce(intersect, list(LASSO_protein(1, batch_corr_diag, group_Gr, TRUE), LASSO_protein(40, batch_corr_diag, group_Gr, TRUE), LASSO_protein(99, batch_corr_diag, group_Gr, TRUE), LASSO_protein(640, batch_corr_diag, group_Gr, TRUE), LASSO_protein(6565, batch_corr_diag, group_Gr, TRUE)))
#LASSO_NonBC <- Reduce(intersect, list(LASSO_protein(1, tifdata_small, group_Gr, TRUE), LASSO_protein(40, tifdata_small, group_Gr, TRUE), LASSO_protein(99, tifdata_small, group_Gr, TRUE), LASSO_protein(640, tifdata_small, group_Gr, TRUE), LASSO_protein(6565, tifdata_small, group_Gr, TRUE)))
#LASSO_Gr <- unique(sort(c(LASSO_BC, LASSO_NonBC)))



# ---------------------------------------------------------------------------------------
# RANDOM FOREST USING DIAGNOSIS
# ---------------------------------------------------------------------------------------



tifdata_small_RF <- t(tifdata_small)
batch_corr_diag_TILs_RF <- t(batch_corr_diag_TILs)
batch_corr_diag_RF <- t(batch_corr_diag)

RF_BC <- Reduce(intersect, list(my_forest(1, batch_corr_diag_TILs_RF, diag), my_forest(40, batch_corr_diag_TILs_RF, diag), my_forest(99, batch_corr_diag_TILs_RF, diag), my_forest(640, batch_corr_diag_TILs_RF, diag), my_forest(6565, batch_corr_diag_TILs_RF, diag)))
RF_BC_conver <- my_forest_conver(1, batch_corr_diag_TILs_RF, diag)

### Converges for TNBC: 0.16 and LumA: 0.05, not for HER2 ###
### "Q7Z3D4" "Q969E4" "Q9BU02" "Q9HB07" ###

RF_NonBC <- Reduce(intersect, list(my_forest(1, tifdata_small_RF, diag), my_forest(40, tifdata_small_RF, diag), my_forest(99, tifdata_small_RF, diag), my_forest(640, tifdata_small_RF, diag), my_forest(6565, tifdata_small_RF, diag)))
RF_NonBC_conver <- my_forest_conver(1, tifdata_small_RF, diag)

### Converges for TNBC: 0.25 and LumA: 0.05, not for HER2 ###
### "Q9HB07", "E9PFK9" "Q9BU02" "Q9HB07" "Q9NYQ6" ###

RF_diag <- unique(sort(c(RF_BC, RF_NonBC)))
# RF_diag <- c("Q9HB07", "E9PFK9", "Q9BU02", "Q9NYQ6", "Q7Z3D4", "Q969E4")

# ---------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------
#write_out(RF_diag, "RF_subtypes_corrected_pool_TILS")



# ---------------------------------------------------------------------------------------
# RANDOM FOREST USING ER, PGR and TILs
# ---------------------------------------------------------------------------------------



# ER
RF_BC <- Reduce(intersect, list(my_forest(1, batch_corr_diag_RF, ER), my_forest(40, batch_corr_diag_RF, ER), my_forest(99, batch_corr_diag_RF, ER), my_forest(640, batch_corr_diag_RF, ER), my_forest(6565, batch_corr_diag_RF, ER)))
RF_BC_conver <- my_forest_conver(1, batch_corr_diag_RF, ER)
### Converges, ERm: 0.13 and ERp: 0.05 ###
### "E9PFK9", "Q7Z3D4" "Q969E4" "Q9BU02", "Q9HB07" "Q9NU22" ###

RF_NonBC <- Reduce(intersect, list(my_forest(1, tifdata_small_RF, ER), my_forest(40, tifdata_small_RF, ER), my_forest(99, tifdata_small_RF, ER), my_forest(640,  tifdata_small_RF, ER), my_forest(6565,  tifdata_small_RF, ER)))
RF_NonBC_conver <- my_forest_conver(1, tifdata_small_RF, ER)
### Converges, ERm: 0.26 and ERp: 0.05 ###
### "Q9NYQ6" "E9PFK9" "P50895" "J3KNL6" "O95049" "Q9BU02" "Q9HB07" ###

RF_ER <- unique(sort(c(RF_BC, RF_NonBC)))
# RF_ER <- c("E9PFK9", "Q7Z3D4", "Q969E4", "Q9BU02", "Q9HB07", "Q9NU22","Q9NYQ6", "P50895", "J3KNL6", "O95049")

# ---------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------
#write_out(RF_ER, "RF_ER_corrected_pool")

# ---------------------------------------------------------------------------------------


# PGR

RF_BC <- Reduce(intersect, list(my_forest(1, batch_corr_diag_RF, PGR), my_forest(40, batch_corr_diag_RF, PGR), my_forest(99, batch_corr_diag_RF, PGR), my_forest(640, batch_corr_diag_RF, PGR), my_forest(6565, batch_corr_diag_RF, PGR)))
RF_BC_conver <- my_forest_conver(1, batch_corr_diag_RF, PGR)
###  Borderline converges, PGRm: 0.16 and PGRp: 0.33 ###
###  "B7ZLZ0" "E7EV62" "P50895", "E9PGQ4", "Q6UXD5" "Q96MY1" ###

RF_NonBC <- Reduce(intersect, list(my_forest(1, tifdata_small_RF, PGR), my_forest(40, tifdata_small_RF, PGR), my_forest(99, tifdata_small_RF, PGR), my_forest(640,  tifdata_small_RF, PGR), my_forest(6565,  tifdata_small_RF, PGR)))
RF_NonBC_conver <- my_forest_conver(1, tifdata_small_RF, PGR)
###  Borderline converges, PGRm: 0.21 and PGRp: 0.33 ###
###  "E9PGQ4", "P50895" , "E9PFK9" "E9PFN5" "Q9BU02" "Q9P206" ###


RF_PGR <- unique(sort(c(RF_BC, RF_NonBC)))
### RF_PGR <- c("B7ZLZ0", "E7EV62", "P50895", "E9PGQ4", "Q6UXD5", "Q96MY1", "E9PFK9", "E9PFN5", "Q9BU02", "Q9P206") ###
### write_out(RF_PGR, "RF_PGR_corrected_pool") ###


# HER2

#RF_BC <- Reduce(intersect, list(my_forest(1, batch_corr_diag_RF, HER2), my_forest(40, batch_corr_diag_RF, HER2), my_forest(99, batch_corr_diag_RF, HER2), my_forest(640, batch_corr_diag_RF, HER2), my_forest(6565, batch_corr_diag_RF, HER2)))
#RF_BC_conver <- my_forest_conver(1, batch_corr_diag_RF, HER2)

###  Borderline converges, TILshigh: 0.14 and TILslow: 0.31 ###
###  "B9A018" "O43505" "Q9NTJ3" "O43570" "P35052" "Q96IZ0" "Q96NL8" "Q9UJA5" ###

#RF_NonBC <- Reduce(intersect, list(my_forest(1, tifdata_small_RF, HER2), my_forest(40, tifdata_small_RF, HER2), my_forest(99, tifdata_small_RF, HER2), my_forest(640,  tifdata_small_RF, HER2), my_forest(6565,  tifdata_small_RF, HER2)))
#RF_NonBC_conver <- my_forest_conver(1, tifdata_small_RF, HER2)

### NOT converged, TILshigh: 0.19 and TILslow: 0.46 ###
###  "E7EVV3" "G3XAP6" "P35052", "B9A018", "O95466", "P16885" ###


#RF_TILs <- unique(sort(c(RF_BC, RF_NonBC)))
#RF_TILs <- c("B9A018", "O43505", "Q9NTJ3", "O43570", "P35052", "Q96IZ0", "Q96NL8", "Q9UJA5", "E7EVV3", "G3XAP6", "O95466", "P16885")


# ---------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------

#write_out(RF_TILs, "RF_TILs_corrected_pool")
# ---------------------------------------------------------------------------------------



# TILS
RF_BC <- Reduce(intersect, list(my_forest(1, batch_corr_diag_RF, TILs), my_forest(40, batch_corr_diag_RF, TILs), my_forest(99, batch_corr_diag_RF, TILs), my_forest(640, batch_corr_diag_RF, TILs), my_forest(6565, batch_corr_diag_RF, TILs)))
RF_BC_conver <- my_forest_conver(1, batch_corr_diag_RF, TILs)

### Borderline converges, TILshigh: 0.14 and TILslow: 0.31 ###
### "B9A018" "O43505" "Q9NTJ3" "O43570" "P35052" "Q96IZ0" "Q96NL8" "Q9UJA5" ###

RF_NonBC <- Reduce(intersect, list(my_forest(1, tifdata_small_RF, TILs), my_forest(40, tifdata_small_RF, TILs), my_forest(99, tifdata_small_RF, TILs), my_forest(640,  tifdata_small_RF, TILs), my_forest(6565,  tifdata_small_RF, TILs)))
RF_NonBC_conver <- my_forest_conver(1, tifdata_small_RF, TILs)

### NOT converged, TILshigh: 0.19 and TILslow: 0.46 ###
### "E7EVV3" "G3XAP6" "P35052", "B9A018", "O95466", "P16885" ###


RF_TILs <- unique(sort(c(RF_BC, RF_NonBC)))
#RF_TILs <- c("B9A018", "O43505", "Q9NTJ3", "O43570", "P35052", "Q96IZ0", "Q96NL8", "Q9UJA5", "E7EVV3", "G3XAP6", "O95466", "P16885") 


# ---------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# ---------------------------------------------------------------------------------------

#write_out(RF_TILs, "RF_TILs_corrected_pool")



# ---------------------------------------------------------------------------------------
# Overlap results DA, LASSO and RF
# ---------------------------------------------------------------------------------------

setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_proteomics/lists_and_plots_10_01_2018/Tables/DA/Subtypes")
HER2_LumA_DA <- read.table("HER2_LumA_DA_corrected_pool.txt", header=TRUE)
HER2_TNBC_DA <- read.table("HER2_TNBC_DA_corrected_pool.txt", header = TRUE)
LumA_TNBC_DA <- read.table("LumA_TNBC_DA_corrected_pool.txt", header = TRUE)
HER2_LumA_DA_BC <- read.table("Best_DA_cand_HER2_LumA.txt", header=TRUE)
HER2_TNBC_DA_BC <- read.table("Best_DA_cand_HER2_TNBC.txt", header = TRUE)
LumA_TNBC_DA_BC <- read.table("Best_DA_cand_LumA_TNBC.txt", header = TRUE)

setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_proteomics/lists_and_plots_10_01_2018/Tables/LASSO")
Subtypes_LASSO <- read.table("LASSO_subtypes_corrected_pool_TILs.txt", header=TRUE)


setwd("~/Desktop/Thilde/MS_MS_TIF_analysis_2014_2015/TIF_proteomics/lists_and_plots_10_01_2018/Tables/RF")
Subtypes_RF <- read.table("RF_subtypes_corrected_pool_TILS.txt", header=TRUE)



intsec.list <- list(as.character(HER2_LumA_DA$Accession), as.character(HER2_LumA_DA_BC$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "BestCandDA", "LASSO", "RF")
Ov_HER2_LumA_DA_LA_RF <- plot_upsetR(intsec.list, c("BestCandDA"), "Best_cand_HER2_LumA", TRUE, FALSE)

intsec.list <- list(as.character(HER2_TNBC_DA$Accession), as.character(HER2_TNBC_DA_BC$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "BestCandDA", "LASSO", "RF")
Ov_HER2_TNBC_DA_LA_RF <- plot_upsetR(intsec.list, c("BestCandDA"), "Best_cand_HER2_TNBC", TRUE, FALSE)

intsec.list <- list(as.character(LumA_TNBC_DA$Accession), as.character(LumA_TNBC_DA_BC$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "BestCandDA", "LASSO", "RF")
Ov_LumA_TNBC_DA_LA_RF <- plot_upsetR(intsec.list, c("BestCandDA"), "Best_cand_LumA_TNBC", TRUE, FALSE)



# ---------------------------------------------------------------------------------------
# PATHWAY ANALYSIS
# ---------------------------------------------------------------------------------------

univ_tif_enz <- u2e.map[u2e.map$uniprot %in% rownames(tifdata_small), ]




# ---------------------------------------------------------------------------------------
# PATHWAY ANALYSIS USING DIAGNOSIS
# ---------------------------------------------------------------------------------------



H_L_up_pw <- enrich_pathway(univ_tif_enz, u2e.map, as.character(rownames(DA_diag$`diagHER2-diagLumA`[[1]])), "H_L_up_pw", TRUE)
H_L_down_pw <- enrich_pathway(univ_tif_enz, u2e.map, character(rownames(DA_diag$`diagHER2-diagLumA`[[2]])), "H_L_down_pw", TRUE)

H_T_up_pw <- enrich_pathway(univ_tif_enz, u2e.map, as.character(rownames(DA_diag$`diagHER2-diagTNBC`[[1]])), "H_T_up_pw", TRUE)
H_T_down_pw <- enrich_pathway(univ_tif_enz, u2e.map, as.character(rownames(DA_diag$`diagHER2-diagTNBC`[[2]])), "H_T_down_pw", TRUE)

L_T_up_pw <- enrich_pathway(univ_tif_enz, u2e.map, as.character(rownames(DA_diag$`diagLumA-diagTNBC`[[1]])), "L_H_up_pw", TRUE)
L_T_down_pw <- enrich_pathway(univ_tif_enz, u2e.map, as.character(rownames(DA_diag$`diagLumA-diagTNBC`[[2]])), "L_H_down_pw", TRUE)



entrez.data <- u2e.map[u2e.map$uniprot %in% as.character(rownames(DA_diag$`diagHER2-diagLumA`[[1]])), ]
my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.05, as.character(univ_tif_enz$entrez), minGSSize = 4, readable = T)





# ---------------------------------------------------------------------------------------
# DE PROTEINS OVERLAP WITH DE CYTOKINES
# ---------------------------------------------------------------------------------------

cytokines <- c("PDGF", "IL1RA", "IL7", "VEGF", "FGF", "IL10", "IL13", "P10", "IL12", "IL1B", "RANTES")
cytokines_accession <- hugo[grep("PDGF|IL1RA|IL7|VEGF|FGF|IL10|IL13|P10|IL12|IL1B|RANTES", hugo$name), ]
cytokines_accession <- sort(as.character(cytokines_accession$Accession))


cytokines_accession2 <- c("P01127", "P04085", "P09619", "P16234", "Q9NRA1",  "Q9GZP0", "P18510", "P13232", "P15692", "P49765", "P35968", "P49767", "Q9NRA1", "O43915", "P17948", "O15520", "P09038", "P55075", "P05230", "Q9GZV9", "P21781", "P31371", "Q92913", "P01584", "O95750", "O76093", "P12034", "P08620", "P22301", "Q6FGW4", 	"P35225", "P29460", "P13501")



ov_d_cyto <- sig_d[sig_d %in% cytokines_accession2]
#ov_pca_d_cyto <- PCA_diag[PCA_diag %in% cytokines_accession2]




# ---------------------------------------------------------------------------------------
# OVERLAP WITH PAM50 GENES
# ---------------------------------------------------------------------------------------


pam50_full <- hugo[hugo$name %in% pam50, ]
setdiff(pam50, unique(sort(pam50_full$name)))
pam50_accession <- c(as.character(pam50_full$Accession), "Q9Y5N6")

# PAM50 GENES IN PROTEOMICS SET 
ident_TIF <- intersect(unique(sort(as.character(tifdata$X1))), pam50_accession)
### 32 ###

ident_TIF_small <- intersect(unique(sort(as.character(rownames(tifdata_small)))), pam50_accession)
### 25 ###


# DIAGNOSIS
pam50_ov_diag_DA <- intersect(pam50_accession, sig_d)
pam50_ov_diag_DA <- hugo[hugo$Accession %in% pam50_ov_diag_DA, ]

#pam50_ov_diag_PCA <- intersect(pam50_accession, PCA_diag)
#pam50_ov_diag_PCA <- hugo[hugo$Accession %in% pam50_ov_diag_PCA, ]

pam50_ov_diag_lasso <- intersect(pam50_accession, rownames(LASSO_diag))
pam50_ov_diag_lasso <- hugo[hugo$Accession %in% pam50_ov_diag_lasso, ]

pam50_ov_diag_RF <- intersect(pam50_accession, rownames(RF_diag))
pam50_ov_diag_RF <- hugo[hugo$Accession %in% pam50_ov_diag_RF, ]



# ---------------------------------------------------------------------------------------
# OVERLAP WITH OTHER DATASETS SUBTYPES
# ---------------------------------------------------------------------------------------

genes_DA_LASSO_RF <- uniprot_to_name(c(sig_d, rownames(LASSO_diag), rownames(RF_diag)))

sig_d_genes <- unique(sort(as.character(genes_DA_LASSO_RF[genes_DA_LASSO_RF$Accession %in% sig_d,]$name)))
LASSO_diag_genes <- unique(sort(as.character(genes_DA_LASSO_RF[genes_DA_LASSO_RF$Accession %in% rownames(LASSO_diag),]$name)))
RF_diag_genes <- unique(sort(as.character(genes_DA_LASSO_RF[genes_DA_LASSO_RF$Accession %in% rownames(RF_diag),]$name)))


# Overlap with Tyanova and Angelo
ov_DA_TyanovaGST <- intersect(TyanovaGST$V1, sig_d_genes)
ov_DA_TyanovaDEG <- intersect(TyanovaDEG$V1, sig_d_genes)
ov_DA_AngeloDEG <- intersect(AngeloDEG$name, sig_d_genes)

ov_LASSO_TyanovaGST <- intersect(TyanovaGST$V1, LASSO_diag_genes)
ov_LASSO_TyanovaDEG <- intersect(TyanovaDEG$V1, LASSO_diag_genes)
ov_LASSO_AngeloDEG <- intersect(AngeloDEG$name, LASSO_diag_genes)

ov_RF_TyanovaGST <- intersect(TyanovaGST$V1, RF_diag_genes)
ov_RF_TyanovaDEG <- intersect(TyanovaDEG$V1, RF_diag_genes)
ov_RF_AngeloDEG <- intersect(AngeloDEG$name, RF_diag_genes)





# ---------------------------------------------------------------------------------------
# CORRELATIONS - DIAG
# ---------------------------------------------------------------------------------------


# Correlations

# TNBC vs LUMA
TL_nwk <- c(as.character(T_L[[1]]$up), as.character(T_L[[2]]$down))
TL_dir <- data.frame(TL_nwk, c(replicate(length(T_L[[1]]$up), "up"), replicate(length(T_L[[2]]$down), "down")))
colnames(TL_dir) <- c("first", "direction")
TL_nwk <- batch_corr_diag_TILs[rownames(batch_corr_diag_TILs) %in% TL_nwk, ]  

# TNBC vs HER2
TH_nwk <- c(as.character(T_H[[1]]$up), as.character(T_H[[2]]$down))
TH_dir <- data.frame(TH_nwk, c(replicate(length(T_H[[1]]$up), "up"), replicate(length(T_H[[2]]$down), "down")))
colnames(TH_dir) <- c("first", "direction")
TH_nwk <- batch_corr_diag_TILs[rownames(batch_corr_diag_TILs) %in% TH_nwk, ]

# LUMA vs HER2
LH_nwk <- c(as.character(L_H[[1]]$up), as.character(L_H[[2]]$down))
LH_dir <- data.frame(LH_nwk, c(replicate(length(L_H[[1]]$up), "up"), replicate(length(L_H[[2]]$down), "down")))
colnames(LH_dir) <- c("first", "direction")
LH_nwk <- batch_corr_diag_TILs[rownames(batch_corr_diag_TILs) %in% LH_nwk, ]


# Create network
T_L_nwk <- my_networks(TL_nwk, TL_dir, 0.75, 0.0001)
T_H_nwk <- my_networks(TH_nwk, TH_dir, 0.75, 0.0001)
L_H_nwk <- my_networks(LH_nwk, LH_dir, 0.75, 0.0001)

# Make a map to get genename of uniprot identifyer
diag_map <- uniprot_to_name(sig_d)

T_L_nwk <- get_nodes(diag_map, T_L_nwk, "TL")
T_H_nwk <- get_nodes(diag_map, T_H_nwk, "TH")
L_H_nwk <- get_nodes(diag_map, L_H_nwk, "LH")



