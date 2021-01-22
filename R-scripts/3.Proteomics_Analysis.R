# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




                                                                                                    ### ANALYSIS OF PROTEOMICS DATA ###




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my.wd <- ""

setwd(paste0(my.wd, "/Data/"))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET UNIPROT ID - GENE ID
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
tif_full_id <- as.character(sort(rownames(tifdata_small)))
tif_full  <- unique(uniprot_to_name(tif_full_id))
tif_full_name <- unique(tif_full$name)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERLAP TIF PROTEINS WITH PLASMA, EXOSOMES, SECRETED AND INDEPENDENT TIF SET
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Gene Symbol
intsec.list <- list(SecretedBoerG, tif_full_name, ExosomesAllG, PlasmaG, PlasmaMPG, TIFYatesG)
names(intsec.list) <- c("Secreted", "TIF", "Exosomes", "Plasma", "PlasmaMP", "TIFYates")

# Make upsetplot
colov <- c("#08605F", "#177E89", "#598381","#8E936D", "#5D737E", "#394053")
highly_expressed <- plot_upsetR(intsec.list, names(intsec.list), "test2", colov, TRUE, FALSE)




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GET OVERLAP BETWEEN FULL TIF, POOLED NIF AND POOLED FIF.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


intsec.list <- list(tif_full_name, nif_pooled_id, fif_pooled_id)
names(intsec.list) <- c("TIF", "NIF", "FIF")

TIF_NIF_FIF <- plot_upsetR(intsec.list, names(intsec.list), "Overlap_TIF_NIF_FIF", colov[c(1,2)], TRUE, FALSE)

# Plot venn diagram
venn <- venn.diagram(list(A=tif_full_id, B=nif_pooled_id, C=fif_pooled_id), category.names = c("TIF", "NIF_pooled", "FIF_pooled"), filename=NULL, lwd = 0.7, fill=rainbow(3), sub.cex = 2, cat.cex= 2, cex=1.5)
grid.draw(venn)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GENE ONTOLOGY ENRICHMENT ANALTSIS OF OVERLAP BETWEEN TIF, NIF AND FIF:
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

tif_nif_fif <- unique(Reduce(intersect, list(tif_full_id, nif_pooled_id, fif_pooled_id)))
tif_nif <- unique(intersect(tif_full_id, nif_pooled_id))
tif_fif <- unique(intersect(tif_full_id, fif_pooled_id))
nif_fif <- unique(intersect(nif_pooled_id, fif_pooled_id)) 
nif_fif <- gsub(pattern = "-[0-9]", "", nif_fif)

# Protein universe from TIF, NIF AND FIF
prot_univ <- unique(sort(c(fif_pooled_id, nif_pooled_id, tif_full_id)))

# GOobject
tif_nif_fif_GOobject <- GOobject("BP", prot_univ, tif_nif_fif, GO_background) 

# Enrich for GO terms
tif_nif_fif_GO <- TOPGO("BP", prot_univ, tif_nif_fif, GO_background, 50) 



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Test significance of enrichment
test.stat = new("elimCount", testStatistic = GOFisherTest, name = "elim test", cutOff = 0.01)
resultElim = getSigGroups(tif_nif_fif_GOobject, test.stat)

showSigOfNodes(tif_nif_fif_GOobject, score(resultElim), firstSigNodes = 5, useInfo ='all')
printGraph(tif_nif_fif_GOobject, resultElim, firstSigNodes = 5, fn.prefix = "tGO", useInfo = "all", pdfSW = TRUE)

# Plots
barplot_func(tif_nif_fif_GO)
corrplot_func(tif_nif_fif_GO)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PATHWAY ANALYSIS - OVERLAP TIF, NIF, FIF
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Convert IDs using biomaRt
ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
univ_enz <- unique(sort(as.character(getBM(attributes=c("uniprotswissprot", "entrezgene_id"), filters="uniprotswissprot", values=prot_univ, mart=ensembl)$entrezgene_id)))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Enrich for pathways
pw_tif_nif_fif <- enrich_pathway(univ_enz, tif_nif_fif, "tif_nif_fif_overlap", my.plot = TRUE)
pw_tif_nif <- enrich_pathway(univ_enz, tif_nif, "tif_nif_overlap", my.plot = TRUE)
pw_tif_fif <- enrich_pathway(univ_enz, tif_fif, "tif_fif_overlap", my.plot = TRUE)
pw_nif_fif <- enrich_pathway(univ_enz, nif_fif, "nif_fif_overlap", my.plot = TRUE)




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BATCH CORRECTION AND MDS PLOT
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# BATCH CORRECTION USING POOL, TUMOR PERCENTAGE AND TUMOR INFILTRATING LYMPHOCYTES.


# MODEL ON DIAGNOSIS and TILs
batch_corr_diag <- combat_corrections(tifdata_small, pool, diag)
batch_corr_TILs <- combat_corrections(tifdata_small, pool, TILs)


# MULTIDIMENSIONAL SCALING, LABELED BY EITHER DIAGNOSIS OR TIL STATUS
myMDSplot(tifdata_small, diag, "", colorcode_diag)
myMDSplot(batch_corr_diag, diag, "", colorcode_diag)
myMDSplot(batch_corr_TILs, TILs, "", brewer.pal(3, "Blues")[-1])



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LABLE CHECK
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# COMPARE THE TWO TYPES OF MODELLING (DIAG; RECEPTOR STATUS):
lab_check_diag <- pvclust(batch_corr_diag_TILs, method.dist="euclidean", method.hclust="ward.D2", nboot=1000)
plot(lab_check_diag, hang=-1, labels=as.factor(diag))
pvrect(lab_check_diag, alpha=0.90)

lab_check_receptor_status <- pvclust(batch_corr_receptor, method.dist="euclidean", method.hclust="ward.D2", nboot=1000)
plot(lab_check_receptor_status, hang=-1, labels=as.factor(receptor_status))
pvrect(lab_check_receptor_status, alpha=0.90)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Kmeans
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Get optimal number of clusters
#optimal_nc(batch_corr_diag)

set.seed(25)
#C2 <- as.factor(paste0("C",as.factor(as.character(kmeans(t(batch_corr_diag), 2, iter.max = 100)$cluster))))
C3 <- paste0("C",as.factor(as.character(kmeans(t(batch_corr_diag), 3, iter.max = 300)$cluster)))
C3 <- factor(C3, levels = c("C3", "C2", "C1"))



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# HIERARCHICAL CLUSTERING
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# heatmap colors in blue
heat.cols <- colorRampPalette(brewer.pal(9, "YlGnBu"))(n = 300)
heat.cols <- viridis(n=100)

# White spacer
my.spacer <- as.matrix(replicate(nrow(tifinfo), "white"))

# Color schemes for heatmap
my.ER <- get_colors(ER, c("#FFFCF7", "grey60"))
my.PGR <- get_colors(PGR, c("#FFFCF7", "grey60"))
my.HER2 <- get_colors(HER2, c("#FFFCF7", brewer.pal(4,"Greys")[2:4]))
my.GR <- get_colors(Gr, c("#032A63","#0071AA"))
my.TILS <- get_colors(TILs, c("#032A63","#0071AA"))
my.ST <- get_colors(diag, c("#FAA275","#DDDC71", "#BF5454"))
my.C3 <- get_colors(C3, viridis(3, begin = 0.2, end =0.9))

# Color sheme combined
my.cols <- cbind(my.ST, my.ER, my.PGR, my.HER2, my.TILS, my.GR)

# Plot heatmap
pdf("proteinclusters.pdf")
heatmap.plus(as.matrix(scale(batch_corr_diag)),col=heat.cols, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", Rowv = NA, labRow="", labCol="", ColSideColors=my.cols, margins = c(14,8), cexCol=1.2, cexRow = 1.3)
legend("bottomleft", legend = c(levels(diag), ""), ncol=1, lty=c(1,1), lwd=c(3,3), cex=0.7, col=c(colorcode_diag, "white"), bty = "n")
legend("bottom", legend = c("TILs & Grade +1/0", "TILs & Grade +3/+2", "", "", "ER & PgR Neg.", "ER & PgR Pos.", "", "", "Her2 Neg.", "Her2 +1", "Her2 +2.", "Her2 +3"), ncol=3, lty=c(1,1), lwd=c(3,3), cex=0.7, col=c("#0071AA", "#032A63", "white", "white", "seashell", "grey60", "white", "white", c("seashell", brewer.pal(4,"Greys")[2:4])), bty = "n")
dev.off()



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Make den dendrogram using the ward.D2 clustering metric
dend <- as.dendrogram(hclust(dist(t(batch_corr_diag)), method = "ward.D2"))

# Plot dendrogram
pdf("Dendogram.pdf", height = 10, width = 12)
par(mar=c(10,1,1,1))
lefcol <- c(rep("#0071AA", 14), rep("#032A63", 20))
dend %>%
  set("labels_col", value = brewer.pal(6, "Blues")[-c(1:3)], k=3) %>%
  set("branches_k_color", value =  brewer.pal(6, "Blues")[-c(1:3)], k = 3) %>%
  set("leaves_pch", 19)  %>% set("leaves_cex", 0.7) %>% set("leaves_col", lefcol) %>% plot()
# Add the colored bar
colored_bars(cbind(my.ST, my.ER, my.PGR, my.HER2, my.TILS, my.GR), dend, rowLabels = c("Subtype", "ER", "PgR", "Her2","TILs","Grade"))
dev.off()




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING DIAGNOSIS.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# CALCULATING DIFFERENCES BETWEEN LINEAR MODELS: LogFC > 1 and FDR < 0.5 OR LogFC < -1 and FDR < 0.5. 

# Diagnosis with pool
diag_design <- model.matrix(~0+diag+pool)
DA_diag1 <- DA_all_contrasts(tifdata_small, diag_design, diag, "diag", 1, 0.05)

# Diagnosis with pool and TIL status
diag_design <- model.matrix(~0+diag+pool+TILs)
DA_diag2 <- DA_all_contrasts(tifdata_small, diag_design, diag, "diag", 1, 0.05)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Writing out tables with DA proteins
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#write_out(rbind(DA_diag$`diagHER2-diagLumA`[[1]],DA_diag$`diagHER2-diagLumA`[[2]]), "HER2_LumA_DA")
#write_out(rbind(DA_diag$`diagHER2-diagTNBC`[[1]],DA_diag$`diagHER2-diagTNBC`[[2]]), "HER2_TNBC_DA")
#write_out(rbind(DA_diag$`diagLumA-diagTNBC`[[1]], DA_diag$`diagLumA-diagTNBC`[[2]]), "LumA_TNBC_DA")



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING HORMONE RECEPTORS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Estrogen
ER_design <- model.matrix(~0+ER+pool)
DA_ER <- DA_all_contrasts(tifdata_small, ER_design, ER, "ER", 1, 0.05)

# Progesterone
PGR_design <- model.matrix(~0+PGR+pool)
DA_PGR <- DA_all_contrasts(tifdata_small, PGR_design, PGR, "PGR", 1, 0.05)

# Androgen
AR_design <- model.matrix(~0+AR+pool)
DA_AR <- DA_all_contrasts(tifdata_small, AR_design, AR, "AR", 1, 0.05)

# HER2 enriched
HER2_design <- model.matrix(~0+HER2Sim+pool)
DA_HER2 <- DA_all_contrasts(tifdata_small, HER2_design, HER2Sim, "HER2Sim", 1, 0.05)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#write_out(rbind(DA_ER$`ERERm-ERERp`[[1]],DA_ER$`ERERm-ERERp`[[2]]), "ER_DA_corrected")
#write_out(rbind(DA_PGR$`PGRPGRm-PGRPGRp`[[1]],DA_PGR$`PGRPGRm-PGRPGRp`[[2]]), "PGR_DA_corrected")
#write_out(rbind(DA_HER2$`HER2H0-HER2H2`[[1]], DA_HER2$`HER2H0-HER2H2`[[2]], DA_HER2$`HER2H1-HER2H2`[[2]]), "H01_H2_DA_corrected")
#write_out(rbind(DA_HER2$`HER2H0-HER2H3`[[1]], DA_HER2$`HER2H0-HER2H3`[[2]]), "H01_H3_DA_corrected_pool")



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING TILS AND TUMOUR GRADE
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#  Make Design Matix
TILs_design <- model.matrix(~0+TILs+pool)
DA_TILs <- DA_all_contrasts(tifdata_small, TILs_design, TILs, "TILs", 1, 0.05)
DA_TILs_permissive <- DA_all_contrasts(tifdata_small, TILs_design, TILs, "TILs", 0.5, 0.05)


#  Make Design Matix
Gr_design <- model.matrix(~0+Gr+pool)
DA_Grp <- DA_all_contrasts(tifdata_small, Gr_design, Gr, "Gr", 0.5, 0.11)

Gr_design <- model.matrix(~0+Gr)
DA_Gr <- DA_all_contrasts(tifdata_small, Gr_design, Gr, "Gr", 0.8, 0.05)

GR_int <- c(intersect(rownames(DA_Gr$`Grhigh-Grlow`[[1]]), rownames(DA_Grp$`Grhigh-Grlow`[[1]])), intersect(rownames(DA_Gr$`Grhigh-Grlow`[[2]]), rownames(DA_Grp$`Grhigh-Grlow`[[2]])))

DA_Gr <- rbind(DA_Gr$`Grhigh-Grlow`[[1]][rownames(DA_Gr$`Grhigh-Grlow`[[1]]) %in% GR_int, ], DA_Gr$`Grhigh-Grlow`[[2]][rownames(DA_Gr$`Grhigh-Grlow`[[2]]) %in% GR_int, ])

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#write_out(rbind(DA_TILs$`TILshigh-TILslow`[[1]],DA_TILs$`TILshigh-TILslow`[[2]]), "TILsHigh_Low_pool_DA")
#write_out(rbind(DA_TILs_permissive$`TILshigh-TILslow`[[1]],DA_TILs_permissive$`TILshigh-TILslow`[[2]]), "TILsHigh_Low_permissive_pool_DA")




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LIMMA DIFFERENTIAL ABUNDANCE ANALYSIS USING CLUSTERS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Three Clusters
C3_design <- model.matrix(~0+C3+pool)
DA_C3 <- DA_all_contrasts(tifdata_small, C3_design, C3, "C3", 1, 0.05)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WRITING OUT TABLES OF DA PROTEINS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#write_out(rbind(DA_C3$`C3C2-C3C1`[[1]],DA_C3$`C3C2-C3C1`[[2]]), "C2_C1_pool_DA")
#write_out(rbind(DA_C3$`C3C3-C3C1`[[1]],DA_C3$`C3C3-C3C1`[[2]]), "C3_C1_pool_DA")
#write_out(rbind(DA_C3$`C3C3-C3C2`[[1]],DA_C3$`C3C3-C3C2`[[2]]), "C3_C2_pool_DA")



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Overlap DA
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GENE ONTOLOGY ENRICHMENT USIGN DIAGNOSIS.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

univ_tif <- rownames(tifdata_small)
univ_tif <- unique(sort(univ_tif))

H_L_GO <- TOPGO(ont = "BP", univ = univ_tif, intpro = unique(c(rownames(DA_diag1$`diagHER2-diagLumA`[[1]]), rownames(DA_diag1$`diagHER2-diagLumA`[[2]]))), GO_background = GO_background, 20)

H_L_GO_up <- GO_visual(rownames(DA_diag1$`diagHER2-diagLumA`[[1]]), univ_tif, GO_background, 30, "H_L_up", TRUE, TRUE)
H_L_GO_down <- GO_visual(rownames(DA_diag1$`diagHER2-diagLumA`[[2]]), univ_tif, GO_background, 30, "H_L_down", TRUE, TRUE)

H_T_GO <- TOPGO(ont = "BP", univ = univ_tif, intpro = unique(c(rownames(DA_diag1$`diagHER2-diagTNBC`[[1]]), rownames(DA_diag1$`diagHER2-diagTNBC`[[2]]))), GO_background = GO_background, 20)

H_T_GO_up <- GO_visual(rownames(DA_diag1$`diagHER2-diagTNBC`[[1]]), univ_tif, GO_background, 30, "H_T_up", TRUE, TRUE)
H_T_GO_down <- GO_visual(rownames(DA_diag1$`diagHER2-diagTNBC`[[2]]), univ_tif, GO_background, 30, "H_T_down", TRUE, TRUE)

L_T_GO <- TOPGO(ont = "BP", univ = univ_tif, intpro = unique(c(rownames(DA_diag1$`diagLumA-diagTNBC`[[1]]), rownames(DA_diag1$`diagLumA-diagTNBC`[[2]]))), GO_background = GO_background, 20)

L_T_GO_up <- GO_visual(rownames(DA_diag1$`diagLumA-diagTNBC`[[1]]), univ_tif, GO_background, 30, "L_T_up", TRUE, TRUE)
L_T_GO_down <- GO_visual(rownames(DA_diag1$`diagLumA-diagTNBC`[[2]]), univ_tif, GO_background, 30, "L_T_down", TRUE, TRUE)




intsec.list <- list(rownames(DA_diag1$`diagHER2-diagLumA`[[2]]), rownames(DA_diag1$`diagHER2-diagTNBC`[[2]]), rownames(DA_diag1$`diagLumA-diagTNBC`[[2]]), rownames(DA_diag1$`diagHER2-diagLumA`[[1]]),rownames(DA_diag1$`diagHER2-diagTNBC`[[1]]),rownames(DA_diag1$`diagLumA-diagTNBC`[[1]]))
names(intsec.list) <- c("HER2LumADown", "HER2TNBCDown", "LumATNBCDown", "HER2LumAUp", "HER2TNBCUp", "LumATNBCUp")

Ov_Subtypes <- plot_upsetR(intsec.list, c("HER2TNBCDown", "HER2LumADown"), "Overlap_Subtypes", TRUE, FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ODDS RATIO DA PROTEINS MICRO VESICLES VS. TIF
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

exo_tif_ov <- intersect(tif_full_name, exosomes)
exo_tif_no_ov <- tif_full_name[!tif_full_name %in% exo_tif_ov]

sig_d_genes <- unique(sort(uniprot_to_name(sig_d)$name)) 

get_stats(sig_d_genes, exo_tif_no_ov, exo_tif_ov)






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# LASSO REGRESSION
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subtype

LASSOdiag <- LASSO_Apply(tifdata_small, diag, FALSE, TRUE, NULL, NULL)
LASSOdiag_B1 <- LASSO_Apply(tifdata_small, diag, FALSE, TRUE, pool, NULL)
LASSOdiag_B2 <- LASSO_Apply(tifdata_small, diag, FALSE, TRUE, pool, TILs)

LASSOdiagAll <- list(LASSOdiag, LASSOdiag_B1, LASSOdiag_B2)
names(LASSOdiagAll) <- c("LASSOdiag", "LASSOdiag_B1", "LASSOdiag_B2")
save(LASSOdiagAll, file = "LASSOdiagAll.Rdata")


#LASSOdiag <- data.frame(Accession =  unique(c(LASSOdiag$Vars, LASSOdiag_B1$Vars, LASSOdiag_B2$Vars)))
#LASSOdiag <- merge(LASSOdiag, uniprot_to_name(LASSOdiag$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOdiag, "LASSO_diag_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Subtype without Her2

data_LumTN <- tifdata_small[,-c(2,5,28)]
group_LumTN <- as.factor(as.character(diag[-c(2,5,28)]))
pool_LumTN <- as.factor(as.character(pool[-c(2,5,28)]))
TILs_LumTN <- as.factor(as.character(TILs[-c(2,5,28)]))

LASSOLumTN <- LASSO_Apply(data_LumTN, group_LumTN, TRUE, FALSE, NULL, NULL)
LASSOLumTN_B1 <- LASSO_Apply(data_LumTN, group_LumTN, TRUE, FALSE, pool_LumTN, NULL)
LASSOLumTN_B2 <- LASSO_Apply(data_LumTN, group_LumTN, TRUE, FALSE, pool_LumTN, TILs_LumTN)

LASSOLumTNAll <- list(LASSOLumTN, LASSOLumTN_B1, LASSOLumTN_B2)
names(LASSOLumTNAll) <- c("LASSOLumTN", "LASSOLumTN_B1", "LASSOLumTN_B2")
save(LASSOLumTNAll, file = "LASSOLumTNAll.Rdata")


#LASSOLumTNAll <- data.frame(Accession = Reduce(intersect, list(LASSOLumTNAll$Vars, LASSOLumTNAll_B1$Vars, LASSOLumTNAll_B2$Vars)))
#LASSOLumTNAll <- merge(LASSOLumTNAll, uniprot_to_name(LASSOLumTNAll$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOLumTNAll, "LASSO_LumTNAll_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tumor Grade

LASSOGR <- LASSO_Apply(tifdata_small, Gr, FALSE, FALSE, NULL, NULL)
LASSOGR_B1 <- LASSO_Apply(tifdata_small, Gr, FALSE, FALSE, pool, NULL)
LASSOGR_B2 <- LASSO_Apply(tifdata_small, Gr, FALSE, FALSE, pool, TILs)

LASSOGRAll <- list(LASSOGR, LASSOGR_B1, LASSOGR_B2)
names(LASSOGRAll) <- c("LASSOGR", "LASSOGR_B1", "LASSOGR_B2")
save(LASSOGRAll, file = "LASSOGRAll.Rdata")


#LASSOGR <- data.frame(Accession = Reduce(intersect, list(LASSOGR$Vars, LASSOGR_B1$Vars, LASSOGR_B2$Vars)))
#LASSOGR <- merge(LASSOGR, uniprot_to_name(LASSOGR$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOGR, "LASSO_GR_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TILs

LASSOTILs <- LASSO_Apply(tifdata_small, TILs, FALSE, FALSE, NULL, NULL)
LASSOTILs_B1 <- LASSO_Apply(tifdata_small, TILs, FALSE, FALSE, pool, NULL)

LASSOTILsAll <- list(LASSOTILs, LASSOTILs_B1)
names(LASSOTILsAll) <- c("LASSOTILs", "LASSOTILs_B1")
save(LASSOTILsAll, file = "LASSOTILsAll.Rdata")


#LASSOTILs <- data.frame(Accession = Reduce(intersect, list(LASSOTILs$Vars, LASSOTILs_B1$Vars)))
#LASSOTILs <- merge(LASSOTILs, uniprot_to_name(LASSOTILs$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOTILs, "LASSO_TILs_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Her2 status

LASSOHer2 <- LASSO_Apply(tifdata_small, HER2Sim, FALSE, FALSE, NULL, NULL)
LASSOHer2_B1 <- LASSO_Apply(tifdata_small, HER2Sim, FALSE, FALSE, pool, NULL)
LASSOHer2_B2 <- LASSO_Apply(tifdata_small, HER2Sim, FALSE, FALSE, pool, TILs)

LASSOHer2All <- list(LASSOHer2, LASSOHer2_B1, LASSOHer2_B2)
names(LASSOHer2All) <- c("LASSOHer2", "LASSOHer2_B1", "LASSOHer2_B2")
save(LASSOHer2All, file = "LASSOHer2All.Rdata")


#LASSOHer2 <- data.frame(Accession = Reduce(intersect, list(LASSOHer2$Vars, LASSOHer2_B1$Vars, LASSOHer2_B2$Vars)))
#LASSOHer2 <- merge(LASSOHer2, uniprot_to_name(LASSOHer2$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOHer2, "LASSO_Her2_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Estrogen

# Split TRUE
LASSOERval <- LASSO_Apply(tifdata_small, ER, TRUE, FALSE, NULL, NULL)
LASSOERval_B1  <- LASSO_Apply(tifdata_small, ER, TRUE, FALSE, pool, NULL)
LASSOERval_B2  <- LASSO_Apply(tifdata_small, ER, TRUE, FALSE, pool, TILs)

LASSOERAllval <- list(LASSOERval, LASSOERval_B1, LASSOERval_B2)
names(LASSOERAllval) <- c("LASSOER", "LASSOER_B1", "LASSOER_B2")
save(LASSOERAllval, file = "LASSOERAllval.Rdata")


# Split FALSE
LASSOER <- LASSO_Apply(tifdata_small, ER, FALSE, FALSE, NULL, NULL)
LASSOER_B1  <- LASSO_Apply(tifdata_small, ER, FALSE, FALSE, pool, NULL)
LASSOER_B2  <- LASSO_Apply(tifdata_small, ER, FALSE, FALSE, pool, TILs)

LASSOERAll <- list(LASSOER, LASSOER_B1, LASSOER_B2)
names(LASSOERAll) <- c("LASSOER", "LASSOER_B1", "LASSOER_B2")
save(LASSOERAll, file = "LASSOERAll.Rdata")

#LASSOER <- data.frame(Accession = Reduce(intersect, list(LASSOER$Vars, LASSOER_B1$Vars, LASSOER_B2$Vars)))
#LASSOER <- merge(LASSOER, uniprot_to_name(LASSOER$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOER, "LASSO_ER_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Progesterone 

# Split TRUE
LASSOPgRval <- LASSO_Apply(tifdata_small, PGR, TRUE, FALSE, NULL, NULL)
LASSOPgRval_B1  <- LASSO_Apply(tifdata_small, PGR, TRUE, FALSE, pool, NULL)
LASSOPgRval_B2  <- LASSO_Apply(tifdata_small, PGR, TRUE, FALSE, pool, TILs)

LASSOPgRAllval <- list(LASSOPgRval, LASSOPgRval_B1, LASSOPgRval_B2)
names(LASSOPgRAllval) <- c("LASSOPgR", "LASSOPgR_B1", "LASSOPgR_B2")
save(LASSOPgRAllval, file = "LASSOPgRAllval.Rdata")

# Split FALSE
LASSOPgR <- LASSO_Apply(tifdata_small, PGR, FALSE, FALSE, NULL, NULL)
LASSOPgR_B1  <- LASSO_Apply(tifdata_small, PGR, FALSE, FALSE, pool, NULL)
LASSOPgR_B2  <- LASSO_Apply(tifdata_small, PGR, FALSE, FALSE, pool, TILs)

LASSOPgRAll <- list(LASSOPgR, LASSOPgR_B1, LASSOPgR_B2)
names(LASSOPgRAll) <- c("LASSOPgR", "LASSOPgR_B1", "LASSOPgR_B2")
save(LASSOPgRAll, file = "LASSOPgRAll.Rdata")


#LASSOPgR <- data.frame(Accession = Reduce(intersect, list(LASSOPgR$Vars, LASSOPgR_B1$Vars, LASSOPgR_B2$Vars)))
#LASSOPgR <- merge(LASSOPgR, uniprot_to_name(LASSOPgR$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(LASSOPgR, "LASSO_PgR_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




# Average Weights LASSO
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightSub1 <- AverageWeight(LASSOdiagAll)
WeightSub2 <- AverageWeight(LASSOLumTNAll)
WeightSub <- AverageWeight(list(WeightSub1[,-5], WeightSub2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightER1 <- AverageWeight(LASSOERAll)
WeightER2 <- AverageWeight(LASSOERAllval)
WeightER <- AverageWeight(list(WeightER1[,-5], WeightER2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightPgR1 <- AverageWeight(LASSOPgRAll)
WeightPgR2 <- AverageWeight(LASSOPgRAllval)
WeightPgR <- AverageWeight(list(WeightPgR1[,-5], WeightPgR2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightHer2 <- AverageWeight(LASSOHer2All)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


WeightTILs <- AverageWeight(LASSOTILsAll)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


WeightGR <- AverageWeight(LASSOGRAll)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# RANDOM FOREST
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# Subtype
RFdiag <- RF_Apply(tifdata_small, diag, FALSE, NULL, NULL)
RFdiag_B1 <- RF_Apply(tifdata_small, diag, FALSE, pool, NULL)
RFdiag_B2 <- RF_Apply(tifdata_small, diag, FALSE, pool, TILs)

RFdiagAll <- list(RFdiag, RFdiag_B1, RFdiag_B2)
names(RFdiagAll) <- c("RFdiag", "RFdiag_B1", "RFdiag_B2")
save(RFdiagAll, file = "RFdiagAll.Rdata")

#RFdiag <- data.frame(Accession =  unique(c(RFdiag$Vars, RFdiag_B1$Vars, RFdiag_B2$Vars)))
#RFdiag <- merge(RFdiag, uniprot_to_name(RFdiag$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFdiag, "RF_diag_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Subtype without Her2

data_LumTN <- tifdata_small[,-c(2,5,28)]
group_LumTN <- as.factor(as.character(diag[-c(2,5,28)]))
pool_LumTN <- as.factor(as.character(pool[-c(2,5,28)]))
TILs_LumTN <- as.factor(as.character(TILs[-c(2,5,28)]))

RFLumTN <- RF_Apply(data_LumTN, group_LumTN, TRUE, NULL, NULL)
RFLumTN_B1 <- RF_Apply(data_LumTN, group_LumTN, TRUE, pool_LumTN, NULL)
RFLumTN_B2 <- RF_Apply(data_LumTN, group_LumTN, TRUE, pool_LumTN, TILs_LumTN)

RFLumTNAll <- list(RFLumTN, RFLumTN_B1, RFLumTN_B2)
names(RFLumTNAll) <- c("RFLumTN", "RFLumTN_B1", "RFLumTN_B2")
save(RFLumTNAll, file = "RFLumTNAll.Rdata")

#RFLumTN <- data.frame(Accession = Reduce(intersect, list(RFLumTN$Vars, RFLumTN_B1$Vars, RFLumTN_B2$Vars)))
#RFLumTN <- merge(RFLumTN, uniprot_to_name(RFLumTN$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFLumTN, "RF_LumTN_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tumor Grade

RFGR <- RF_Apply(tifdata_small, Gr, FALSE, NULL, NULL)
RFGR_B1 <- RF_Apply(tifdata_small, Gr, FALSE, pool, NULL)
RFGR_B2 <- RF_Apply(tifdata_small, Gr, FALSE, pool, TILs)

RFGRAll <- list(RFGR, RFGR_B1, RFGR_B2)
names(RFGRAll) <- c("RFGR", "RFGR_B1", "RFGR_B2")
save(RFGRAll, file = "RFGRAll.Rdata")


#RFGR <- data.frame(Accession = Reduce(intersect, list(RFGR$Vars, RFGR_B1$Vars, RFGR_B2$Vars)))
#RFGR <- merge(RFGR, uniprot_to_name(RFGR$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFGR, "RF_GR_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tumor infiltrating-lymphocytes

RFTILs <- RF_Apply(tifdata_small, TILs, FALSE, NULL, NULL)
RFTILs_B1 <- RF_Apply(tifdata_small, TILs, FALSE, pool, NULL)

RFTILsAll <- list(RFTILs, RFTILs_B1)
names(RFTILsAll) <- c("RFTILs", "RFTILs_B1")
save(RFTILsAll, file = "RFTILsAll.Rdata")

#RFTILs <- data.frame(Accession = Reduce(intersect, list(RFTILs$Vars, RFTILs_B1$Vars)))
#RFTILs <- merge(RFTILs, uniprot_to_name(RFTILs$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFTILs, "RF_TILs_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Her2 status

RFHer2 <- RF_Apply(tifdata_small, HER2Sim, FALSE, NULL, NULL)
RFHer2_B1 <- RF_Apply(tifdata_small, HER2Sim, FALSE, pool, NULL)
RFHer2_B2 <- RF_Apply(tifdata_small, HER2Sim, FALSE, pool, TILs)

RFHer2All <- list(RFHer2, RFHer2_B1, RFHer2_B2)
names(RFHer2All) <- c("RFHer2", "RFHer2_B1", "RFHer2_B2")
save(RFHer2All, file = "RFHer2All.Rdata")


#RFHer2 <- data.frame(Accession = Reduce(intersect, list(RFHer2$Vars, RFHer2_B1$Vars, RFHer2_B2$Vars)))
#RFHer2 <- merge(RFHer2, uniprot_to_name(RFHer2$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFHer2, "RF_Her2_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Estrogen

# Split TRUE
RFER <- RF_Apply(tifdata_small, ER, TRUE, NULL, NULL)
RFER_B1  <- RF_Apply(tifdata_small, ER, TRUE, pool, NULL)
RFER_B2  <- RF_Apply(tifdata_small, ER, TRUE, pool, TILs)

RFERAllval <- list(RFER, RFER_B1, RFER_B2)
names(RFERAllval) <- c("RFER", "RFER_B1", "RFER_B2")
save(RFERAllval, file = "RFERAllval.Rdata")

# Split FALSE
RFER <- RF_Apply(tifdata_small, ER, FALSE, NULL, NULL)
RFER_B1  <- RF_Apply(tifdata_small, ER, FALSE, pool, NULL)
RFER_B2  <- RF_Apply(tifdata_small, ER, FALSE, pool, TILs)

RFERAll <- list(RFER, RFER_B1, RFER_B2)
names(RFERAll) <- c("RFER", "RFER_B1", "RFER_B2")
save(RFERAll, file = "RFERAll.Rdata")

#RFER <- data.frame(Accession = Reduce(intersect, list(RFER$Vars, RFER_B1$Vars, RFER_B2$Vars)))
#RFER <- merge(RFER, uniprot_to_name(RFER$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFER, "RF_ER_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Progesterone

# Split TRUE
RFPgR <- RF_Apply(tifdata_small, PGR, TRUE, NULL, NULL)
RFPgR_B1  <- RF_Apply(tifdata_small, PGR, TRUE, pool, NULL)
RFPgR_B2  <- RF_Apply(tifdata_small, PGR, TRUE, pool, TILs)

RFPgRAllval <- list(RFPgR, RFPgR_B1, RFPgR_B2)
names(RFPgRAllval) <- c("RFPgR", "RFPgR_B1", "RFPgR_B2")
save(RFPgRAllval, file = "RFPgRAllval.Rdata")

# Split FALSE
RFPgR <- RF_Apply(tifdata_small, PGR, FALSE, NULL, NULL)
RFPgR_B1  <- RF_Apply(tifdata_small, PGR, FALSE, pool, NULL)
RFPgR_B2  <- RF_Apply(tifdata_small, PGR, FALSE, pool, TILs)

RFPgRAll <- list(RFPgR, RFPgR_B1, RFPgR_B2)
names(RFPgRAll) <- c("RFPgR", "RFPgR_B1", "RFPgR_B2")
save(RFPgRAll, file = "RFPgRAll.Rdata")


#RFPgR <- data.frame(Accession = Reduce(intersect, list(RFPgR$Vars, RFPgR_B1$Vars, RFPgR_B2$Vars)))
#RFPgR <- merge(RFPgR, uniprot_to_name(RFPgR$Accession), by = "Accession", all.x = TRUE, all.y = FALSE)
#write.table(RFPgR, "RF_PgR_corrected.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



# Average Weights
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightSub1 <- AverageWeight(RFdiagAll)
WeightSub2 <- AverageWeight(RFLumTNAll)
WeightSub <- AverageWeight(list(WeightSub1[,-5], WeightSub2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


WeightER1 <- AverageWeight(RFERAll)
WeightER2 <- AverageWeight(RFERAllval)
WeightER <- AverageWeight(list(WeightER1[,-5], WeightER2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


WeightPgR1 <- AverageWeight(RFPgRAll)
WeightPgR2 <- AverageWeight(RFPgRAllval)
WeightPgR <- AverageWeight(list(WeightPgR1[,-5], WeightPgR2[,-5]))

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightHer2 <- AverageWeight(RFHer2All)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightTILs <- AverageWeight(RFTILsAll)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

WeightGR <- AverageWeight(RFGRAll)

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVerlap DA proteins
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
setwd(paste0(my.wd, "/Results/Tables/DA/Subtypes/Corrected_for_Pool_permissive"))

HER2_TNBC_DA <- read.delim("HER2_TNBC_DA_corrected_pool.txt", header = TRUE)
HER2_LumA_DA <- read.delim("HER2_LumA_DA_corrected_pool.txt", header = TRUE)
LumA_TNBC_DA <- read.delim("LumA_TNBC_DA_corrected_pool.txt", header = TRUE)

intsec.list <- list(unique(as.character(HER2_TNBC_DA[HER2_TNBC_DA$dir == "up",]$Accession)),
                    unique(as.character(HER2_LumA_DA[HER2_LumA_DA$dir == "up",]$Accession)), 
                    unique(as.character(LumA_TNBC_DA[LumA_TNBC_DA$dir == "up",]$Accession)), 
                    unique(as.character(HER2_LumA_DA[HER2_LumA_DA$dir == "down",]$Accession)),
                    unique(as.character(HER2_TNBC_DA[HER2_TNBC_DA$dir == "down",]$Accession)),
                    unique(as.character(LumA_TNBC_DA[LumA_TNBC_DA$dir == "down",]$Accession)))

names(intsec.list) <- c("Up in Her2 vs TNBC", "Up in Her2 vs Lum", "Up in Lum vs TNBC", "Up in Lum vs Her2", "Up in TNBC vs Her2", "Up in TNBC vs Lum")


coloST <- c("#FAA275", "#FAA275", "#FFFD82", "#FFFD82", "#BF5454", "#BF5454")
DAplot <- plot_upsetR(intsec.list, names(intsec.list), "Figure4", coloST, TRUE, FALSE)






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Overlap results DA, LASSO and RF
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd("/Results/Tables/DA/Subtypes/Corrected_for_Pool_permissive/")
HER2_LumA_DA <- read.table("HER2_LumA_DA_corrected_pool.txt", header=TRUE)
HER2_TNBC_DA <- read.table("HER2_TNBC_DA_corrected_pool.txt", header = TRUE)
LumA_TNBC_DA <- read.table("LumA_TNBC_DA_corrected_pool.txt", header = TRUE)


setwd("/Results/Tables/LASSO")
Subtypes_LASSO <- read.table("LASSO_subtypes_corrected.txt", header=TRUE)


setwd("/Results/Tables/RF")
Subtypes_RF <- read.table("RF_subtypes_corrected.txt", header=TRUE)


my.white <- c("white", "white", "white")
my.DLRcol <- c("#08605F", "#BFBDC1",  "#545E75")

my.DLRcol <-c("#6279B8", "#CDE7B0", "#34435E")

intsec.list <- list(as.character(HER2_LumA_DA$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "LASSO", "RF")
Ov_HER2_LumA_DA_LA_RF <- plot_upsetR(intsec.list,  names(intsec.list), "Best_cand_HER2_LumA", colov[c(1,2,3)], FALSE, FALSE)

pdf("Ov_HER2_LumA_DA_LA_RF.pdf", height = 11, width = 14)
Overlap_GOplot(Ov_HER2_LumA_DA_LA_RF, my.white, my.DLRcol)
dev.off()


intsec.list <- list(as.character(HER2_TNBC_DA$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "LASSO", "RF")
Ov_HER2_TNBC_DA_LA_RF <- plot_upsetR(intsec.list,  names(intsec.list), "Best_cand_HER2_TNBC", colov[c(1,2,3)], FALSE, FALSE)

pdf("Ov_HER2_TNBC_DA_LA_RF.pdf", height = 9, width = 10)
Overlap_GOplot(Ov_HER2_TNBC_DA_LA_RF, my.white, my.DLRcol)
dev.off()

intsec.list <- list(as.character(LumA_TNBC_DA$Accession), as.character(Subtypes_LASSO$Accession), as.character(Subtypes_RF$Accession))
names(intsec.list) <- c("DA", "LASSO", "RF")
Ov_LumA_TNBC_DA_LA_RF <- plot_upsetR(intsec.list, names(intsec.list), "Best_cand_LumA_TNBC", colov[c(1,2,3)], FALSE, FALSE)

pdf("Ov_LumA_TNBC_DA_LA_RF.pdf", height = 11, width = 14)
Overlap_GOplot(Ov_LumA_TNBC_DA_LA_RF, my.white, my.DLRcol)
dev.off()





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PATHWAY ANALYSIS USING DIAGNOSIS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
univ_tif_enz <- unique(sort(as.character(getBM(attributes=c("uniprotswissprot", "entrezgene_id"), filters="uniprotswissprot", values=tif_full_id, mart=ensembl)$entrezgene_id)))


H_L_up_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagHER2-diagLumA`[[1]])), "H_L_up_pw", TRUE)
H_L_down_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagHER2-diagLumA`[[2]])), "H_L_down_pw", TRUE)

H_T_up_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagHER2-diagTNBC`[[1]])), "H_T_up_pw", TRUE)
H_T_down_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagHER2-diagTNBC`[[2]])), "H_T_down_pw", TRUE)

L_T_up_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagLumA-diagTNBC`[[1]])), "L_H_up_pw", TRUE)
L_T_down_pw <- enrich_pathway(univ_tif_enz, as.character(rownames(DA_diag$`diagLumA-diagTNBC`[[2]])), "L_H_down_pw", TRUE)



entrez.data <- u2e.map[u2e.map$uniprot %in% as.character(rownames(DA_diag$`diagHER2-diagLumA`[[1]])), ]
my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.05, as.character(univ_tif_enz$entrez), minGSSize = 4, readable = T)




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DE PROTEINS OVERLAP WITH DE CYTOKINES
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cytokines <- read.delim("Cytokines.txt", header=TRUE)
cytokines_accession <- sort(as.character(cytokines$UniprotID))



ov_d_cyto <- sig_d[sig_d %in% cytokines_accession]
#ov_pca_d_cyto <- PCA_diag[PCA_diag %in% cytokines_accession2]




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERLAP WITH PAM50 GENES
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# OVERLAP WITH OTHER DATASETS SUBTYPES
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------








# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CORRELATIONS - DIAG
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# IHC and MS plots
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

setwd(paste0(my.wd, "Results/Tables/IHC/"))

# Load Immunohistochemistry data
IHC <- read.delim("IHCscoresFull.txt", header = TRUE)
proteins <- rownames(IHC)
IHC$protein <- proteins
IHC <- melt(IHC)
IHC$protein <- factor(IHC$protein, levels=proteins)

  
# Heatmap
ggplot(IHC, aes(protein, variable, fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Load Currated and transformed Mass-spec data
MS <- read.delim("MSscoresFull.txt", header = TRUE)
MS$protein <- proteins
MS <- melt(MS)
MS$protein <- factor(MS$protein, levels=proteins)

# Heatmap
ggplot(MS, aes(protein, variable, fill = value)) + geom_tile() + scale_fill_viridis(begin =0, end =1, option = "viridis", direction = -1) + theme_bw() + coord_equal() + theme(axis.text.x = element_text(angle = 90, hjust = 1))

# -----------------------------------------------------------------------------
# Trend Plots for proteins

MSIHC <- cbind(MS, IHC[,3])
colnames(MSIHC) <- c("protein", "STP", "MS", "IHC")
MSIHC$STP <- as.factor(as.character(gsub(pattern = ".[^.]*$", "", MSIHC$STP)))
MSIHC <- MSIHC[-which(is.na(MSIHC$IHC)),]
MSIHC$IHC <- as.factor(as.character(MSIHC$IHC))


# Plot (Box + scatter + trend)
ggplot(MSIHC, aes(x=IHC, y=MS)) + geom_boxplot(alpha = 0.7) + 
  geom_point(aes(color=STP)) + 
  geom_smooth(aes(group = STP, color=STP), method = "lm", se=FALSE, size =0.8) + 
  scale_color_manual(values = c("#FAA275","#DDDC71", "#BF5454")) +
  theme_light() +
  theme(axis.text.x = element_text(color = "grey40", size = 10, face = "plain"),
        axis.text.y = element_text(color = "grey40", size = 10, face = "plain"),
        axis.title.x = element_text(color = "grey40", size = 11, face = "plain"),
        axis.title.y = element_text(color = "grey40", size = 11, face = "plain"),
        legend.title=element_text(size=12), legend.text=element_text(size=11), strip.background = element_rect(fill = "white", linetype = "solid", color="black"), strip.text.x = element_text(color = "black")) +
  facet_wrap(~ protein, ncol = 5) + theme(strip.text.x = element_text(size = 11, face = "bold"))

# -----------------------------------------------------------------------------

my.cols <- c("#FCB97D", "#FBB7C0", "#87BCDE", "#70A9A1", "#2B3A67", "#3A5743", "#DB4552")

MSIHC$IHC <- as.numeric(as.character(MSIHC$IHC))
MSIHC$IHC <- ifelse(MSIHC$IHC == 0.5, 0.0, MSIHC$IHC)
MSIHC$IHC <- ifelse(MSIHC$IHC == 1.5, 1.0, MSIHC$IHC)
MSIHC$IHC <- ifelse(MSIHC$IHC == 2.5, 2.0, MSIHC$IHC)
MSIHC$IHC <- as.factor(MSIHC$IHC)

my.cols <- c("#FCB97D", "#FBB7C0", "#87BCDE", "#70A9A1", "#2B3A67", "#3A5743", "#DB4552")
my.cols <- c("#FCB97D","#70A9A1","#87BCDE", "#2B3A67")

# Logistisk Regression Plots
s.MSIHC <- split(MSIHC, f = MSIHC$protein)

polr.MSIHC <- lapply(s.MSIHC, function(x) polr(IHC ~ MS, data = x, Hess=TRUE))
#confint.MSIHC <- lapply(polr.MSIHC, function(x) confint(x))
test.MSIHC <- lapply(polr.MSIHC, function(x) OLR_pval(x))
new.MSIHC <- lapply(s.MSIHC, function(x) data.frame(MS = rep(seq(from = min(x$MS), to = max(x$MS), length.out = 100), 4)))

n.MSIHC <- names(polr.MSIHC)
prob.MSIHC <- list()

for (idx in 1:length(polr.MSIHC)) {
  m.prob <- cbind(new.MSIHC[[idx]], predict(polr.MSIHC[[idx]], new.MSIHC[[idx]], type = "probs"))
  m.prob$Protein <- rep(n.MSIHC[[idx]], nrow(m.prob))
  prob.MSIHC[[idx]] <- m.prob
}

names(prob.MSIHC) <- n.MSIHC


lp.MSIHC <- lapply(prob.MSIHC, function(x) melt(x, id.vars = c("MS", "Protein"), variable.name = "Level", value.name="Probability"))
plot.MSIHC <- do.call(rbind, lp.MSIHC)

remove <- sort(c(which(plot.MSIHC$Protein == "TMEM51" & plot.MSIHC$MS < -1.3), which(plot.MSIHC$Protein == "MIEN1" & plot.MSIHC$MS > 2), which(plot.MSIHC$Protein == "PIP4K2B" & plot.MSIHC$MS > 1.3)))
plot.MSIHC <- plot.MSIHC[-remove,]

plot.MSIHC$Protein<- as.factor(plot.MSIHC$Protein)
plot.MSIHC$variable <- as.factor(as.character(plot.MSIHC$variable))


ggplot(plot.MSIHC, aes(x = MS, y = value, colour = variable)) + 
  geom_line(size = 1) + theme_light() + scale_colour_manual(values = my.cols) +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(size = 10, face = "plain"),
        axis.text.y = element_text(size = 10, face = "plain"),
        axis.title.x = element_text(size = 11, face = "plain"),
        axis.title.y = element_text(size = 11, face = "plain"),
        legend.title=element_text(size=12), legend.text=element_text(size=11),
        strip.background = element_rect(fill = "white", linetype = "solid", color="black"),
        strip.text.x = element_text(color = "black")) +
  facet_wrap(~Protein, scales = "free", ncol = 5) + 
  theme(strip.text.x = element_text(size = 11, face = "bold"))



# -----------------------------------------------------------------------------



# Regression by Breast Cancer Subtype
s.MSIHC <- split(MSIHC, f = MSIHC$STP)

# Linear Regression LumA
LumA <- s.MSIHC$LumA
s.LumA <- split(LumA, f = LumA$protein)
l.LumA <- unlist(lapply(s.LumA, function(x) summary(lm(IHC ~ MS, data = x))$coefficients[2,4]))

# Linear Regression TNBC
TNBC <- s.MSIHC$TNBC
s.TNBC <- split(TNBC, f = TNBC$protein)
l.TNBC <- unlist(lapply(s.TNBC, function(x) summary(lm(IHC ~ MS, data = x))$coefficients[2,4]))

# Logistic Regression LumA
s.LumA$NAT1 <- NULL
polr.LumA <- lapply(s.LumA, function(x) polr(IHC ~ MS, data = x, Hess=TRUE))
test.LumA <- lapply(polr.LumA, function(x) OLR_pval(x))

# Logistic Regression TNBC
s.TNBC$TMEM51 <- NULL
polr.TNBC <- lapply(s.TNBC, function(x) polr(IHC ~ MS, data = x, Hess=TRUE))
test.TNBC <- lapply(polr.TNBC, function(x) OLR_pval(x))

df.pval <- data.frame("All" = p.adjust(l.MSIHC, "BH"), "LumA" = p.adjust(l.LumA, "BH"), "TNBC" = p.adjust(l.TNBC, "BH"))


