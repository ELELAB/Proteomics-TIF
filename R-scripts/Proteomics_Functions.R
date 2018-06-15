# LOAD PACKAGES
library(pamr)
library(limma)
library(sva)
library(openxlsx)
library(ggplot2)
library(dendextend)
library(heatmap.plus)
library(reshape)
library(VennDiagram)
library(gdata) 
library(topGO)
library(biomaRt)
library(GOSim)
library(corrplot)
library(pvclust)
library(stringr)
library(plyr)
library(gplots)
library(ReactomePA)
library(glmnet)
library(fitdistrplus)
library(caTools)
#library(rgl)
#library(plot3D)
library(randomForest)
library(glmnet)
library(e1071)
library(caret)
library(varSelRF)
library(UpSetR)
library(RColorBrewer)
#library(gage)
#library(pathview)




# ------------------------------------------------------------------------------------------
# GO TERMS ENRICHMENT
# ------------------------------------------------------------------------------------------

# Function which generates GO2Gene object. The function takes the argument: humandata =  HUGO human gene/protein database. 

read_in_data <- function(humandata) {
  hugo <- read.delim(humandata, header=F)
  colnames(hugo) <- c("database", "Accession", "name", "ID1", "GO_ID", "GO_ID_ref", "ID2", "Pubmed", "ID3", "Description", "full", "type", "taxon", "ID4", "where", "ID5", "ID6") # column names
  # create go2gene object
  goids <- as.character(unique(hugo$GO_ID))
  grps <- lapply(goids, FUN=function(x) as.character(hugo$Accession)[grep(x, as.character(hugo$GO_ID))])
  names(grps) <- goids
  return(grps)
}



# --------------------------------------------------------------------------------------

# Functions which creates GOobject. 
# The function takes arguments: ont = the ontology to enrich for, univ = a protein univers, intprot = a set of proteins of interest and
# GO_background = GO-background (dataframe where proteins are coupled to GO-terms)

GOobject <-function(ont, univ, intpro, GO_background){
  if(is.null(univ)){
    univ=(unique(unlist(GO_background)))
  }
  # factorise genelist
  geneList = factor(as.integer(univ%in%intpro))
  names(geneList) = univ
  #new topgo object
  GOdata = new("topGOdata", ontology = as.character(ont), allGenes = geneList, annot = annFUN.GO2genes, GO2genes=GO_background)
  return(GOdata)
}  

# --------------------------------------------------------------------------------------


# Function which creates GOobject and identifies significant GO-terms. 
# The function takes the arguments: ont = the ontology to enrich for, univ = a protein univers, intprot = a set of proteins of interest,
# GO_background = GO-background (list object where proteins are coupled to GO-terms) and nterms = number of GO hits to display.

TOPGO <- function(ont, univ, intpro, GO_background, nterms) {
  
  GOdata <- GOobject(ont, univ, intpro, GO_background)
  test.stat = new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
  resultFisher = getSigGroups(GOdata, test.stat)
  test.stat = new("elimCount", testStatistic = GOFisherTest, name = "Fisher test", cutOff = 0.01)
  resultElim = getSigGroups(GOdata, test.stat)
  
  res <- GenTable(GOdata, fisher = resultFisher, elim = resultElim, orderBy = "elim", ranksOf = "elim", topNodes = nterms)
  
  # pvalue correction
  FDR_fisher <- data.frame(p.adjust(res$fisher, method = "fdr", n = nrow(res)))
  colnames(FDR_fisher) <- "FDR_fisher"
  FDR_elim <- data.frame(p.adjust(res$elim, method = "fdr", n = nrow(res)))
  colnames(FDR_elim) <- "FDR_elim"
  
  res=cbind(res, FDR_fisher, FDR_elim)
  return(res)
}


# --------------------------------------------------------------------------------------


# Function that creates barplot. Input is results of TOPGO

barplot_func <- function(res) {
  temp <- data.frame(res$Term, log2(1/as.numeric(res$fisher)))
  colnames(temp) <- c("Term", "Inverse.Rank")
  temp$Term <- factor(temp$Term, levels=as.vector(temp$Term))
  ggplot(data=temp, aes(x=Term, y=Inverse.Rank)) + geom_bar(stat="identity", fill="#00BCD8") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
}


# --------------------------------------------------------------------------------------

# Function that creates correlation plot. Input is results of TOPGO

corrplot_func <- function(res) {
  term_corr <- getTermSim(res$GO.ID)
  rownames(term_corr) <- res$Term
  colnames(term_corr) <- NULL
  term_corr[which(term_corr==Inf)] <- 0 
  corrplot(term_corr, type = "lower", is.corr = FALSE, order = "hclust", tl.cex=1, tl.col = "black")
}





# ------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANLYSIS
# ------------------------------------------------------------------------------------------


# Batch correction of data. 
# The function takes arguments: my.data = abundance/expression data, my.group = patient groups,
# my.pool = vector with batches pool(e.g batch), my.tp = tumor percentage of sample and my.tils = TILs (level of infiltrating lymphocytes in sample).

combat_corrections <- function(my.data, my.group, my.pool, my.tils=NULL){
  mod_design <- model.matrix(~as.factor(my.group))
  batch_corr <- ComBat(my.data, as.factor(my.pool), mod_design, par.prior=TRUE,prior.plots=FALSE)
  if (!is.null(my.tils)) {
    batch_corr <- ComBat(batch_corr, as.factor(my.tils), mod_design, par.prior=TRUE,prior.plots=FALSE)
  }
  return(batch_corr)
}



# ----------------------------------------------------------------------------------------


# Multidimensional scaling plot for plotting data. 
# The function takes arguments: my.data = abundance/expression data, my.group = vector of patient groups for coloring,
# my.labels = vector of patient groups for labeling and my.colors = a vector for choice of colors (must be same length as levels of my.groups).



myMDSplot <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}





# ----------------------------------------------------------------------------------------
# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT GLYCANS:
# Takes as arguments;
# a contrast between groups of interest
# a dataframe, a design matrix with all comparisons
# cutoffs for logFC and FDR
# if blocking than a vector of patient IDs
# ----------------------------------------------------------------------------------------


DA_proteins <- function(my.contrast, my.data, my.design, coLFC, coFDR, my.block=NULL) {
  if(is.null(my.block)) {
    fit3 <- eBayes(contrasts.fit(lmFit(object = my.data, design = my.design), my.contrast))
  }
  else {
    corfit <- duplicateCorrelation(my.data, my.design, block=my.block) 
    fit3 <- eBayes(contrasts.fit(lmFit(object = my.data, design = my.design, block = my.block, correlation=corfit$consensus), my.contrast))
  }
  tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(my.data))
  
  up <- tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]
  down <- tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]
  
  up$dir <- rep("up", nrow(up))
  down$dir <- rep("down", nrow(down))
  
  final <- list(up, down)
  return(final)
}




# ----------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE:
# Takes as arguments;
# all contrasts between groups of interest
# a dataframe
# a design matrix with all comparisons
# cutoffs for logFC and FDR
# if blocking than a vector of patient IDs
# TRUE/FALSE statment specifying output format, if TRUE the function return a vector of glycan IDs only
# ----------------------------------------------------------------------------------------


DA_proteins_apply <- function(my.contrasts, my.data, my.design, coLFC, coFDR, my.vector, my.block=NULL) {
  my.proteins.l <- apply(my.contrasts, 2, function(x) DA_proteins(x, my.data, my.design, coLFC, coFDR, my.block)) 
  if(my.vector == TRUE) {
    my.proteins <- do.call(rbind, lapply(my.proteins.l, function(x) do.call(rbind, x)))
    my.proteins <- unique(do.call(rbind, strsplit(rownames(my.proteins), "[.]"))[,2])
    return(my.proteins)
  }
  else {
    return(my.proteins.l)
  }
}



# ----------------------------------------------------------------------------------------
# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS THE "DA_glycan_apply" FROM ABOVE.
# Takes as arguments;
# a dataframe
# a vector of groups do perform contrasts on (same length as ncol(dataframe))
# a cutoff for logFC and FDR
# if remove is different from NULL, a vector of indices to remove must be supplied
# ----------------------------------------------------------------------------------------



DA_all_contrasts <- function(my.data, my.design, my.group, my.group.name, my.logFC, my.FDR, my.block=NULL, my.remove=NULL) {
  if (!is.null(my.remove)) {
    my.data <- my.data[, -my.remove]
    my.group <- my.group[-my.remove]
  }
  combinations<- data.frame(t(combn(paste0(my.group.name, levels(my.group)), 2)))
  combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(my.design))))
  my.DA <- DA_proteins_apply(contrast.matrix, my.data, my.design, my.logFC, my.FDR, FALSE, my.block)
  return(my.DA)
}




# -----------------------------------------------------------------------------------------------------------------------------
# FUNCTION THAT GENERATES A COLOR MATRIX FOR HEATMAP.PLUS
# -----------------------------------------------------------------------------------------------------------------------------

get_colors <- function(my.truestatus, my.cols) {
  hm_col <- data.frame(levels(as.factor(as.character(my.truestatus))), my.cols)
  colnames(hm_col) <- c("status", "mycolor")
  true_status <- data.frame(my.truestatus)
  myorder <- 1:nrow(true_status)
  true_status$order <- myorder
  colnames(true_status) <- c("status", "order")
  col <- merge(true_status, hm_col, by="status", all.x =TRUE)
  col <- col[order(col$order),]
  col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor))
  return(as.matrix(col$mycolor))
}



# --------------------------------------------------------------------------------------
# Function for H-clustering of experssion/abundance data.
# The function takes arguments: dataset = abundance/expression data, nc = number of clusters desired,
# colorcode = vector of colors for labels and colorby = vector of lables/groups for coloring.


tree <- function(dataset, nc, colorcode, colorby) {
  d <- dist(t(dataset), method = "euclidean")
  fit_tree <- hclust(d, method="ward.D2") 
  dend <- as.dendrogram(fit_tree)
  #dend <- color_branches(dend, k=2, col=c("#F8766D", "#00BCD8"))
  par(cex=0.7, mar=c(18, 4, 2, 1))
  labels_colors(dend) <- colorcode[as.factor(colorby)][order.dendrogram(dend)]
  plot(dend, type="triangle", edgePar = list(lwd = 2))
  groups <- cutree(fit_tree, nc) 
  rect.hclust(fit_tree, nc, border="grey")
}



# --------------------------------------------------------------------------------------
# Function for visualization of Gene Ontology by barplots and corrplot

GO_visual <- function(my.proteins, my.univers, GO.background, nterms, filename, my.barp, my.corrp) {
  p <- sort(as.character(my.proteins))
  GO <- TOPGO("BP", my.univers, p, GO.background, nterms)
  if (my.barp == TRUE) {
    png(filename = paste0(filename, "_barplot.png"), height = 600, width = 1000)
    barplot_func(GO)
    dev.off()
    }
  if (my.corrp == TRUE) {
    png(filename = paste0(filename, "_corrplot.png"), height = 800, width = 1200)
    corrplot_func(GO) 
    dev.off()
    }
  return(GO)
}




GO_map_to <- function(humanGO, my.GOs, my.proteins.vector, my.gene.vector) {
  mappedunip <- lapply(my.GOs$GO.ID, function(x) intersect(sort(as.character(my.proteins.vector)), sort(as.character(humanGO[humanGO$V5 %in% x,]$V2))))
  mappedunip <- lapply(mappedunip, function(x) paste(x, collapse = ","))
  my.GOs$uniprot <- mappedunip
  mappedgene <- lapply(my.GOs$GO.ID, function(x) intersect(sort(as.character(my.gene.vector)), sort(as.character(humanGO[humanGO$V5 %in% x,]$V3))))
  mappedgene <- lapply(mappedgene, function(x) paste(x, collapse = ","))
  my.GOs$symbol <- mappedgene
  return(my.GOs)
}




# --------------------------------------------------------------------------------------
# PRINCIPAL COMPONENT ANALYSIS
# --------------------------------------------------------------------------------------


# Function that creates a rotation object

rotation_func <- function(my.data) {
  pca <- prcomp(t(my.data), center = TRUE, scale=TRUE)
  rotation <- data.frame(pca$rotation)
  return(rotation)
}


# --------------------------------------------------------------------------------------

# Function that creates barplot of PC values.
# Function takes arguments; my.data = abundance/expresssion data, my.group = vector of groups/lables and
# my.colors = vector of desired colors (must be the same length as levels of my.group)

boxplots <- function(my.data, my.group, my.colors) {
  pca <- prcomp(t(my.data), center = TRUE, scale=TRUE)
  pc <- data.frame(predict(pca, newdata = t(my.data)))
  data <- melt(pc[1:ncol(pc),], var='PC')
  group_reps <- data.frame(c(replicate(ncol(pc), my.group)))
  bp_data<- data.frame(group_reps, data)
  colnames(bp_data) <- c("group", "PC", "values")
  ggplot(data=bp_data, aes(x=PC, y=values, fill=group)) +  theme_bw() + geom_boxplot(stat="boxplot", position="dodge") + scale_fill_manual(values = my.colors) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=10, face="bold")) + theme(axis.text = element_text(colour = "black"))
}




# --------------------------------------------------------------------------------------
# PATHWAY ENRICHMENT ANALYSIS
# --------------------------------------------------------------------------------------



enrich_pathway <- function(background, my.u2e.map, my.proteins, my.name, my.plot) {
  entrez.data <- my.u2e.map[my.u2e.map$uniprot %in% my.proteins, ]
  my.pathway <- enrichPathway(as.character(entrez.data$entrez), organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.05, as.character(background$entrez), minGSSize = 3, readable = T)
  if (my.plot == TRUE) {
    png(filename = paste0(my.name,".png"), height = 800, width = 1200)
    cnetplot(my.pathway, categorySize="pvalue")
    dev.off()
  }
  return(my.pathway)
} 



# --------------------------------------------------------------------------------------
# Random forest
# --------------------------------------------------------------------------------------


my_forest_conver <- function(my.seed, my.data, my.groups) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=3000)
  plot(rf)
  return(rf)
}


my_forest <- function(my.seed, my.data, my.groups) {
  set.seed(my.seed)
  rfsel <- varSelRF(my.data, my.groups, ntree=3000, ntreeIterat=1500, vars.drop.frac=0.2)
  res <- rfsel$selected.vars
  return(res)
}





# --------------------------------------------------------------------------------------
# LASSO FUNCTION - Takes arguments, my.data = data marix of values, my.group = vector of integers specifying group.
# IF my.multinorm=TRUE, then analysis will be multinomial, IF my.multinorm=FALSE, then then analysis will be binomial
# --------------------------------------------------------------------------------------



LASSO_protein <- function(my.seed, my.data, my.group, my.multinorm=TRUE) {
  # orginal seed was 1011
  if(my.multinorm == TRUE) {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.group, family="multinomial", type.multinomial = "grouped", nfolds = 10)
    my.coef <- coef(my.fit, s=my.fit$lambda.min)
    my.ma <- as(my.coef$`1`, "matrix")
    rm(my.fit)
    rm(my.coef)
  }
  else {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.group, family = "binomial", type.measure = "class", nfolds = 10)
    my.coef <- coef(my.fit, s=my.fit$lambda.min)
    my.ma <- as(my.coef, "matrix")
    rm(my.fit)
    rm(my.coef)
  }
  my.ma <- names(my.ma[my.ma[,1] != 0, ])
  #my.ma <- my.data[rownames(my.data) %in% my.ma, ]
  return(my.ma[-1])
}




# --------------------------------------------------------------------------------------
# UpsetR plot - Intersection
# --------------------------------------------------------------------------------------



plot_upsetR <- function(list.of.sets, my.intersection, my.name, my.cols, my.plot, write.ids) {
  full.set <- data.frame(unique(sort(c(unlist(list.of.sets)))))
  colnames(full.set) <- "protein"
  for (name in  names(list.of.sets)) {
    full.set <- data.frame(full.set, ifelse(full.set$protein %in% as.character(list.of.sets[[name]]), 1, 0))
  }
  colnames(full.set) <- c("protein", names(list.of.sets))
  metadata <- data.frame("sets" = colnames(full.set)[-1], "sets2" = colnames(full.set)[-1])
  if (my.plot==TRUE) {
    pdf(paste0(my.name, ".pdf"), height = 6, width = 10)
    upset(full.set, sets=colnames(full.set)[2:ncol(full.set)], sets.bar.color = my.cols, set.metadata = list(data = metadata, plots = list(list(type="matrix_rows", column = "sets", colors = c("TIF" = my.cols[1], "Nglyco" = my.cols[2], "Plasma" = my.cols[3], "Secreted"=my.cols[4], "Exosomes"=my.cols[5]), alpha = 0.5))), order.by = "freq", text.scale = 1.7, keep.order = TRUE) 
    dev.off()
  }
  if (write.ids == TRUE) {
    idx <- which(names(list.of.sets) %in% my.intersection)
    write_out(Reduce(intersect, list.of.sets[idx]), my.name)
  } 
  return(full.set)
}



# --------------------------------------------------------------------------------------
# Uniprot accession to ensamble gene names
# --------------------------------------------------------------------------------------



uniprot_to_name <- function(my.vector) {
  ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
  annot1 <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol", "uniprot_gn"), filters="uniprotswissprot", values=my.vector, mart=ensembl)
  annot1$uniprot_gn <- NULL
  colnames(annot1) <- c("Accession", "name")
  annot2 <- hugo[hugo$Accession %in% my.vector, ]
  annot2[3] <- NULL
  annot <- rbind(annot1, annot2)
  annot <- unique(annot[order(annot$Accession),])
  return(annot) 
}



# --------------------------------------------------------------------------------------
# Write out as .txt
# --------------------------------------------------------------------------------------

write_out <- function(my.proteins, my.name) {
  if (class(my.proteins) == "data.frame") {
    my.proteins$Accession <- rownames(my.proteins)
    my.proteins$my.order <- 1:nrow(my.proteins)
  } else {
    my.proteins <- data.frame(my.proteins)
    colnames(my.proteins) <- "Accession"
    my.proteins$my.order <- 1:nrow(my.proteins)
  }
  my.genes <- uniprot_to_name(my.proteins$Accession)
  my.merged <- merge(my.proteins, my.genes, by ="Accession", all.x=TRUE)
  my.merged <- my.merged[order(my.merged$my.order),]
  my.merged$my.order <- NULL
  write.table(my.merged, paste0(my.name,".txt"), quote = FALSE, row.names = FALSE, sep = "\t")
  #return(my.merged)
}




# --------------------------------------------------------------------------------------
# CORRELATION NETWORKS
# --------------------------------------------------------------------------------------




my_networks <- function(data, direction, cuofcorr, cuoffdr) {
  combinations <- combs(as.character(rownames(data)), 2)
  
  df1_names <- data.frame(combinations[,1])
  colnames(df1_names) <- "first"
  df2_names <- data.frame(combinations[,2])
  colnames(df2_names) <- "second"
  
  
  temp_data <- data
  names <- rownames(temp_data)
  first <- rownames(temp_data)
  rownames(temp_data) <- NULL
  temp_data <- cbind(temp_data, names, first)
  
  
  df1 <- merge(df1_names, temp_data, by="first", sort = FALSE)
  df1 <- cbind(df2_names, df1)
  df1 <- with(df1, df1[order(second) , ])
  
  df2_names <- data.frame(df1$second)
  colnames(df2_names) <- "names"
  df1_names <- data.frame(df1$first)
  colnames(df1_names) <- "first"
  
  df1$first <- NULL
  df1$second <- NULL
  df1$names <- NULL
  
  
  df2 <- merge(df2_names, temp_data, by="names", sort = FALSE)
  
  df2$first <- NULL
  df2$names <- NULL
  
  
  df1 <- as.matrix(df1)
  df2 <- as.matrix(df2)
  
  spear_corr <- sapply(1:nrow(df1), function(i) cor(df1[i,], df2[i,], method = "spearman"))
  spear_corr <- data.frame(spear_corr)
  colnames(spear_corr) <- "cor_coef"
  
  # spearson correlation p-values
  spear_p_val <- sapply(1:nrow(df1), function(i) cor.test(df1[i,], df2[i,], method = "spearman")$p.value)
  spear_p_val <- data.frame(spear_p_val)
  colnames(spear_p_val) <- "pval"
  
  # correction for multiple testing with fdr
  fdr <- data.frame(p.adjust(spear_p_val$pval, method = "fdr"))
  colnames(fdr) <- "fdr"
  
  spear_corr_full <- cbind(df1_names, df2_names, spear_corr, fdr)
  
  spear_corr_full$cor_dir <- ifelse(spear_corr_full$cor_coef >= 0, "+", "-")
  spear_corr_full$cor_coef <- abs(spear_corr_full$cor_coef)
  
  spear_sig <- spear_corr_full[spear_corr_full$cor_coef >= cuofcorr & spear_corr_full$fdr < cuoffdr, ]
  
  max_min <- range(spear_sig$fdr[spear_sig$fdr != 0])
  spear_sig$fdr <- ifelse(spear_sig$fdr == 0, max_min[1], spear_sig$fdr)
  
  spear_sig$inverse_fdr <- -1*log(spear_sig$fdr)
  
  spear_sig <- merge(spear_sig, direction, by ="first")
  spear_sig$dir <- ifelse(spear_sig$cor_dir == '+' & spear_sig$direction == "up", "u", "o")
  spear_sig$dir <- ifelse(spear_sig$cor_dir == '+' & spear_sig$direction == "down", "d", spear_sig$dir)
  return(spear_sig)
}



get_nodes <- function(my.map, my.network, name) {
  nwk <- merge(my.network, my.map, by.x="first", by.y="Accession", all.x=TRUE)
  nwk <- merge(nwk, my.map, by.x="names", by.y="Accession", all.x=TRUE)
  nwk$node1 <- ifelse(is.na(nwk$name.x), as.character(nwk$first),  nwk$name.x)
  nwk$node2 <- ifelse(is.na(nwk$name.y), as.character(nwk$names),  nwk$name.y)
  nwk$name.x <- NULL
  nwk$name.y <- NULL
  write.table(nwk, paste0(name,"network.txt"), sep="\t", quote = FALSE, row.names = FALSE)
  return(nwk)
}




# --------------------------------------------------------------------------------------
# CHI-SQUARED TEST AND ODDS RATIO
# --------------------------------------------------------------------------------------

get_stats <- function(l1, l2, l3) {
  mytab <- data.frame(c(length(intersect(l1, l2)), length(l2)-length(intersect(l1, l2))), c(length(intersect(l1, l3)), length(l3)-length(intersect(l1, l3))))
  #chi <- chisq.test(mytab)$p.val
  #OR <- (mytab[1,1]*mytab[2,2])/(mytab[1,2]*mytab[2,1])
  #return(c(chi, OR))
  fisher <- fisher.test(mytab)
  return(fisher)
  }




# --------------------------------------------------------------------------------------
# 3Dplot
# --------------------------------------------------------------------------------------


my_3Dplot <- function(my.data, my.colors, my.r, my.names, is.PCAobejct) {
  if(is.PCAobejct == TRUE) {
    pca <- prcomp(t(my.data), center = TRUE, scale=TRUE)
    pc <- data.frame(predict(pca, newdata = t(my.data)))
  }
  else {
    pc <- data.frame(t(my.data))
    colnames(pc) <- c("PC1", "PC2", "PC3")
  }
  rgl.bg(color = "antiquewhite")
  rgl.spheres(pc$PC1, pc$PC2, pc$PC3, r=my.r, col= my.colors)
  lim <- function(x){c(-max(abs(x)), max(abs(x)))}
  rgl.lines(lim(pc$PC1), c(0,0), c(0,0), color= "grey50")
  rgl.lines(c(0,0), lim(pc$PC2),  c(0,0), color= "grey50")
  rgl.lines(c(0,0), c(0,0),  lim(pc$PC3),  color= "grey50")
  ellips <- ellipse3d(cov(cbind(pc$PC1, pc$PC2, pc$PC3)), centre = c(mean(pc$PC1), mean(pc$PC2), mean(pc$PC3)), level = 0.95)
  wire3d(ellips, col="grey50", lit=FALSE)
  aspect3d(1,1,1)
  axes <- rbind(c(lim(pc$PC1)[2], 0, 0), c(0, lim(pc$PC2)[2], 0), c(0, 0, lim(pc$PC3)[2]))
  rgl.texts(axes, text = c(my.names[1], my.names[2], my.names[3]), color = "grey50", adj = c(0.5, -0.8), size = 2)
}


# coo <- 1:720
# for(i in 1:length(coo)){
#   rgl.viewpoint(coo[i])
# }


my_3Dscatter <- function(my.data, my.colors, my.group, my.names, is.PCAobejct) {
  if(is.PCAobejct == TRUE) {
    pca <- prcomp(t(my.data), center = TRUE, scale=TRUE)
    pc <- data.frame(predict(pca, newdata = t(my.data)))
  }
  else {
    pc <- data.frame(t(my.data))
    colnames(pc) <- c("PC1", "PC2", "PC3")
  }
  greycol <- rep("grey60", nrow(my.data))
  with(pc, text3D(PC1, PC2, PC3, col = my.colors, theta = 10, phi=15, labels = my.group, cex = 0.8, bty = "b2", d = 3, clab = c("Urban","Pop"), adj = 0.6, font = 2, xlab=my.names[1], ylab=my.names[2], zlab=my.names[3]))
  with(pc, scatter3D(PC1, PC2, PC3, col=greycol, type = "h", add = TRUE))
}





# ----------------------------------------------------------------------------------------------------------------------------
# FUNCTION WHICH GENERATES AN OBJECT WITH ENTREZ GENE IDs AND ASSOCIATED LOGFC FROM A DE-TABLE
# ----------------------------------------------------------------------------------------------------------------------------



# Takes as arguments a DE_table from the DE_limma function and a u2e-map (uniprot-to-entrezID-map)

get_entrez <- function(my.table, my.u2e.map) {
  my.table <- data.frame(rownames(my.table), my.table$logFC)
  colnames(my.table) <- c("uniprot", "logFC")
  final.table <- merge(my.table, my.u2e.map, by = "uniprot")
  final.table <- final.table[order(final.table$entrez), ]
  logFC <- final.table$logFC
  names(logFC) <- final.table$entrez
  return(logFC)
}






# ---------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR KEGG PATHWAY ENRICHMENT ANALYSIS
# ---------------------------------------------------------------------------------------------------------------------------


keggs_entrez<- function(my.table) {
  ensembl <- useMart('ensembl', dataset="hsapiens_gene_ensembl")
  annot <- getBM(attributes=c("uniprot_swissprot", "entrezgene", "uniprot_genename"), filters="uniprot_swissprot", values=rownames(my.table), mart=ensembl)
  annot <- annot[match(rownames(my.table), annot$uniprot_swissprot),]
  annot <- cbind(my.table, annot)
  annot <- annot[order(as.numeric(annot$entrezgene)), ]
  logFC <- annot$logFC
  names(logFC) <- annot$entrezgene
  return(logFC)
}



my.keggs <- function(my.geneset, my.keggset) {
  keggres <- gage(my.geneset, gsets=my.keggset, same.dir=FALSE)
  keggpathways <- data.frame(id=rownames(keggres$greater), keggres$greater)
  keggpathways <- keggpathways[!is.na(keggpathways$p.val),]
  keggpathways <- as.character(keggpathways[keggpathways$p.val< 0.05,]$id)
  keggids <- substr(keggpathways, start=1, stop=8)
  pv.out.list <- sapply(keggids, function(pid) pathview(gene.data=both, pathway.id=pid, species="hsa", gene.idtype="KEGG", both.dirs = TRUE))
  return(pv.out.list)
}





# --------------------------------------------------------------------------------------
# FUNCTION FOR CLUSTER ANALYSIS
# optimal_nc function takes as arguments:
# a dataframe of expression/abundance values
# plot_clusters function takes as arguments:
# a dataframe of expression/abundance values and number of clusters
# --------------------------------------------------------------------------------------


# Calculating number of optimal clusters - clustGap 

optimal_nc <- function (my.dataframe) {
  pam1 <- function(my.dataframe,k) list(cluster = pam(my.dataframe,k, cluster.only=TRUE))
  optimal_number <- clusGap(t(my.dataframe), FUN = pam1, K.max = 15, B = 500)
  plot(optimal_number)
}


# Plotting clusters with patient IDs

plot_clusters <- function(my.dataframe, my.clusters) {
  d <- dist(t(my.dataframe))
  clust <- kmeans(d, my.clusters)
  library(fpc)
  library(cluster)
  clusplot(as.matrix(d), clust$cluster, color=TRUE, shade=TRUE, col.txt=col, col.clus = c("dodgerblue3", "magenta","springgreen2"), col.p=col, labels=2, cex=0.6, lines=0, main="CLUSPLOT")
}
