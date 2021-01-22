# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


# LOAD PACKAGES
library(pamr)
library(openxlsx)
library(tidyverse)
library(reshape)
library(stringr)
library(plyr)

library(limma)
library(sva)

library(gdata) 
library(topGO)
library(biomaRt)
library(GOSim)
#library(pvclust)

library(glmnet)
library(fitdistrplus)
library(caTools)
library(randomForest)
library(glmnet)
library(e1071)
library(caret)
library(varSelRF)

library(ggplot2)
library(UpSetR)
library(corrplot)
library(gplots)
library(dendextend)
library(heatmap.plus)
library(VennDiagram)
library(RColorBrewer)

#library(gage)
#library(pathview)



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GO OBJECT FUNCTION:
# Takes as arguments;
    # univ = a protein univers
    # intprot = a set of proteins of interest and
    # GO_background = GO-background (dataframe where proteins are coupled to GO-terms)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



GOobject <-function(ont, univ, intpro, GO_background){
  if(is.null(univ)){
    univ=(unique(unlist(GO_background)))
  }
  # factorise genelist
  geneList = factor(as.integer(univ%in%intpro))
  names(geneList) = univ
  #new topgo object
  GOdata = new("topGOdata", ontology = as.character(ont), nodeSize = 5, allGenes = geneList, annot = annFUN.GO2genes, GO2genes=GO_background)
  return(GOdata)
}  



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# TOPGO FUNCTION:
# Takes as arguments;
    # ont = the ontology to enrich for
    # univ = a protein univers
    # intprot = a set of proteins of interest,
    # GO_background = GO-background (list object where proteins are coupled to GO-terms) and nterms = number of GO hits to display.

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


TOPGO <- function(ont, univ, intpro, GO_background, nnodes) {
  
  GOdata <- GOobject(ont, univ, intpro, GO_background)
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
  resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
  allRes <- GenTable(GOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, topNodes=nnodes)
  
  return(allRes)
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function that creates barplot. Input is results of TOPGO

barplot_func <- function(res) {
  temp <- data.frame(res$Term, log2(1/as.numeric(res$classicFisher)))
  colnames(temp) <- c("Term", "Inverse.Rank")
  temp$Term <- factor(temp$Term, levels=as.vector(temp$Term))
  ggplot(data=temp, aes(x=Term, y=Inverse.Rank)) + geom_bar(stat="identity", fill="#00BCD8") + theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=12))
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function that creates correlation plot. Input is results of TOPGO

corrplot_func <- function(res) {
  term_corr <- getTermSim(res$GO.ID)
  rownames(term_corr) <- res$Term
  colnames(term_corr) <- NULL
  term_corr[which(term_corr==Inf)] <- 0 
  corrplot(term_corr, type = "lower", is.corr = FALSE, order = "hclust", tl.cex=1, tl.col = "black")
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# GENE ONTOLOGY FUNCTION:
# takes as arguments:
# my.proteins = list of uniprot IDs
# my.univers = a protein univers
# GO_background = GO-background (list object where proteins are coupled to GO-terms) and nterms = number of GO hits to display.
# nterms = number of terms to return
# filename = name of ourput
# if my.barp = TRUE and/or my.corrp = TRUE a barplot and a corrplot will be generated
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



GO_visual <- function(my.proteins, my.univers, GO.background, nterms, filename, my.barp, my.corrp) {
  p <- sort(as.character(my.proteins))
  GO <- TOPGO("BP", my.univers, p, GO.background, nterms)
  if (my.barp == TRUE) {
    pdf(filename = paste0(filename, "_barplot.pdf"), height = 600, width = 1000)
    barplot_func(GO)
    dev.off()
  }
  if (my.corrp == TRUE) {
    pdf(filename = paste0(filename, "_corrplot.pdf"), height = 800, width = 1200)
    corrplot_func(GO) 
    dev.off()
  }
  return(GO)
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

GO_map_to <- function(humanGO, my.GOs, my.proteins.vector, my.gene.vector) {
  mappedunip <- lapply(my.GOs$GO.ID, function(x) intersect(sort(as.character(my.proteins.vector)), sort(as.character(humanGO[humanGO$V5 %in% x,]$V2))))
  mappedunip <- do.call(rbind,lapply(mappedunip, function(x) paste(x, collapse = ",")))
  my.GOs$uniprot <- as.character(mappedunip[,1])
  mappedgene <- lapply(my.GOs$GO.ID, function(x) intersect(sort(as.character(my.gene.vector)), sort(as.character(humanGO[humanGO$V5 %in% x,]$V3))))
  mappedgene <- do.call(rbind,lapply(mappedgene, function(x) paste(x, collapse = ",")))
  my.GOs$symbol <- as.character(mappedgene[,1])
  return(my.GOs)
}
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DATA CORRECTION AND VISUALIZATION
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Batch correction of data:
# Takes as arguments;
    # my.data = abundance/expression data
    # my.group = patient groups,
    # my.batch = technical batches


combat_corrections <- function(my.data, my.batch, my.group=NULL){
  if(!is.null(my.group)) {
    mod_design <- model.matrix(~as.factor(my.group))
    batch_corr <- ComBat(dat = as.matrix(my.data), batch = as.factor(my.batch), mod = mod_design, par.prior=TRUE,prior.plots=FALSE)
  } else {
    batch_corr <- ComBat(dat = as.matrix(my.data), batch = as.factor(my.batch), par.prior=TRUE,prior.plots=FALSE)
  }
}


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# MULTIDIMENTIONAL SCALING PLOT:
# Takes as arguments;
# my.data = a dataframe of expression/abundance counts
# my.group, my.labels = a vector of IDs for coloring and labeling (may be the same or different, length should be equal to ncol(dataframe))
# my.cols = a vector of colors (one color for each group)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


myMDSplot <- function(my.data, my.group, my.labels, my.cols) {
  d<-dist(t(my.data))
  fit <- cmdscale(d,eig=TRUE, k=2)
  res<-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
  p <- ggplot(data=res)
  p + geom_point(aes(x=M1,y=M2,color=my.group), size=2) + geom_text(aes(x=M1,y=M2, label= my.labels, color=my.group)) + scale_color_manual(values  = my.cols) +
    coord_cartesian(xlim=c(min(res$M1)*1.4,max(res$M1)*1.4)) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DIFFERENTIAL EXPRESSION ANALYSIS
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# FUNCTION TO OBTAIN DIFFERENTIALLY ABUNDANT PROTEINS:
# Takes as arguments;
# a contrast between groups of interest
# a dataframe, a design matrix with all comparisons
# cutoffs for logFC and FDR
# if blocking than a vector of patient IDs
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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




# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO APPLY DIFFERENTIALLY ABUNDANCE ANALYSIS TO ALL COMPARISONS AT ONCE:
# Takes as arguments;
# all contrasts between groups of interest
# a dataframe
# a design matrix with all comparisons
# cutoffs for logFC and FDR
# if blocking than a vector of patient IDs
# TRUE/FALSE statment specifying output format, if TRUE the function return a vector of protein IDs only
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DA ANALYSIS WITH CLINICAL PARAMETERS. THE FUNCTION CALLS THE DA_proteins_apply FROM ABOVE.
# Takes as arguments;
# a dataframe
# a vector of groups do perform contrasts on (same length as ncol(dataframe))
# a cutoff for logFC and FDR
# if remove is different from NULL, a vector of indices to remove must be supplied
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



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



# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR GETTING HEATMAP COLORS
# Takes as arguments:
# my.truestatus = a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
# my.cols = a vector with colors to use (a character vector with the length of the number of groups/levels).
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BOXPLOT VISUALIZATION:
# takes as arguments:
    # my.data = abundance/expresssion data
    # my.group = vector of groups/lables and
    # my.colors = vector of desired colors (must be the same length as levels of my.group)
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

boxplots <- function(my.data, my.group, my.colors) {
  pca <- prcomp(t(my.data), center = TRUE, scale=TRUE)
  pc <- data.frame(predict(pca, newdata = t(my.data)))
  data <- melt(pc[1:ncol(pc),], var='PC')
  group_reps <- data.frame(c(replicate(ncol(pc), my.group)))
  bp_data<- data.frame(group_reps, data)
  colnames(bp_data) <- c("group", "PC", "values")
  ggplot(data=bp_data, aes(x=PC, y=values, fill=group)) +  theme_bw() + geom_boxplot(stat="boxplot", position="dodge") + scale_fill_manual(values = my.colors) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + guides(colour = guide_legend(override.aes = list(size=6))) + theme(legend.position = "top") + theme(axis.text=element_text(size=10, face="bold")) + theme(axis.text = element_text(colour = "black"))
}






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR LASSO REGRESSION
# Takes as arguments:
# my.seed = A random seed
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# my.validation = TRUE or FALSE
# my.multinorm = TRUE or FALSE
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

LASSO_Features <- function(my.seed, my.data, my.group, my.validation=FALSE, my.multinorm=TRUE) {
  
  if (my.validation == TRUE) {
    
    ll <- list()
    llev <- levels(as.factor(my.group))
    
    for (idx in 1:length(llev)) {
      pos <- which(my.group == as.character(llev[idx]))
      ll[[idx]] <- pos
    }
    
    set.seed(5)
    my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))
    
    
    my.data.val <- my.data[,my.samp]
    my.group.val <- my.group[my.samp]
    
    my.data <- my.data[,-my.samp]
    my.group <- my.group[-my.samp]
  }
  
  if(my.multinorm == TRUE) {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.group, family="multinomial", type.multinomial = "grouped", nfolds = 10, alpha = 0.9, parallel=TRUE, keep = TRUE)
    my.coef <- coef(my.fit, s=my.fit$lambda.1se)
    my.ma <- as(my.coef[[1]], "matrix")
  } else {
    set.seed(my.seed)
    my.fit <- cv.glmnet(x = t(my.data), y = my.group, family = "binomial", type.measure = "class", nfolds = 10, alpha = 0.9, keep = TRUE)
    my.coef <- coef(my.fit, s=my.fit$lambda.1se)
    my.ma <- as(my.coef, "matrix")
  }
  
  my.err <- as.numeric(my.fit$cvm[my.fit$lambda == my.fit$lambda.min])
  my.imp <- my.ma
  my.vars <- names(my.ma[my.ma[,1] != 0, ])
  
  if (my.validation == TRUE) {
    my.pred <- predict(my.fit, t(my.data.val), s=my.fit$lambda.1se, type="class")
    my.group.val <- as.factor(my.group.val)
    my.pred <- factor(my.pred[,1], levels = levels(my.group.val))
    my.pred <- confusionMatrix(my.pred, my.group.val)
    res <- list(my.fit, my.vars, my.imp, my.err, my.pred)
    
  } else {
    res <- list(my.fit, my.vars, my.imp, my.err, NA)
  }
  return(res)
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR APPLYING LASSO OVER DIFFERENT SEEDS.
# Takes as arguments:
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# my.validation = TRUE or FALSE
# my.multinorm = TRUE or FALSE
# b1 and b2 = two batch variables
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

LASSO_Apply <- function(my.data, my.group, my.validation, my.multinorm, b1=NULL, b2=NULL){
  
  if(is.null(b1) & is.null(b2)) {
    my.data <- my.data
  } else {
    my.data <- combat_corrections(my.data, my.group, b1, b2)
  }
  
  seeds <- c(1, 40, 99, 640, 6565)
  #seeds <- c(1, 41, 100, 641, 6566)
  
  LAvars <- list()
  LAimp <- list()
  LAerr <- list()
  #LAclass <- list()
  LAacc <- list()
  
  
  for (idx in 1:length(seeds)) {
    LA <- LASSO_Features(seeds[[idx]], my.data, my.group, my.validation, my.multinorm)
    LAvars[[idx]] <-  LA[[2]]
    imp <-  data.frame(LA[[3]])
    imp[,1] <- abs(imp[,1])
    LAimp[[idx]] <- imp
    LAerr[[idx]] <-  LA[[4]]
    
    if(my.validation == TRUE) {
      LAacc[[idx]] <-  LA[[5]]$overall[c(3,1,4)]
      #LAclass[[idx]] <-  as.numeric(LA[[5]]$class)
    } else {
      #LAclass[[idx]] <-  NA
      LAacc[[idx]] <- NA
    }
    print(paste0(idx, " out of ", length(seeds), " runs completed."))
  }
  
  
  all.df <- LAimp[[1]]
  
  for (idx in 2:length(LAimp)) {
    all.df <- all.df+LAimp[[idx]]
  }
  
  LAimp <- all.df/length(LAimp)
  LAimp$Accession <- rownames(LAimp)
  
  LAvars <- unique(unlist(LAvars))[-1]
  
  LAimp <- LAimp[LAimp$Accession %in% LAvars,]
  LAimp <- LAimp[order(LAimp$X1, decreasing = TRUE),]
  
  LAres <- list(LAvars, LAimp, mean(unlist(LAerr)), do.call(rbind, LAacc))
  names(LAres) <- c("Vars", "Importance", "CV.Error", "Accuracy")
  #LAres <- list(LAvars, LAimp, mean(unlist(LAerr)), unlist(LAclass), mean(unlist(LAclass)), do.call(rbind, LAacc))
  #names(LAres) <- c("Vars", "Importance", "CV.Error", "Class.Err", "Class.Err.Mean", "Accuracy")
  return(LAres)
}    




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR AVERAGING VARIABLES ACROSS LASSO RUNS.
# Takes as arguments:
# list.of.lassos = list of LASSO regressions
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


AverageWeight <- function(list.of.lassos) {
  
  if(class(list.of.lassos[[1]])=="list") {
    all.lassos <- list.of.lassos[[1]]$Importance
    
    for (idx in 2:length(list.of.lassos)) {
      all.lassos <- merge(all.lassos, list.of.lassos[[idx]]$Importance, by = "Accession", all.x = TRUE, all.y =TRUE)
    }
  } else {
    all.lassos <- list.of.lassos[[1]]
    
    for (idx in 2:length(list.of.lassos)) {
      all.lassos <- merge(all.lassos, list.of.lassos[[idx]], by = "Accession", all.x = TRUE, all.y =TRUE)
    }
  }
  all.lassos$weight <- rowMeans(all.lassos[,-1], na.rm = TRUE)
  all.lassos <- all.lassos[order(all.lassos$weight, decreasing = TRUE),]
  return(all.lassos)
}




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR RANDOM FOREST
# Takes as arguments:
# my.seed = A random seed
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# my.validation = TRUE or FALSE
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Forest_Features <- function(my.seed, my.data, my.group, my.validation=FALSE) {
  
  if (my.validation == TRUE) {
    
    ll <- list()
    llev <- levels(as.factor(my.group))
    
    for (idx in 1:length(llev)) {
      pos <- which(my.group == as.character(llev[idx]))
      ll[[idx]] <- pos
    }
    
    set.seed(5)
    my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))
    
    my.data.val <- t(my.data[,my.samp])
    my.group.val <- my.group[my.samp]
    
    my.data <- t(my.data[,-my.samp])
    my.group <- my.group[-my.samp]
    
    set.seed(my.seed)
    RFvars <- varSelRF(my.data, my.group, ntree=3000, ntreeIterat=1500, vars.drop.frac=0.2)
    set.seed(my.seed)
    RFforest <- randomForest(my.data, my.group, ntree=3000, ntreeIterat=1500, vars.drop.frac=0.2)
    pred <- predict(RFforest, newdata = my.data.val)
    RFpred <- confusionMatrix(pred, my.group.val)
    
    res <- list(RFvars, RFpred, RFforest)
    
  } else {
    my.data <- t(my.data)
    set.seed(my.seed)
    RFvars <- varSelRF(my.data, my.group, ntree=3000, ntreeIterat=1500, vars.drop.frac=0.2)
    
    res <- list(RFvars, NA)
  }
  return(res)
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR APPLYING RANDOM FOREST OVER DIFFERENT SEEDS.
# Takes as arguments:
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# my.validation = TRUE or FALSE
# b1 and b2 = two batch variables
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

RF_Apply <- function(my.data, my.group, my.validation, b1=NULL, b2=NULL){
  
  if(is.null(b1) & is.null(b2)) {
    my.data <- my.data
  } else if(!is.null(b1) & is.null(b2)) {
    my.data <- combat_corrections(my.data, my.group, b1)
  } else if(is.null(b1) & !is.null(b2)) {
    my.data <- combat_corrections(my.data, my.group, b2)
  } else {
    my.data <- combat_corrections(my.data, my.group, b1, b2)
  }
  
  seeds <- c(1, 40, 99, 640, 6565)
  #seeds <- c(1, 41, 100, 641, 6566)
  
  RFvars <- list()
  RFobb <- list()
  RFclass <- list()
  RFacc <- list()
  RFimp <- list()
  
  for (idx in 1:length(seeds)) {
    RF <- Forest_Features(seeds[[idx]], my.data, my.group, my.validation)
    RFvars[[idx]] <-  RF[[1]]$selected.vars
    RFobb[[idx]] <-  median(RF[[1]]$firstForest$err.rate[,1])
    RFclass[[idx]] <-  RF[[1]]$firstForest$confusion
    RFimp[[idx]] <-  RF[[1]]$firstForest$importance
    
    if(class(RF[[2]]) == "confusionMatrix") {
      RFacc[[idx]] <- RF[[2]]$overall[c(1,3,4)]
    } else {
      RFacc[[idx]] <- NA
    }
    print(paste0(idx, " out of ", length(seeds), " runs completed."))
  }
  
  
  all.df <- RFimp[[1]]
  
  for (idx in 2:length(RFimp)) {
    all.df <- all.df+RFimp[[idx]]
  }
  
  RFimp <- all.df/length(RFimp)
  
  
  RFres <- list(unlist(RFvars), unlist(RFobb), RFclass, do.call(rbind, RFacc), RFimp)
  names(RFres) <- c("Vars", "Obb", "Class", "Accuracy", "Importance")
  return(RFres)
}    



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR AVERAGING VARIABLES ACROSS RF RUNS.
# Takes as arguments:
# list.of.rfs = list of random forests
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

AverageWeight <- function(list.of.rfs) {
  
  if(class(list.of.rfs[[1]])=="list") {
    
    list.of.dfs <- lapply(list.of.rfs, function(x) x$Importance)
    list.of.rfs <- list()
    
    for (idx in 1:length(list.of.dfs)) {
      df <- data.frame(list.of.dfs[[idx]])[]
      df$Accession <- rownames(df)
      n <- ncol(df)
      list.of.rfs[[idx]] <- df[, c(n,(n-1))]
    }
  }
  all.rfs <- list.of.rfs[[1]]
  for (idx in 2:length(list.of.rfs)) {
    all.rfs <- merge(all.rfs, list.of.rfs[[idx]], by="Accession", all.x = TRUE, all.y =TRUE)
  }
  all.rfs$weight <- rowMeans(all.rfs[,-1], na.rm = TRUE)
  all.rfs <- all.rfs[order(all.rfs$weight, decreasing = TRUE),]
  return(all.rfs)
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PATHWAY ENRICHMENT ANALYSIS:
# takes as arguments:
    # background = protein background for enrichment
    # my.u2e.map = a map of uniprot IDs to entrez IDs
    # my.proteins = list of uniprot IDs of interest
    # my.name =  name of pdf file
    # if my.plot = TRUE pdf will be generated
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



enrich_pathway <- function(background, my.proteins, my.name, my.plot) {
  ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
  entrez.data <- unique(sort(as.character(getBM(attributes=c("uniprotswissprot", "entrezgene_id"), filters="uniprotswissprot", values=my.proteins, mart=ensembl)$entrezgene_id)))
  my.pathway <- enrichPathway(entrez.data, organism = "human", pvalueCutoff = 0.05, pAdjustMethod = "fdr", qvalueCutoff = 0.05, background, minGSSize = 3, readable = T)
  if (my.plot == TRUE) {
    pdf(paste0(my.name,".pdf"), height = 8, width = 12)
    p <- cnetplot(my.pathway, categorySize="pvalue")
    print(p)
    dev.off()
  }
  return(my.pathway)
} 



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR RANDOM FOREST CONVERGENCE
# Takes as arguments:
# my.seed = A random seed
# my.data = matrix of countes/expression/abundance
# my.groups = vector of group IDs, must be as.factor()
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



my_forest_conver <- function(my.seed, my.data, my.groups) {
  set.seed(my.seed)
  rf <- randomForest(x=my.data, y=my.groups, ntree=3000)
  plot(rf)
  return(rf)
}


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Plot UpSetR
# Takes as arguments;
# list.of.sets = list with sets to plot
# my.intersection = the names of the sets to intersect
# my.name = name of output plot
# my.cols = colors vector with as, one color per set
# if my.plot= TRUE a pdf is writen out, og write.ids = TRUE, write out intersection of all sets
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plot_upsetR <- function(list.of.sets, my.intersection, my.name, my.cols, my.plot, write.ids) {
  full.set <- data.frame(unique(sort(c(unlist(list.of.sets)))))
  colnames(full.set) <- "Accession"
  for (name in  names(list.of.sets)) {
    full.set <- data.frame(full.set, ifelse(full.set$Accession %in% as.character(list.of.sets[[name]]), 1, 0))
  }
  colnames(full.set) <- c("Accession", names(list.of.sets))
  metadata <- data.frame("sets" = colnames(full.set)[-1], "sets2" = colnames(full.set)[-1])
  if (my.plot==TRUE) {
    pdf(paste0(my.name, ".pdf"), height = 6, width = 10)
    my.combination <- my.cols
    names(my.combination) <- my.intersection
    p <- upset(full.set, sets=colnames(full.set)[2:ncol(full.set)], sets.bar.color = my.cols, set.metadata = list(data = metadata, plots = list(list(type="matrix_rows", column = "sets", colors = my.combination, alpha = 0.5))), order.by = "freq", text.scale = 1.7, keep.order = TRUE) 
    print(p)
    dev.off()
  }
  if (write.ids == TRUE) {
    idx <- which(names(list.of.sets) %in% my.intersection)
    write_out(Reduce(intersect, list.of.sets[idx]), my.name)
  } 
  return(full.set)
}


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# UNIPROT ACCESSION TO GENE NAME
# Takes as arguments;
    # my.vector = vector uniprot IDs
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



uniprot_to_name <- function(my.vector) {
  ensembl <- useMart('ensembl', dataset = "hsapiens_gene_ensembl")
  annot1 <- getBM(attributes=c("uniprotswissprot", "hgnc_symbol", "uniprot_gn_id"), filters="uniprotswissprot", values=my.vector, mart=ensembl)
  annot1$uniprot_gn_id <- NULL
  colnames(annot1) <- c("Accession", "name")
  annot2 <- hugo[hugo$Accession %in% my.vector, ]
  annot2[3] <- NULL
  annot <- rbind(annot1, annot2)
  annot <- unique(annot[order(annot$Accession),])
  return(annot) 
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# WRITE OUT
# Takes as arguments;
    # my.proteins = vector of uniprot IDs OR dataframe of expression/abundance values with IDs as rownames
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# miRNA-miRNA co-expression network with spearman analysis.
# Takes as arguments:
# data = a dataframe with expression/abundance counts
# direction = vector specifying if the given miRNA is up -or down regulated in the given comparison ("up" or "down")
# coCorr, coFDR = cut-off for significant fdr and correlation
# name = name of output file (as string)
# IF genes=NULL then only miRNAs will be used, or the user may provide a list of miRNA-gene targets.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




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




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# CHI-SQUARED TEST AND ODDS RATIO
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

get_stats <- function(l1, l2, l3) {
  mytab <- data.frame(c(length(intersect(l1, l2)), length(l2)-length(intersect(l1, l2))), c(length(intersect(l1, l3)), length(l3)-length(intersect(l1, l3))))
  #chi <- chisq.test(mytab)$p.val
  #OR <- (mytab[1,1]*mytab[2,2])/(mytab[1,2]*mytab[2,1])
  #return(c(chi, OR))
  fisher <- fisher.test(mytab)
  return(fisher)
  }



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ENTREZ GENE IDs AND ASSOCIATED LOGFC FROM A DE-TABLE
# Takes as arguments;
    # my.table = DE_table from the DE_limma function
    # my.u2e.map = uniprot to entrezID map
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


get_entrez <- function(my.table, my.u2e.map) {
  my.table <- data.frame(rownames(my.table), my.table$logFC)
  colnames(my.table) <- c("uniprot", "logFC")
  final.table <- merge(my.table, my.u2e.map, by = "uniprot")
  final.table <- final.table[order(final.table$entrez), ]
  logFC <- final.table$logFC
  names(logFC) <- final.table$entrez
  return(logFC)
}






# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR KEGG PATHWAY ENRICHMENT ANALYSIS:
# Takes as arguments;
    # my.table = DE_table from the DE_limma function
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR CLUSTER ANALYSIS
# optimal_nc function takes as arguments:
# a dataframe of expression/abundance values
# plot_clusters function takes as arguments:
# a dataframe of expression/abundance values and number of clusters
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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



Overlap_GOplot <- function(my.set, my.lfc.cols, my.rib.cols) {
  genenames <- uniprot_to_name(my.set$Accession)
  my.set <- merge(my.set, genenames, by ="Accession", all.x=TRUE, all.y=FALSE)
  my.set <- my.set[!is.na(my.set$name),]
  rownames(my.set) <- my.set$name
  my.set$name <- NULL
  my.set$Accession <- NULL
  my.set$logFC <- rnorm(nrow(my.set), 0, 2)
  remove <- grep("^$", rownames(my.set))
  if(length(remove) > 0) {
    my.set <- my.set[-remove,]
  }
  GOChord(as.matrix(my.set), lfc.col = my.lfc.cols, ribbon.col = my.rib.cols, border.size = 0.1)
}


# Get p-values of Ordinal Logistic Regression
OLR_pval <- function(model) {
  ctable <- coef(summary(model))
  p <- pnorm(abs(ctable[, "t value"]), lower.tail = FALSE) * 2
  #p <- p.adjust(p, method = "fdr", n = 10)
  ctable <- cbind(ctable, "p value" = p)
  return(ctable)
}
