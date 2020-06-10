
my.wd <- ""

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR DOWNLOADING STRING DATABASE
# Take arguments:
# my.geneIDs = An approved gene-ID for conversion. Options are "ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id" or "uniprotswissprot".
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


DownloadPPInt <- function(my.geneIDs, my.version = "11.0") {
  approvedGeneIDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
  print("\nDownloading and preparing string database for protein-protein interactions - this may take a couple of minutes!\n")
  download.file(paste0("https://stringdb-static.org/download/protein.links.v",my.version,"/9606.protein.links.v",my.version,".txt.gz"), paste0("9606.protein.links.v", my.version, ".txt.gz"))
  stringDB <- data.table::fread(paste0("9606.protein.links.v", my.version, ".txt.gz"))
    
  stringDB$protein1 <- gsub("9606.", "", stringDB$protein1)
  stringDB$protein2 <- gsub("9606.", "", stringDB$protein2)
  
  uqprot <- unique(c(stringDB$protein1, stringDB$protein2))
  uqprot <- split(uqprot, ceiling(seq_along(uqprot)/10000))
  print("Calling IDs from BiomaRt")
  ensemblMT <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  #map <- unique(getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = uqprot, mart = ensemblMT))
  map <- lapply(uqprot, function(x) getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = x, mart = ensemblMT))
  map <- do.call("rbind", map)
  colnames(map) <- c("protein1", "ID1")
  stringDB <- merge(stringDB, map, by = "protein1")
  colnames(map) <- c("protein2", "ID2")
  stringDB <- merge(stringDB, map, by = "protein2")
  stringDB <- stringDB[,-c(1,2)]
  
  stringDB <- stringDB[stringDB$combined_score > as.numeric(quantile(stringDB$combined_score)[2]),]
  save(stringDB, file = "stringDB.Rdata")
  return(stringDB)
}





# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR SUBSETTING DATAFRAME FROM DAA INTO A LIST OF DATAFRAMES BY COMPARISONS.
# Take arguments:
# my.data = a dataframe with results of differential expression analysis.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DFtoList <- function(my.data) {
  
  my.data$name <- gsub("_", "-", my.data$name)
  comparisons <- levels(as.factor(my.data$comparison))
  df.list <- list()
  
  for (idx in 1:length(comparisons)) {
    df <- my.data[my.data$comparison == as.character(comparisons[idx]),]
    df.list[[idx]] <- df[,-c(1,3:4,6)]
  }
  
  names(df.list) <- comparisons
  
  df.lengths <- as.numeric(unlist(lapply(df.list, function(x) nrow(x))))
  df.lengths.min <- which(df.lengths < 2)
  if (length(df.lengths.min) > 0) {
    df.list <- df.list[-df.lengths.min]
  }
  return(df.list)
}




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR EXTRACTING PROTEIN-PROTEIN INTERACTIONS FOR DE/DA GENES/PROTEINS 
# Take arguments:
# my.DB = output of the function "DownloadPPInt", e.g. a curated string database.
# my.Geneset = a dataframe with results of differential expression analysis.
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



GetDeaPPInt <- function(my.DB, my.Genedat) {
  
  df.list <- DFtoList(my.Genedat)
  nodes.list <- list()
  
  for (idx in 1:length(df.list)) {
    # Merging StringDB and datasets
    df <- df.list[[idx]]
    df$ID1 <- df$name
    nodes <- merge(my.DB, df, by = "ID1")
    df$ID2 <- df$name
    nodes <- unique(merge(nodes, df, by = "ID2"))
    nodes <- nodes[order(nodes$combined_score, decreasing = TRUE),]
    df <- as.data.frame(nodes[, c(2,1,3,4,5,6,9,10,11,13)])
    colnames(df) <- c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "logFC.node2", "fdr.node2", "dir.node2", "comparison")
    df <- df[!duplicated(apply(df,1,function(x) paste(sort(x),collapse=''))),]
    df$dir <- rep(NA, nrow(df))
    df$dir <- ifelse(df$dir.node1 == "up" & df$dir.node2 == "up", "up", df$dir)
    df$dir <- ifelse(df$dir.node1 == "down" & df$dir.node2 == "down", "down", df$dir)
    df$dir <- ifelse(is.na(df$dir), "inverse", df$dir)
    nodes.list[[idx]] <- df
  }
  
  names(nodes.list) <- names(df.list)
  return(nodes.list)
}




# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION WHICH EXTRACTS NODE INFORMATION
# Take arguments:
# my.trimmed.list = List of trimmed networks, i.e. output of the function "TrimWriteInt".
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



NodeInfo <- function(my.nodes.list) {
  
  node.lengths <- as.numeric(unlist(lapply(my.nodes.list, function(x) nrow(x))))
  node.lengths.min <- which(node.lengths < 2)
  
  if (length(node.lengths.min) > 0) {
    my.nodes.list <- my.nodes.list[-node.lengths.min]
  }
  
  set.names <- names(my.nodes.list)
  trimlist <- list()
  
  for (idx in 1:length(my.nodes.list)) {
    p.nodes <- my.nodes.list[[idx]]
    p.nodes$myrank <- rowMeans(apply(cbind(abs(p.nodes$logFC.node1), abs(p.nodes$logFC.node2)), 2, rank))
    p.nodes <- p.nodes[order(p.nodes$myrank, decreasing = TRUE),]
    
    p.info <- data.table(c(as.character(p.nodes$node1), as.character(p.nodes$node2)), c(p.nodes$logFC.node1, p.nodes$logFC.node2))
    colnames(p.info) <- c("name", "logFC")
    p.info<- data.frame(p.info[, .(Freq = .N), by = .(name, logFC)])
    p.info <- p.info[order(p.info$Freq, decreasing = TRUE),]
    p.info$group <- 1:nrow(p.info)
    vircls <- viridis(2, direction = -1 , end = 0.9, option = "cividis")
    #vircls <- viridis(2, end = 0.6, direction = -1, option = "magma")
    p.info$calfb <- ifelse(p.info$logFC > 0, vircls[1], "grey50")
    p.info$calfb <- ifelse(p.info$logFC < 0, vircls[2], as.character(p.info$calfb))

    myorder <- unique(as.vector(t(data.frame(p.nodes[,1:2]))))
    p.info <- p.info[match(myorder, p.info$name),]
    
    write.table(p.nodes, paste0(set.names[[idx]],"_AllInteractions.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
    write.table(p.info, paste0(set.names[[idx]],"_InfoNodes.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)
    
    trimmed <- list(p.nodes, p.info)
    names(trimmed) <- c("p.nodes", "p.info")
    trimlist[[idx]] <- trimmed
  }
  names(trimlist) <- set.names
  return(trimlist)
}

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FUNCTION FOR PLOTTING INTERACTION NETWORKS
# Take arguments:
# my.trimmed.list = List of trimmed networks, i.e. output of the function "TrimWriteInt".
# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



PlotInt <- function(my.trimmed.list) {
  
  trim.lengths <- as.numeric(unlist(lapply(my.trimmed.list, function(x) nrow(x))))
  trim.lengths.min <- which(trim.lengths < 10)
  
  if (length(trim.lengths.min) > 0) {
    my.trimmed.list <- my.trimmed.list[-trim.lengths.min]
  }
  
  trimmed.names <- names(my.trimmed.list)
  
  for (idx in 1:length(my.trimmed.list)) {
    
    p.nodes <- my.trimmed.list[[idx]]$p.nodes
    p.info <- my.trimmed.list[[idx]]$p.info
    
    final.nodes <- as.matrix(p.nodes[, 1:2])
    degrees <- as.numeric(p.info$Freq)
    names(degrees) <- p.info$name
    meta <- data.frame(p.info$group, degrees, p.info$name, ind=1:nrow(p.info))
    my.order <- as.integer(meta[order(meta$degrees, decreasing = TRUE),]$ind)
    
    tiff(paste0(trimmed.names[[idx]],"_TopInteracions.tiff"), width = 14, height = 10, units = 'in', res = 600)
    
    #cex.nodes = 2^(log2(abs(p.info$logFC))/10)
    
    arcplot(final.nodes, ordering=my.order, labels=p.info$name, cex.labels=0.6,
            show.nodes=TRUE, col.nodes=p.info$calfb, bg.nodes=p.info$calfb,
            cex.nodes = (log(degrees)/2)+0.2, pch.nodes=21,
            lwd.nodes = 2, line=-0.5,
            col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = log2(p.nodes$score/100)*1.5)
    dev.off()
  }
}



# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Load sets of differentially abundant proteins

H_L <- read.table("my.wd/HER2_LumA_DA_corrected_pool.txt", header = TRUE)
H_T <- read.table("my.wd/HER2_TNBC_DA_corrected_pool.txt", header = TRUE)
L_T <- read.table("my.wd/LumA_TNBC_DA_corrected_pool.txt", header = TRUE)

ER <- read.table("my.wd/ERp_ERn_DA_corrected_pool.txt", header = TRUE)[-38,]
PgR <- read.table("my.wd/PGRp_PGRn_DA_corrected_pool.txt", header = TRUE)
Her2 <- read.table("my.wd/H23_H01_corrected_pool.txt", header = TRUE)[,-3]


GR <- read.table("my.wd/GR_High_Low.txt", header = TRUE)[,-3]
TILs <- read.table("my.wd/TILs_High_Low_DA_corrected_pool.txt", header = TRUE)[-6,]


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Bind sets together

DA.df <- rbind(H_L, H_T, L_T)
DA.df$comparison <- c(rep("Her2 vs Lum", nrow(H_L)),rep("Her2 vs TNBC", nrow(H_T)), rep("Lum vs TNBC", nrow(L_T)))

DA.df <- ER
DA.df$comparison <- rep("ER+ vs ER-", nrow(ER))

DA.df <- rbind(TILs, GR)
DA.df$comparison <- c(rep("TILs High vs TILs Low", nrow(TILs)), rep("Grade High vs Grade Low", nrow(GR)))


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Generate protein-protein networks and node information

PPIs <- GetDeaPPInt(stringDB, DA.df)
NoIf <- NodeInfo(PPIs)


DA.df <- rbind(H_L, H_T, L_T)
DA.df$comparison <- c(rep("Her2 vs Lum", nrow(H_L)),rep("Her2 vs TNBC", nrow(H_T)), rep("Lum vs TNBC", nrow(L_T)))
