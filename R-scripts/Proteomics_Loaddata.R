#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Copyright (C) 2018, Thilde Bagger Terkelsen <thilde@cancer.dk> <thildebate@gmail.com>
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

my.wd = ""

# ---------------------------------------------------------------------------------------
# LOADING DATA
# ---------------------------------------------------------------------------------------

setwd(paste0(my.wd,"/Data"))


# ---------------------------------------------------------------------------------------
# LOODING POOLED DATA (TIF, NIF, FIF):
# ---------------------------------------------------------------------------------------

tif_pooled <- read.xlsx("pooled_TIF.xlsx", sheet = 1, colNames = TRUE, rowNames = FALSE)
nif_pooled <- read.xlsx("pooled_NIF.xlsx", sheet = 1, colNames = TRUE, rowNames = FALSE)
fif_pooled <- read.xlsx("pooled_FIF.xlsx", sheet = 1, colNames = TRUE, rowNames = FALSE)

# ORDERING PROTEIN ACCESSION NUMBER
tif_pooled_id <- as.character(sort(tif_pooled$Accession))
nif_pooled_id <- as.character(sort(nif_pooled$Accession))
fif_pooled_id <- as.character(sort(fif_pooled$Accession))



# ---------------------------------------------------------------------------------------
# LOADING PROTEOMICS TIF DATA
# ---------------------------------------------------------------------------------------

tifdata <- read.xlsx("TIFproteomics.xlsx", sheet = 2, startRow = 4, colNames = FALSE, rowNames = FALSE)
tifinfo <- read.xlsx("TIFsampleinfo.xlsx")

# ONLY PATIENT SPECIFIC PROTEOMICS
tifdata_small <- tifdata[,12:46]
colnames(tifdata_small) <- tifinfo$patient
rownames(tifdata_small) <- tifdata$X1


# DROP MISCALSSIFIED
tifdata_small[22] <- NULL
tifinfo <- tifinfo[-22,]

TILs <- as.factor(as.character(ifelse(tifinfo$TILs %in% c("T0", "T1"), "low", "high")))
diag <-  as.factor(as.character(tifinfo$diag))
tp <-  as.factor(as.character(tifinfo$tp))
receptor_status <-  as.factor(as.character(tifinfo$receptor_status))
pool <-  as.factor(as.character(paste0("p",tifinfo$pool)))
Gr <- as.factor(as.character(ifelse(tifinfo$Gr == "3", "high", "low")))
ER <- as.factor(as.character(tifinfo$ER))
PGR <- as.factor(as.character(tifinfo$PGR))
AR <- as.factor(as.character(tifinfo$AR))
HER2 <- as.factor(as.character(tifinfo$HER2))

# NUMBER OF NA's PER ROW - REMOVE THE PROTEINS WITH MORE THAN 12 NA's.
numNAs <- apply(tifdata_small, 1, function(x) sum(is.na(x)))
tifdata_small <- tifdata_small[!(numNAs > 12),]


# ADD PSEUDOCOUND (MEAN ABUNDANCE VALUE)
tif_means <- as.vector(rowMeans(tifdata_small, na.rm = TRUE))

for(i in 1:nrow(tifdata_small)){
  tifdata_small[i, is.na(tifdata_small[i,])] <- tif_means[i]
}


# ---------------------------------------------------------------------------------------
# LOADING PLASMA PROTEOMICS FROM Farrah et. al.
# ---------------------------------------------------------------------------------------

plasma <- read.xlsx("plasma.xlsx", sheet = 2, colNames = TRUE, rowNames = FALSE)
plasma <- sort(plasma$identifier)


# ---------------------------------------------------------------------------------------
# LOADING BREAST CANCER SECRETED PROTEOMICS FROM Mann et. al.
# ---------------------------------------------------------------------------------------

secreted <- read.table("secretome.txt", header=TRUE)
secreted <- unique(sort(as.character(secreted$uniprot)))



# ---------------------------------------------------------------------------------------
# LOADING EXSOSOMES Willms et. al.
# ---------------------------------------------------------------------------------------

exosomes <- as.character(read.table("exosomes.txt")$V1)


# ---------------------------------------------------------------------------------------
# LOADING BREAST CANCER N-GLYCOSYLATED PROTEOMICS FROM Mann et. al.
# ---------------------------------------------------------------------------------------
nglyco <- read.table("N_glyco.txt",header = TRUE)
nglyco <- unique(sort(as.character(nglyco$uniprot)))
nglyco <- nglyco[nglyco %in% secreted]


# --------------------------------------------------------------------------------------
# LOADING PROT2GO.Rdata (Already mapped protein to GO terms dataframe).
# --------------------------------------------------------------------------------------

load("prot2go.RData")
GO_background <- grps
rm(grps)



# ---------------------------------------------------------------------------------------
# LOADING HUMAN MAP DATA FOR CONVERSION OF UNIPROT TO ENTREZ
# ---------------------------------------------------------------------------------------
load("human_map_26Aug15.RData")



# ---------------------------------------------------------------------------------------
# LODING HUGO DATA FOR INFO ON GENES
# ---------------------------------------------------------------------------------------

hugo_original <- read.delim("gene_association.goa_human", header=F)
hugo <- unique(data.frame(hugo_original$V2, hugo_original$V3, hugo_original$V10)) 
colnames(hugo) <- c("Accession", "name", "Description")


# ---------------------------------------------------------------------------------------
# LODING HUMAN KEGGS AND GENESETS
# ---------------------------------------------------------------------------------------

# human KEGG pathways and genesets
load("kegg.hsa.sigmet.gsets.RData")
load("go.hs.gsets.RData")


# ---------------------------------------------------------------------------------------
# LODING DATASET FOR DE COMPARISON
# ---------------------------------------------------------------------------------------

TyanovaGST <- read.table("All_genes_from_Tyanova_Geneset.txt", header = FALSE)
TyanovaDEG <- read.table("DE_gene_from_Tyanova_Geneset.txt", header = FALSE)

# AngeloDEG <- read.table("Angelo_DE_genes.txt", header = TRUE)


# ---------------------------------------------------------------------------------------
# LODING PAM50 genes
# ---------------------------------------------------------------------------------------


pam50 <- c("ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20", "CDC6", "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1", "FGFR4", "FOXA1", "FOXC1", "GPR160", "GRB7", "KIF2C", "KRT14", "KRT17", "KRT5", "MAPT", "MDM2", "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80", "NUF2", "ORC6L", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS", "UBE2C", "UBE2T")


# ---------------------------------------------------------------------------------------
# VECTORS WITH COLOR CODES FOR DIAG AND RECEPTOR STATUS 
# ---------------------------------------------------------------------------------------

colorcode_receptor <- c("#F8766D", "#FD61D1", "#B983FF", "#00BA38", "#00BCD8")
colorcode_diag <- c("#FD61D1","#00B0F6","#00BF7D")
#colorcode_diag_new <- c("#FD61D1","#00B0F6","#00BF7D", "#B983FF")


# myprots <- c("Q8TD06", "P50895", "Q9NYQ6", "Q9BRT3", "P18440", "P78356", "Q15437", "Q9BU02", "Q9NW97", "Q9BZM5")
