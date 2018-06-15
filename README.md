# Proteomics-TIF

Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark.

Repository associated with the publication:

###

corresponding author: Elena Papaleo, elenap@cancer.dk

r-scripts were created by: Thilde Bagger Terkelsen, thilde@cancer.dk

This repository contains proteomics data from tumour interstitial fluids obtained from a cohort of patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used for analysis in relation to the publication.

The repository contains two folders:

(1) Proteomics data and patient metadata. These are the data used as the starting point for our analyses.

(2) A collection of R-scripts that recapitulate our work:
                                

Requirements:

R version 3.3.1 or higher
Rstudio version 1.1.383 or higher        
Bioconductor version 3.6 or higher	

Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages for manual installation:

CRAN:

openxlsx
ggplot2
dendextend
heatmap.plus
reshape
gdata
plyr
data.table
RColorBrewer
squash
survminer
car
scales
statmod

Bioconductor:

limma
sva
survcomp               
