# Proteomics-TIF

Computational Biology Laboratory, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark.

Repository associated with the publication:


High-throughput proteomics of breast cancer interstitial fluid: identification of tumor subtype-specific serologically relevant biomarkers.
Terkelsen T, Pernemalm M, Gromov P, Børresen-Dale AL, Krogh A, Haakensen VD, Lethiö J, Papaleo E, Gromova I.
Mol Oncol. 2021 Feb;15(2):429-461. doi: 10.1002/1878-0261.12850.

Reference contact for repository: Elena Papaleo, elenap-at-cancer.dk

The program files (R scripts - files with extension .R) in this repository are licensed under the terms of the GNU General Public License (see LICENSE file). You can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. These programs are distributed in the hope that they will be useful,but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
![image](https://user-images.githubusercontent.com/12696167/156525794-07795c9f-e0c6-4a61-bbb0-ad9826626aa5.png)



This repository contains proteomics data from tumour interstitial fluids obtained from a cohort of patients with breast cancer. The repository was made with intent of openly sharing both data and R-scripts used for analysis in relation to the publication.

The repository contains X folders:

    (1) Data: Proteomics data and patient metadata. These are the data used as the starting point for our analyses.
    (2) Backgrounds_and_Databases: Databases and files used for the analysis, KEGG pathways maps, i2d database etc. 
    (3) R-scripts: A collection of R-scripts that recapitulate the work.
                                

Requirements:

    R version 4.0.0 or higher
    Rstudio version 1.2.5019 or higher        	

Although R-packages should automatically be installed and errors raised if they cannot be, we here provide the user with the list of required packages for manual installation:

CRAN:

    pamr
    openxlsx
    ggplot2
    dendextend
    heatmap.plus
    reshape
    VennDiagram
    RColorBrewer
    stringr
    caTools
    plyr
    gplots
    corrplot
    fitdistrplus
    randomForest
    glmnet
    e1071
    caret
    varSelRF

Bioconductor:

    limma
    sva
    topGO
    GOSim
    gdata
    UpSetR
    biomaRt
    WGCNA
    viridis
    gprofiler2

