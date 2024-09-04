# DNA methylation Analysis
This project is part of the DNA methylation pipeline developed at Center of Epigenetics Research (CER), MSKCC.

It will process bismark output, perform exploratory analyses, and conduct differential methylation. Diagnositc plots and a bsseq object will be generated from the exploratory analyses in "explore" mode. Differential DMCs and DMRs will be identified using 'DSS' and annotated with 'annotatr'. The resulted tables can be used for further downstream analyses.  

The code has been executed with R4.0.4 and depends on the following packages:
```
library(data.table)
library(BiocParallel)
library(bsseq)
library(aod)
library(DelayedMatrixStats)
library(ggplot2)
library(dendextend)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggcorrplot)
library(rafalib)
library(factoextra)
library(ggpubr)
library(httr)
set_config(config(ssl_verifypeer = 0L)) # solve the biomart connection problem that may happen sometimes
library(annotatr)
```
