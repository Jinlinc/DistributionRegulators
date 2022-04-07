# Abiotic and biotic regulators of species distribution

## Overview
This repository provides the data, analysis code for the manuscript "Factors determining distributions of rainforest Drosophila shift from interspecific competition to high temperature with decreasing elevation"

We investigated thermal tolerances and interspecific competition as causes of species turnover in the nine most abundant species of Drosophila along elevational gradients in the Australian Wet Tropics. Specifically, we 1) analyzed the distribution patterns of the studies Drosophila species; 2) fitted thermal performance curves; 3) tested the correlation between multiple thermal traits and distribution patterns; 4) fitted the Beverton-Holt model to describe the single-generation intra- and inter-specific competition effect; 5) examined the long-term effect of competition and temperature on the population size of a pair of Drosophila species.


## Authors information
Jinlin Chen (Department of Zoology, University of Oxford)
Owen T. Lewis (Department of Zoology, University of Oxford)

Authors contribution: JC and OTL both contributed to the development of ideas. JC designed and conducted the experimental work. JC analyzed the results and led the writing of the manuscript. OTL contributed to the writing.

Correspondence: Jinlin Chen (chen.jin.lin@hotmail.com)


## Layout

The repository is split into 3 main directories, many of which have subdirectories. Each section is described below. 

### **`Analysis`** 
Where all of the *executed* analysis lives. This includes two scripts, Data, stanFits and StanModels. 
 * **`Data`**: all the data for this manuscript. See more details in the **`README.md`** file within directory. 
 * **`Main_results_analysis.R`**: all code to analyze the data and generate tables, figures, statisticla results in the main text and the SI. 
 * **`NegBFunction_1treatment_6Pairs_5block`**: A script to generate the stanModel for the pairwise competition analysis.
 * **`StanFits`**: Fitted stan models. The actually fitted models were not included as they are two large for github repository.
 * **`StanModels`**: The ".stan" files used to fit the TPC and competition models.

### **`Results`** 
All table and figures in the main text and most of supplementary materials generated from the code are included in this section. 
 * **`Main`**: Two tables and five figures from the main text.
 * **`Supp`**: Two tables, eleven figures from SI and one phylogenetic tree file.

### **`Writing`** 
The old and current versions of the munuscript and related documents for submission to scientific journals. The folder is still empty now.


## Version control
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Catalina 10.15.7

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] brms_2.14.4          Rcpp_1.0.6           ape_5.4-1            bayesplot_1.8.0     
 [5] rstan_2.21.1         StanHeaders_2.21.0-7 cowplot_1.1.1        ggrepel_0.9.1       
 [9] lme4_1.1-26          Matrix_1.2-18        forcats_0.5.1        stringr_1.4.0       
[13] dplyr_1.0.4          purrr_0.3.4          readr_1.4.0          tidyr_1.1.2         
[17] tibble_3.0.6         ggplot2_3.3.3        tidyverse_1.3.0     

loaded via a namespace (and not attached):
  [1] cubature_2.0.4.1     minqa_1.2.4          colorspace_2.0-0     ellipsis_0.3.1      
  [5] ggridges_0.5.3       rsconnect_0.8.16     markdown_1.1         corpcor_1.6.9       
  [9] base64enc_0.1-3      fs_1.5.0             rstudioapi_0.13      DT_0.17             
 [13] mvtnorm_1.1-1        lubridate_1.7.9.2    xml2_1.3.2           bridgesampling_1.0-0
 [17] codetools_0.2-16     splines_4.0.3        shinythemes_1.2.0    projpred_2.0.2      
 [21] jsonlite_1.7.2       nloptr_1.2.2.2       broom_0.7.4          dbplyr_2.1.0        
 [25] shiny_1.6.0          compiler_4.0.3       httr_1.4.2           backports_1.2.1     
 [29] assertthat_0.2.1     fastmap_1.1.0        cli_2.3.0            later_1.1.0.1       
 [33] htmltools_0.5.1.1    prettyunits_1.1.1    tools_4.0.3          igraph_1.2.6        
 [37] coda_0.19-4          gtable_0.3.0         glue_1.4.2           reshape2_1.4.4      
 [41] V8_3.4.0             cellranger_1.1.0     vctrs_0.3.6          nlme_3.1-149        
 [45] crosstalk_1.1.1      tensorA_0.36.2       ps_1.5.0             rvest_0.3.6         
 [49] miniUI_0.1.1.1       mime_0.9             lifecycle_0.2.0      gtools_3.8.2        
 [53] statmod_1.4.35       MASS_7.3-53          zoo_1.8-8            scales_1.1.1        
 [57] colourpicker_1.1.0   hms_1.0.0            promises_1.1.1       Brobdingnag_1.2-6   
 [61] parallel_4.0.3       inline_0.3.17        shinystan_2.5.0      gamm4_0.2-6         
 [65] yaml_2.2.1           curl_4.3             gridExtra_2.3        loo_2.4.1           
 [69] stringi_1.5.3        dygraphs_1.1.1.6     boot_1.3-25          pkgbuild_1.2.0      
 [73] rlang_0.4.10         pkgconfig_2.0.3      matrixStats_0.58.0   lattice_0.20-41     
 [77] rstantools_2.1.1     htmlwidgets_1.5.3    processx_3.4.5       tidyselect_1.1.0    
 [81] plyr_1.8.6           magrittr_2.0.1       R6_2.5.0             generics_0.1.0      
 [85] DBI_1.1.1            pillar_1.4.7         haven_2.3.1          withr_2.4.1         
 [89] mgcv_1.8-33          xts_0.12.1           abind_1.4-5          modelr_0.1.8        
 [93] crayon_1.4.1         grid_4.0.3           readxl_1.3.1         MCMCglmm_2.32       
 [97] callr_3.5.1          threejs_0.3.3        reprex_1.0.0         digest_0.6.27       
[101] xtable_1.8-4         httpuv_1.5.5         RcppParallel_5.0.2   stats4_4.0.3        
[105] munsell_0.5.0        shinyjs_2.0.0    


# License Information
To the extent possible under law, *Jinlin Chen* has waived all copyright to use the included code as a template for related analysis. Please cite the archived or published version of the paper if adapting the code for your research. Please contact *Jinlin Chen* for potential collaboration if you want to use any of the datasets involved in this project. 

Copyright (c) 2022 JINLIN CHEN
