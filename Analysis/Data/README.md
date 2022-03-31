# Datasets Description

## Overview
This repository provides the data for the manuscript "Factors determining distributions of rainforest Drosophila shift from interspecific competition to high temperature with decreasing elevation"

We investigated thermal tolerances and interspecific competition as causes of species turnover in the nine most abundant species of Drosophila along elevational gradients in the Australian Wet Tropics. Specifically, we 1) analyzed the distribution patterns of the studies Drosophila species; 2) fitted thermal performance curves; 3) tested the correlation between multiple thermal traits and distribution patterns; 4) fitted the Beverton-Holt model to describe the single-generation intra- and inter-specific competition effect; 5) examined the long-term effect of competition and temperature on the population size of a pair of Drosophila species.


## Authors information
Jinlin Chen (Department of Zoology, University of Oxford)
Owen T. Lewis (Department of Zoology, University of Oxford)

Authors contribution: JC and OTL both contributed to the development of ideas. JC designed and conducted the experimental work. JC analyzed the results and led the writing of the manuscript. OTL contributed to the writing.

Correspondence: Jinlin Chen (jinlinchenn@gmail.com)


## Layout
The repository contains 7 datasets. The Protocol to collect these data were described in detail in the related publication: XXX. 

### **`pupaeSamplingCore.csv`** 
This dataset is to analyze the distribution pattens of the Drosophila species along elevation gradients. Drosophila pupae were sampled using bottle traps baited with fermented banana from 11th March – 12th April 2016 for three sites at elevation of 70m, 350-390m, and 730-880m on the two mountain ranges. Variables that were used in the analysis of this manuscript were described below:
 * **`Transect`**: Two mountain ranges where the sampling was conducted  
 * **`Site`**: The name of the sampling sites. The first character indicate the transect and the following number indicate the elevation of the particular site. For example, "P880" indicates that the sample was collected on Paluma transect on 880 metre elevation. 
 * **`Host`**: The species identity of the sample. 
 
 ### **`longtermClimate_formated.csv`** 
 This dataset is the hourly temperautre and humidity of survey sites on Paluma and Kirama, Queensland, Australia, recorded from Apri 2016 to March 2017. Variables that were used in the analysis of this manuscript were described below:
  * **`Site2016`**: It indicates the mountain range (Kirama or Paluma) and the elevation of the site. 
  * **`Celsius`**: The recorded temperature.
  * **`year_month`**: The year and month that the climate data was recorded.
  * **`justTime.c`**: The corrected time of the day that the climate data was recorded.
  * **`day`**: The date of the corresponding month that the climate data was recorded.

### **`TPC_data_without na.csv`** 
This dataset is to analyze the reproductive thermal performance of the Drosophila species in the laboratory. We exposed flies to temperatures ranging from 14°C to 32°C and measured how their reproductive success changes with temperature. Variables were described below:
 * **`sp`**: The name of the species tested  
 * **`round`**: The biological repeat block. Two blocks in total, namely "1" and "2". 
 * **`pos`**: The position of the tested Drosophila vial in the vial racks (10 rows per rack)
 * **`period`**: The period that the parent flies were exposed to temperature treatment. Day count starts from the first day that all parent flies are exposed to temperature treatment. 
 * **`rep`**: Replication. 4 replicates per round.
 * **`temp`**: The set temperature of the water bath.
 * **`temp.corr`**: The realized average temperature of the water bath that carry to specific Drosophila vials. 
 * **`EL`**: Whether (Y) or not (N) that eggs or larvae were observed by eye inspection.
 * **`pupa`**: Whether (Y) or not (N) that pupae were observed by eye inspection.
 * **`adultF`**: The number of adult female offspring produced in the particular Drosophila vial
 * **`adultM`**: The number of adult male offspring produced in the particular Drosophila vial

### **`cold_tolerance_data.csv`** 
This dataset documents the knockdown and recovery time by extreme coldness of the Drosophila species in the laboratory. Resistance to extreme cold temperature was measured as knockdown time for each individual at 5°C and the time for recovery of mobility after a 30-minute exposure to 5°C. Variables were described below:
 * **`treatment`**: The experimental temperature that the flies were exposed to
 * **`species`**: The name of the species tested  
 * **`distribution`**: The distribution type of the tested species
 * **`gender`**: The gender of the adult flies
 * **`age`**: The age of the adult flies, starting from the emergence
 * **`round`**: The biological repeat block. Three blocks in total, namely "A", "B" and "C". 
 * **`position`**: The position of the tested Drosophila vial on the observation rack (3*3 grid)
 * **`order`**: The order of the 7 tubes in one cell on the observation rack
 * **`rep`**: Replication. 7 replicates per round.
 * **`kd.t`**: Knockdown time (min)
 * **`rc.t`**: Recovery time (min)
 
### **`hot_tolerance_data.csv`** 
This dataset documents the knockdown by extreme heat of the Drosophila species in the laboratory. Resistance to extreme cold temperature was measured as knockdown time for each individual at 40°C. Variables were described below:
 * **`treatment`**: The experimental temperature that the flies were exposed to
 * **`species`**: The name of the species tested  
 * **`distribution`**: The distribution type of the tested species
 * **`gender`**: The gender of the adult flies
 * **`age`**: The age of the adult flies, starting from the emergence
 * **`round`**: The biological repeat block. Three blocks in total, namely "A", "B" and "C". 
 * **`position`**: The position of the tested Drosophila vial on the observation rack (3*3 grid)
 * **`order`**: The order of the 7 tubes in one cell on the observation rack
 * **`rep`**: Replication. 7 replicates per round.
 * **`kd.t`**: Knockdown time (min)
 
### **`sixPairs_cold.csv`** 
This dataset documents the short-term competition outcome in cold treatment. Specific numbers of one or two Drosophila species were placed in one vial in a experimental temperature. The outcome were measured by the number of offspring of each species. Variables were described below:
 * **`tubeID`**: A unique index for each tube to measure competition
 * **`type`**: Whether the vial contains one (intra) or two (inter) Drosophila species 
 * **`temp`**: The experimental temperature that the flies were exposed to
 * **`block`**: Two block from the SEP experiment, and three blocks from the DEC experiment. 
 * **`replicates`**: 3-5 replicates per block.
 * **`tray`**: The tray that held the experimental tube.
 * **`pair`**: The pair or single species that were in the tube.
 * **`i`**: The name of the focal species
 * **`i.ind`**: The species index of the focal species
 * **`FocalSp_Den`**: The starting density of the focal species
 * **`j`**: The name of the competing species. If the tube contains only one species, the name is NaN
 * **`j.ind`**: The species index of the competing species. If the tube contains only one species, the index is 0
 * **`OtherSp_Den`**: The starting density of the competing species
 * **`densityT`**: The density combination of the two species
 * **`obs_count`**: The number of offspring after one generation of the focal species

### **`sixPairs_hot.csv`** 
This dataset documents the short-term competition outcome in hot treatment. A lowland species, D. pandora, and an upland species, D. pallidifrons, were reared in monoculture and mixed-culture environments for multiple generations in lowland and upland temperature regimes. The outcome were measured by the number of offspring of each species. Variables were described below:
 * **`tubeID`**: A unique index for each tube to measure competition
 * **`type`**: Whether the vial contains one (intra) or two (inter) Drosophila species 
 * **`temp`**: The experimental temperature that the flies were exposed to
 * **`block`**: Two block from the SEP experiment, and three blocks from the DEC experiment. 
 * **`replicates`**: 3-5 replicates per block.
 * **`tray`**: The tray that held the experimental tube.
 * **`pair`**: The pair or single species that were in the tube.
 * **`i`**: The name of the focal species
 * **`i.ind`**: The species index of the focal species
 * **`FocalSp_Den`**: The starting density of the focal species
 * **`j`**: The name of the competing species. If the tube contains only one species, the name is NaN
 * **`j.ind`**: The species index of the competing species. If the tube contains only one species, the index is 0
 * **`OtherSp_Den`**: The starting density of the competing species
 * **`densityT`**: The density combination of the two species
 * **`obs_count`**: The number of offspring after one generation of the focal species

### **`long_term_competition_PAN-PAL.csv`** 
This dataset documents the long-term competition outcome in cold and hot temperautre treatments. Variables were described below:
 * **`vialID`**: A unique index for each tube to measure competition
 * **`week`**: The week of when the census was taken since the start of the experiment
 * **`treatment`**: The experimental temperature that the flies were exposed to
 * **`block`**: Two block starting one day apart: "A" and "B".
 * **`intra`**: Whether (TRUE) or not (FALSE) that this is the monoculture of one species.
 * **`species`**: The identity of the species that are examined in this tube
 * **`rep`**: 4 replicates for monoculture and 8 replicated for mix-species culture in total (two blocks combined) 
 * **`count.type`**: The type of counted flies. "survial" indicates these flies are the surviving flies from the most current bottle. "reproduction" indicates these flies are the newly emerged adults from the other old bottles.
 * **`count`**: The number of flies of the focal species defined by **`species`**.


# License Information
To the extent possible under law, *Jinlin Chen* has waived all copyright to use the included code as a template for related analysis. Please cite the archived or published version of the paper if adapting the code for your research. Please contact *Jinlin Chen* for potential collaboration if you want to use any of the datasets involved in this project. 

Copyright (c) 2022 JINLIN CHEN
