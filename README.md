# Predicting crop yield
### Lutz et al
This is the home of all the bioinformatic code used in our paper: LINK  
Calculations were performed at sciCORE (http://scicore.unibas.ch/) scientific computing center at University of Basel.

#### 01: Demultiplexing
This folder contains all files used for converting raw data to CCS sequences and demultiplexing.

#### 02: ASV dadapipe
This folder contains all code to run the dadapipe. With **dada2** sequences were oriented, quality filtered, truncated, dereplicated and denoised. A count table was created and taxonomy was assigned. **DECIPHER** was used to cluster sequences by similarity. The lock-file contains the state of our project library and can be reproduced with the **renv** package.

#### 03: ASV tables
This folder contains all output files created by the dadapipe. These are count tables, sequences and taxonomy of ASVs with different similarity thresholds and files to check quality. 


#### 04: MGR prediction
This folder contains the statistical analysis for the MGR prediction.
