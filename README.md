# Leveraging three-dimensional chromatin architecture for effective reconstruction of enhancer-target gene regulatory interactions

The 3D-ETG repository contains all the data and source scripts required to reproduce the results presented in [*Salviato et al. (2021)*](https://doi.org/10.1101/2021.03.01.432687). A schematic workflow of the methodological framework is available at Figure 1 of the paper.

The complete list of candidate enhancer-promoter pairs (annotated with the HC score, corrected and uncorrected p-values, validations according to multiple reference datasets) is available in the [Result](https://github.com/ElisaSalviato/3D-ETG/tree/main/Results) folder and consultable via [shiny app](https://bioinformatics.ifom.eu/3D-ETG).


## /Data

The data required for the analysis are:
1.	Matrices with signals for enhancers and promoters to quantify the synchronized activity;
2.	Genomic regions that describe domains called by a TADs caller at multiple level of resolution. 


### 1. /Data/Roadmap/Matrix/
The [Data/Roadmap/Matrix](https://github.com/ElisaSalviato/3D-ETG/tree/main/Data/Roadmap/Matrix) folder contains the original functional data used in [*Salviato et al. (2021)*](https://doi.org/10.1101/2021.03.01.432687), from Epigenome Roadmaps for the reference set of enhancers and promoters. Each line contains the raw maximum activity signal measured across cell and tissue types within the genomic region of enhancers or promoters, according to the name file. 

At the end of each matrix the enhancer or promoter region must be specified. Namely, the following columns are mandatory:
-	`chr`, `start`, `end`, for enhancers;
-	`chr`, `start`, `end`, `symbol`, `strand`: for promoters.

Files must be named as `activity_type_key.tsv`, where:
-	`activity`: the type of data used to measure the activity of enhancers and promoters (e.g., DNase, H3K4me3, H3K27ac).
-	`type`: enhancer or promoter;
-	`key`: other info to identify the files.

Customize the script **ComputeCanonicalCorrelationTAD.R**:
1.	*line 26-27*: specify activity and key for enhancers.
2.	*line 29-30*: specify activity and key for promoters.

Note: the matrices will be split by chromosome, pseudo-log2 transformed and normalized by quantile. Change or comment *line:115-126* to customize this option. 



