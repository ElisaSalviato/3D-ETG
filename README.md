# Leveraging three-dimensional chromatin architecture for effective reconstruction of enhancer-target gene regulatory interactions

The 3D-ETG repository contains all the data and source scripts required to reproduce the results presented in [*Salviato et al. (2021)*](https://doi.org/10.1101/2021.03.01.432687). 
<p align="center">
  <img width="30%" src="https://github.com/ElisaSalviato/3D-ETG/blob/main/Images/Figure1_revised_withcaption_small.png">
</p>

The complete list of candidate enhancer-promoter pairs (annotated with the HC score, corrected and uncorrected p-values, validations according to multiple reference datasets) is available in the [Result](https://github.com/ElisaSalviato/3D-ETG/tree/main/Results) folder and consultable via [shiny app](https://bioinformatics.ifom.eu/3D-ETG).


## /Data

The data required for the analysis are:
1.	Matrices with signals for enhancers and promoters to quantify the synchronized activity;
2.	Genomic regions that describe domains called by a TADs caller at multiple level of resolution. 


### 1. Data/Roadmap/Matrix/
The [Data/Roadmap/Matrix](https://github.com/ElisaSalviato/3D-ETG/tree/main/Data/Roadmap/Matrix) folder contains data from Epigenome Roadmaps for the reference set of enhancers and promoters. 

Each line contains the raw maximum activity signal measured across cell and tissue types (columns) within the genomic region of enhancers or promoters, according to the name file. At the end of each matrix the enhancer or promoter region must be specified. Namely, the following columns are mandatory:
-	`chr`, `start`, `end`, for enhancers;
-	`chr`, `start`, `end`, `symbol`, `strand`: for promoters.

Files must be named as `activity_type_key.tsv`, where:
-	`activity`: the type of data used to measure the activity of enhancers and promoters (e.g., DNase, H3K4me3, H3K27ac).
-	`type`: enhancer or promoter;
-	`key`: other info to identify the files.

Customize the script [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R):
-	*line 26-27*: specify `key.enhancer` and `activity.enhancer`, according to `key` and `activity` information in the matrix file names for enhancers;
-	*line 29-30*: specify `key.promoter` and `activity.promoter`, according to `key` and `activity` information in the matrix file names for promoter.

**Note**: the matrices will be split by chromosome, pseudo-log2 transformed and normalized by quantile. Change or comment *line:115-126* in [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R) to customize this option. 


### 2. Data/TAD/
The [Data/TAD/](https://github.com/ElisaSalviato/3D-ETG/tree/main/Data/TAD) folder contains Topological Associating Domains (TADs) called for eleven Hi-C datasets covering nine different cell and tissue types, using the Local Score Differentiator ([LSD](https://www.bioconductor.org/packages/release/bioc/vignettes/HiCBricks/inst/doc/IntroductionToHiCBricks.html#call-topologically-associated-domains-with-local-score-differentiator-lsd)) algorithm. TADs are used to calculate the Hierarchical Contact (HC) score, a score proportional to the likelihood of enhancer-promoter (EP) pairs co-localization.

Each file contains TADs called for one specific chromosome, organized in a list object:
- level 1: list of Hi-C matrix resolutions (`10000` bp);
- level 2: list of di-windows (`5`,`10`,`20`,`50`), also called TADs hierarchy (from lower, to higher;
- level 3: GenomicRanges object contains domains coordinates.

Files must be names as `TAD_method_chr_binsize_info.RData`, where:
- `TAD`: is a mandatory flag to reconize files containing domains;
- `method`: the algorithm used to call TADs;
- `chr`: the reference chromosome;
- `binsize`: the resolution of the Hi-C matrix;
- `info`: other info to identify the Hi-C experiment.

The preferable TADs caller can be used, as long as the above described lists and file names structure are preserved. Algorithm that are able to call different hierarchy of TADs are advisable.

**Note**: the level 2 names of the di-windows list will be converted as numeric and used to compute a scaling factor for the HC score: be sure to provide valid names for the conversion from string to numeric. To avoid that the scaling factor is computed set `use.resolution=FALSE` at line:100 in [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R).





