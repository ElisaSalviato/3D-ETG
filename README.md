# Leveraging three-dimensional chromatin architecture for effective reconstruction of enhancer-target gene regulatory interactions

The 3D-ETG repository contains all the data and source scripts required to reproduce the results presented in [*Salviato et al. (2021)*](https://doi.org/10.1101/2021.03.01.432687). The complete list of candidate enhancer-promoter pairs (annotated with the HC score, corrected and uncorrected p-values, validations according to multiple reference datasets) is available in the [Result](https://github.com/ElisaSalviato/3D-ETG/tree/main/Results) folder and consultable via [shiny app](https://bioinformatics.ifom.eu/3D-ETG).

3D-ETG is a general framework for the definition of enhancer-target gene (ETG) pairs leveraging the current biological knowledge on chromatin 3D architecture and integrating heterogeneous functional genomics data into a rigorous statistical framework. 
Its three key features are:
1. **Statistical framework for quantifying enhancer-promoter pairs synchronization** ([Figure 1A](https://github.com/ElisaSalviato/3D-ETG/blob/main/Images/Figure1_revised_withcaption_small.png)). The method is flexible in terms of input, as it starts from user-defined sets of i) enhancer and promoter regions and ii) functional genomics data to quantify their activity . This flexibility is ensured by the use of Canonical-Correlation Analysis (CCA) to quantify the synchronization of enhancer-promoter (EP) pairs activity across cell types. Moreover, it is designed to leverage multiple types of functional genomics data, also accounting for the correlation within sets of features.
2. **Hierarchical Contact (HC) score** ([Figure 1B](https://github.com/ElisaSalviato/3D-ETG/blob/main/Images/Figure1_revised_withcaption_small.png)). It incorporates chromatin architecture as experimentally measured by Hi-C, to compute the HC score accounting for ETG pairs 3D co-localization. Differently from previous methods, we leverage biological knowledge on TADs multi-scale hierarchical organization and their conservation across cell types.
3. **Chromatin 3D architecture and functional genomics data integration** ([Figure 1C](https://github.com/ElisaSalviato/3D-ETG/blob/main/Images/Figure1_revised_withcaption_small.png)). The information on chromatin 3D architecture is used to increase the statistical power to detect ETG pairs synchronization, while controlling false discoveries. This is the first time that chromatin 3D architecture is directly integrated as side information in the statistical model for defining ETG pairs.

<p align="center">
  <img width="60%" src="https://github.com/ElisaSalviato/3D-ETG/blob/main/Images/Figure1_revised_withcaption_small.png">
</p>


## /Data: input data

The data required for the analysis are:
1.	two sets of matrices containing activity status for enhancers and promoters;
2.	genomic regions that describe hierarchy of chromatin structural domains derived from a set of Hi-C datasets. 


### 1. Data/Roadmap/Matrix/
The [Data/Roadmap/Matrix](https://github.com/ElisaSalviato/3D-ETG/tree/main/Data/Roadmap/Matrix) folder contains DNase-seq and ChIP-seq enrichment profiles quantifing the activity of our reference set of enhancers and promoters. Namely, we downloaded H3K27ac, H3K4me3, and DNase-seq consolidated fold-change enrichment signal tracks (bigwig format) from the Roadmap Epigenomic consortium [web portal](https://egg2.wustl.edu/roadmap/data/byFileType/signal/consolidated/macs2signal/foldChange/) for all the cell and tissue types for which all the three epigenetics marks were available. 

Each line of each matrix contains the maximum signal measured across cell and tissue types (columns) within the genomic region of enhancers or promoters, according to the name file. The maximum signal has been computed using [rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html) R package.
At the end of each matrix the enhancer or promoter region must be specified. Namely, the following columns are mandatory:
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
The [Data/TAD/](https://github.com/ElisaSalviato/3D-ETG/tree/main/Data/TAD) folder contains Topological Associating Domains (TADs) called  using the Local Score Differentiator (LSD) TAD borders calling procedure ([LSD](https://www.bioconductor.org/packages/release/bioc/vignettes/HiCBricks/inst/doc/IntroductionToHiCBricks.html#call-topologically-associated-domains-with-local-score-differentiator-lsd)) algorithm, for eleven Hi-C datasets covering different cell lines and primary tissues from a compendium of public datasets. 

We defined TADs as regions between two consecutive domain boundaries. LSD is based on the directionality index (DI) score originally proposed by [*Dixon et al. (2012)*](https://pubmed.ncbi.nlm.nih.gov/22495300/). Among the user defined parameters in this algorithm, the DI-window (i.e., the number of up-stream and down-stream bins over which the DI score is computed) influences the scale of the TAD domains that are identified: the larger the DI-window, the larger the average resulting TAD size. TADs are used to calculate the Hierarchical Contact (HC) score, a score proportional to the likelihood of enhancer-promoter (EP) pairs co-localization.

Each file contains TADs called for one specific chromosome, organized in a list object:
- level 1: list of Hi-C matrix resolutions (`10000` bp);
- level 2: list of di-windows (`5`,`10`,`20`,`50`), also called TADs hierarchy (from lower to higher);
- level 3: a [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) object containing domains coordinates.

Files must be named as `TAD_method_chr_binsize_info.RData`, where:
- `TAD`: is a mandatory flag to reconize files containing domains;
- `method`: the algorithm used to call TADs;
- `chr`: the reference chromosome;
- `binsize`: the resolution of the Hi-C matrix;
- `info`: other info to identify the Hi-C experiment.

The preferable TADs caller can be used, as long as the above described lists and file names structure are preserved. Algorithm that are able to call different hierarchy of TADs are advisable.

**Note**: the level 2 names of the di-windows list will be converted as numeric and used to compute a scaling factor for the HC score: be sure to provide valid names for the conversion from string to numeric. To avoid that the scaling factor is computed set `use.resolution=FALSE` at *line:100* in [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R).


## /Rscript: run the analysis
The [Rscript](https://github.com/ElisaSalviato/3D-ETG/tree/main/Rscript) folder contais the scripts required to perform the 3D-ETG analysis. Namely:
1. [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R): it performs the CCA (as implemented in the [ccaPP](https://cran.rstudio.com/web/packages/ccaPP/index.html) R package) to quantify the strenght of coordinated activity for each EP pair and calculate the HC score. The function will return a [data.table](https://cran.r-project.org/web/packages/data.table/index.html) object (called `Cca.chr.res`) for each chromosome, that will be automatically saved in the [Results/3D-ETG](https://github.com/ElisaSalviato/3D-ETG/tree/main/Results/3D-ETG) folder (RData format). Together with the result table, to [GenomicRanges](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html) objects with enhancers (`E.chr`) and promoters (`P.chr`) coordinates of each pair are returned.
2. [AdjustPvaluesCca_20200305.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/AdjustPvaluesCca_20200305.R): it estimates Bayes-optimal p-value rejection threshold based on the 3D co-localization information encoded in the HC score, as implemented in the [adaptMT](https://cran.r-project.org/web/packages/adaptMT/index.html) R package. It allows to prioritize hypothesis that are more likely to be false. The function will read the [data.table](https://cran.r-project.org/web/packages/data.table/index.html) objects saved in the [Results/3D-ETG](https://github.com/ElisaSalviato/3D-ETG/tree/main/Results/3D-ETG) folder and will updat them with additional columns reporting the adjusted p-values. 
3. [UtilityFunction_20200305.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/UtilityFunction_20200305.R): utility functions that are sourced by the two main functions described above.

A representative example of the final expected output for chromsome 19 is provided ([ETG_Cca_LSD_chr19_10000.Rdata](https://github.com/ElisaSalviato/3D-ETG/blob/main/Results/3D-ETG/ETG_Cca_LSD_chr19_10000.Rdata)).

**Note**: the CCA analysis can be performed for EP pairs localized within the same TAD in at least one of the hierarchy levels in at least one Hi-C dataset (i.e., HC score ≥ 1). To ensure the robustness of the downstream results, the script automatically discarded poorly-supported EP pairs (HC score ≤ `th.weight`, where `th.weight` is equal to the number of provided Hi-C datasets), as they may be the consequence of noise in the data depending on technical variables (e.g. coverage). To avoid this filter set `th.weight=0` at *line:107* in [ComputeCanonicalCorrelationTAD.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/Rscript/ComputeCanonicalCorrelationTAD_20200305.R).


## /Results
The [Formatted_Candidate_ETG_20200305_TAD10Kb.tsv.gz](https://github.com/ElisaSalviato/3D-ETG/blob/main/Results/Formatted_Candidate_ETG_20200305_TAD10Kb.tsv.gz) file provide the complete list of candidate enhancer-target pairs (ETG) inferred in [*Salviato et al. (2021)*](https://doi.org/10.1101/2021.03.01.432687). The results are interactively consultable via shiny app at this [link](https://bioinformatics.ifom.eu/3D-ETG) or running the [app.R](https://github.com/ElisaSalviato/3D-ETG/blob/main/app.R) script using:

```
shiny::runApp("3D-ETG/")
```

Each line corresponds to an ETG pair, where:
- `chr`, `start`, `end`: enhancer region (hg19). Promoter proximal elements are not considered (3.5 kb upstream and 1.5 kb downstream of coding and non-coding TSS).
symbol: target-gene symbol of protein coding genes, based on RefSeq annotations in UCSC. Target genes with overlapping promoter regions are merged as a single pairs (separated by semicolon).
- `distance`: distance between mid-points of gene promoter and enhancer regions. Promoters are defined as 1.5kb upstream and 0.5kb downstream regions of TSS of coding genes.
- `cca`: first canonical correlation coefficient. It is based on the enrichment of DNase-seq and H3K27ac for enhancers and DNase-seq, H3K27ac and H3K4me3 for promoters.
- `HC`: Hierarchical Contact (HC) score. HC score is proportional to the likelihood of enhancer-promoter pairs co-localization within hierarchy of TADs across multiple Hi-C matrices.
- `Pval.raw`: -log10 p-value obtained by testing the canonical correlation coefficient. It quantifies the amount of evidence provided by the data for the presence of synchronized activity between the enhancer and promoter gene.
- `Pval.adapt`: -log10 adjusted p-values by considering the 3D co-localization information encoded in the HC score.
- `GTEx`, `PancanQTL`, `pcHiC`: 1 if supported by the Genotype-Tissue Expression project, the pan-cancer eQTL analysis, or capture Hi-C experiments, respectively.

For any further details, please refer to the *Materials and methods* section of the paper.
