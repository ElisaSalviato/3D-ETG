## Script to correct pvalues using score information: last version ##
# read cca files
# calculate correction
# add pvalues columns

require(adaptMT) #devtools::install_github("lihualei71/adaptMT")
require(splines)
require(data.table)

### ------------ load utility functions ------------ ###
source("/storage/home/esalviat/MatrixETG/CanonicalCorrelation_Version_20190627/UtilityFunction_CanonicalCorrelation.R")
# ---------------------------------------------------- #

## INPUTS:
dir.result<-"Results/3D-ETG/"
# Key words to identify chromosome-wise ETG files
key.cca.files<-c("ETG","Cca")

## ----- AdaPT model parameters for Generalized Linear Models ----- ##
# Five natural cubic splines
alphas<-seq(0.01,0.1,by=0.005)
formulas <- paste0("ns(x, df = ", c(5:10), ")")


## ------------------- #########

cat("START \n\n")
name.cca<-find.Files(dir.result,key.cca.files)
CHR<-sapply(strsplit(name.cca,split="_"),function(x) x[4])

for(chr in CHR){
	
	cat(chr,": \n\tload..\n")
	i<-grep(paste0("_",chr,"_"),name.cca)

	# load: "Cca.chr.res" "E.chr" "P.chr" 
	load(paste0(dir.result,name.cca[i]))

	if( sum(names(Cca.chr.res) %in% c("Pval.adapt") )> 0){
		cat(chr," pvalues already adjusted.. skipe!\n\n")
		next
	}

	## Delete NA Cca:  pairs with no enough variability to calculate canonical correlation (flat signal) ##
	ind.na<-which(!is.na(Cca.chr.res$Cca))
	Cca.chr.res<-Cca.chr.res[ind.na, ]
	E.chr<-E.chr[ind.na,]
	P.chr<-P.chr[ind.na,]

	## ----- adapt correction ----- ##
	# Predictor (x): HC score (i.e., Weight column)
	pvals.chr<-Cca.chr.res$Pval
	x.chr<-data.frame(x=Cca.chr.res$Weight)

	# -- Time consuming step -- #
	start.time <- Sys.time()
	res.glm.chr <- adaptMT::adapt_glm(x = x.chr, pvals = pvals.chr, pi_formulas = formulas, mu_formulas = formulas,alphas=alphas)
	end.time <- Sys.time()

	time.chr<-end.time-start.time
	## --------------------------- ##
	cat("\t time: ",round(time.chr,2), attr(time.chr,"units"),"\n")
	
	Cca.chr.res<-data.table::data.table(
		Cca.chr.res,
		Pval.adapt=pmin(res.glm.chr$qval,1),
		Pval.bonf=p.adjust(pvals.chr,method="bonferroni"),
		Pval.bh=p.adjust(pvals.chr,method="BH")
	)
	cat("\t saved!\n")
	
	save(Cca.chr.res,E.chr,P.chr, file=paste0(dir.result,name.cca[i]))
}


