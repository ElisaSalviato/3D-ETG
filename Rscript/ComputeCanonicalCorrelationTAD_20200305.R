require(GenomicRanges)
require(data.table)
require(ccaPP)
require(Matrix)
require(limma)
require(pbapply)
require(candisc)

# Working directory: 3D-ETG/

#### INPUT parameters: ###
dir.Matrix<-"Data/Roadmap/Matrix/"
dir.result<-"Results/3D-ETG/"


dir.TAD<-"Data/TAD/"
method.tad<-"LSD"
bin.size<-10000
# key word to identify TAD files
key.TAD<-c("TAD",method.tad)



 
## Use the key words to find the name of file
key.enhancer<-c("enhancer_consolidated","max","20200305")
enhancer.activity<-c("DNase","H3K27ac")

key.promoter<-c("promoter_consolidated","max","20200305")
promoter.activity<-c("DNase","H3K27ac","H3K4me3")

CHR<-rev(c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX"))



#########################################################
cat("START \n\n")


### ------------ load utility functions ------------ ###
source("Rscript/UtilityFunction_20200305.R")



E.files<-find.Files(dir.Matrix,key.enhancer)
names(E.files)<-sapply(E.files,function(x) strsplit(x,split="_")[[1]][1],USE.NAMES=FALSE)
P.files<-find.Files(dir.Matrix,key.promoter)
names(P.files)<-sapply(P.files,function(x) strsplit(x,split="_")[[1]][1],USE.NAMES=FALSE)

## read and split by chromosome (definition names are extra columns)
E.data<-lapply(E.files,read.ChipSeq.Matrix,path.directory=dir.Matrix,definition.name=c("chr","start","end"))
E.data<-E.data[enhancer.activity]

P.data<-lapply(P.files,read.ChipSeq.Matrix,path.directory=dir.Matrix,definition.name=c("chr","start","end","symbol","strand"))
P.data<-P.data[promoter.activity]

celltype<-names(which(table(unlist(lapply(E.data,function(x) colnames(mcols(x[[1]])) )))==length(E.data)))

## -------- ##



cat("\n\n\n\n\n")
### ------------- Compute CANONICAL CORRELATION and HC score ------------------ ###
###  For each chromosome ###

for(chr in CHR){

	cat("--> ",chr,": \n",sep="")

	## Check if file already exist ##
	# if exist: do you want to update the validation columns?
	name.file.result<-paste0(paste("ETG","Cca",method.tad,chr,bin.size,sep="_"),".Rdata")
	if(name.file.result %in% dir(dir.result)) next


	####### Compute Weights ############
	nE<-length(E.data[[1]][[chr]])
	nP<-length(P.data[[1]][[chr]])
	## Mid point enhancers promoters
	GR.E.mid<-GenomicRanges::GRanges(
			seqnames=Rle(E.data[[1]][[chr]]@seqnames),
			ranges=IRanges(start=E.data[[1]][[chr]]@ranges@start+(E.data[[1]][[chr]]@ranges@width/2),width=1),
			strand=Rle(rep("*",length(E.data[[1]][[chr]]))) 
			)

	GR.P.mid<-GenomicRanges::GRanges(
			seqnames=Rle(P.data[[1]][[chr]]@seqnames),
			ranges=IRanges(start=P.data[[1]][[chr]]@ranges@start+(P.data[[1]][[chr]]@ranges@width/2),width=1),
			strand=Rle(rep("*",length(P.data[[1]][[chr]]))) 
			)

	## Read TAD file 
	TAD.files<-find.Files(dir.TAD,c(key.TAD,paste0("_",chr,"_"),paste0("_",bin.size,"_") ))
	names(TAD.files)<-sapply(strsplit(TAD.files,split="[_.]"),function(x) paste(x[-c(1:4,length(x))],collapse="_") )
	

	## Compute the HC score matrices for each HiC matrix dataset and each level of hierarchy
	cat("HC score.. ",sep="")
	D.weight<-lapply(TAD.files,dataset.chr.tad.weight,path.directory=dir.TAD,E.list=GR.E.mid,P.list=GR.P.mid,bin.size=bin.size,use.resolution=TRUE)
	# Merged all the matrix toghether
	W.chr<-Reduce("+",D.weight)


	### Discarded poorly-supported EP pairs (HC score ≤ th.weight) ########
	# by default equal to the number of HiC matrices
	th.weight<-length(TAD.files)
	W.chr[W.chr<=th.weight]<-0

	########### Compute Canonical correlation ############
	index.candidate<-data.table::data.table(summary(W.chr))
	colnames(index.candidate)<-c("E.ind","P.ind","Weight")


	### Quantile nromalization by chromosome ###
	E.chr.data<-lapply(E.data,function(x,chr,cell) { 
		mat<-log2(as.matrix( mcols(x[[chr]])[cell]  )+1)
		matN<-limma::normalizeQuantiles(mat)
		return(matN)
	},chr=chr,cell=celltype)
	
	P.chr.data<-lapply(P.data,function(x,chr,cell) { 
		mat<-log2(as.matrix( mcols(x[[chr]])[cell]  )+1)
		matN<-limma::normalizeQuantiles(mat)
		return(matN)
	},chr=chr,cell=celltype)

	## check matched E.definition/P.definition
	E.code<-data.table::data.table(sapply(E.data,function(x) paste0(x[[chr]]) ))
	pairs<-t(combn(1:length(E.data),2))
	E.res<-apply(pairs,1,function(i) all.equal(E.code[[i[1]]],E.code[[i[2]]] ))

	P.code<-data.table::data.table(sapply(P.data,function(x) paste0(x[[chr]]) ))
	pairs<-t(combn(1:length(P.data),2))
	P.res<-apply(pairs,1,function(i) all.equal(P.code[[i[1]]],P.code[[i[2]]] ))

	if(!all(c(E.res,P.res))) stop("Enhancers/Promoters in different chromatine marks don't match")
	### ---------------------------------------- ###

	cat("correlation [",nrow(index.candidate),"].. ",sep="")

	start.time <- Sys.time()
	Cca.chr<-canonical.correlation.chr(Index= index.candidate,E.chr.data=E.chr.data,P.chr.data=P.chr.data,common.samples=celltype)
	end.time <- Sys.time()

	E.chr<-granges(E.data[[1]][[chr]])[index.candidate$E.ind]
	P.chr<-granges(P.data[[1]][[chr]])[index.candidate$P.ind]

	distance<-GR.E.mid@ranges@start[index.candidate$E.ind]-GR.P.mid@ranges@start[index.candidate$P.ind]
	EP.code<-paste(E.code[[1]][index.candidate$E.ind],P.code[[1]][index.candidate$P.ind],sep=";")

	Cca.chr.res<-data.table::data.table(EP.code=EP.code,Cca.chr,Distance=distance)

	## Add attributes ##
	attr(Cca.chr.res,"time")<-list(start=start.time,end=end.time,elapsed=end.time-start.time)
	attr(Cca.chr.res,"TAD")<-TAD.files
	attr(Cca.chr.res,"Activity")<-list(Enhancer=enhancer.activity,Promoter=promoter.activity)

	save(Cca.chr.res,E.chr,P.chr,file=paste0(dir.result,name.file.result))

	cat("\n")
}


cat("DONE! \n\n")














