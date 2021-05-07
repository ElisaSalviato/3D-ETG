######### UTILITY FUNCTION ############

find.Files<-function(path.directory,key.words){
	directory.files<-dir(path.directory)
	index<-lapply(key.words,function(x,names) grep(x,names),names=directory.files)
	return(directory.files[Reduce("intersect",index)])
}


read.ChipSeq.Matrix<-function(name.file,path.directory,definition.name){
	Di<-data.table::fread(file=paste0(path.directory,name.file))
	
	name.data<-setdiff(colnames(Di),definition.name)
	name.data<-c(definition.name[-c(1:3)],name.data)
	Di.data<-Di[,..name.data]

	## correct name ##
	ind<-which(colnames(Di.data) %in% c("seqnames", "ranges",
    "strand", "seqlevels", "seqlengths", "isCircular", "start", "end",
    "width", "element"))
    if(length(ind)>0) colnames(Di.data)[ind]<-paste0("G.",colnames(Di.data)[ind])

	GR.Di<-GRanges(
			seqnames=Rle(Di[[ definition.name[1] ]]),
			ranges=IRanges( start=as.numeric(Di[[ definition.name[2] ]]),end=as.numeric(Di[[ definition.name[3] ]]) ),
			strand=Rle(rep("*",nrow(Di)) )
	)
	mcols(GR.Di)<-Di.data

	GR.Di<-split(GR.Di,seqnames(GR.Di))
	return(GR.Di)
}


dataset.chr.tad.weight<-function(name.file,path.directory,E.list,P.list,bin.size=5000,di.window=NULL,use.resolution=TRUE){

	# @E.list: GenomicRanges object with midpoint of enhancer regions
	# @P.list: GenomicRanges object with midpoint of promoter regions
	# @use.resolution: if TRUE we use a weight for the di.window 

	tad<-get(load(paste0(path.directory,name.file)))[[as.character(bin.size)]]
	if(!is.null(di.window)){ tad<-tad[paste(di.window)] }

	di.window<-names(tad)
	nE<-length(E.list)
	nP<-length(P.list)

	if(use.resolution) balance<- sqrt(max(as.numeric(di.window))/as.numeric(di.window))
	else balance<- rep(1,length(di.window))
	names(balance)<-di.window

	Result<-vector("list",length(di.window))
	names(Result)<-di.window
	for(w in 1:length(di.window)){

		#w<-di.window[i]
		nTw<-length(tad[[w]])
		E.mat<-Matrix::Matrix(0,nrow=nE,ncol=nTw)
		P.mat<-Matrix::Matrix(0,nrow=nP,ncol=nTw)

		ind.w.e<-GenomicRanges::findOverlaps(E.list,tad[[w]],type="within")
		ind.w.p<-GenomicRanges::findOverlaps(P.list,tad[[w]],type="within")

		E.mat[matrix(c(ind.w.e@from,ind.w.e@to),ncol=2)]<-1
		P.mat[matrix(c(ind.w.p@from,ind.w.p@to),ncol=2)]<-1

		Result[[w]]<- E.mat %*% t(P.mat) * balance[w]
	}
	
	Result<-Reduce("+",Result)

}


canonical.correlation.chr<-function(Index,E.chr.data,P.chr.data,common.samples){


	Res.candisc<-pbapply::pbsapply(1:nrow(Index),function(i,index,Ed,Pd,cs){

		ind<-index[i,]

		x<-sapply(Ed,function(z,j,cs) z[j,cs],j=ind[[1]],cs=cs)
		y<-sapply(Pd,function(z,j,cs) z[j,cs],j=ind[[2]],cs=cs)
		

		## Some times there are too much 0's ##
		mod<-tryCatch(
			candisc::cancor(x=x,y=y),
			error=function(err){ return(NA) })

		if(class(mod)=="cancor"){
			pval<-candisc::Wilks(mod)[[6]][1]

			interaction<-c(pval,mod$cancor[1])
		} else { interaction<-c(NA,NA) }

		names(interaction)<-c("pvalue","cancor")
		return(interaction)

	},index=Index,Ed=E.chr.data,Pd=P.chr.data,cs=common.samples)


	Result.chr<-data.table(Pval=Res.candisc[1,],Cca=Res.candisc[2,],Weight=Index$Weight)
	return(Result.chr)
}
