require(shiny)
require(DT) # graphical render: A Wrapper of the JavaScript Library 'DataTables'
require(data.table)
require(R.utils)
require(GenomicRanges)

#setwd("/Users/esalviat/Desktop/ShinyApp/DataTable/")


## load data.tabe predictions
Data<-data.table::fread(file="Results/Formatted_Candidate_ETG_20200305_TAD10Kb.tsv.gz")

GR.E<-GenomicRanges::GRanges(
		seqnames=Rle(Data$chr),
		ranges=IRanges(start=as.numeric(Data$start),end=as.numeric(Data$end)),
		strand=Rle(rep("*",nrow(Data))) 
	)

ui<-fluidPage(

	titlePanel("Basic DataTable: 3D ETG prediction"),



  	fluidRow(
  		column(5,
  			sliderInput(
  				inputId="pval",
  				label="FDR threshold:",
  				min=0,max=1,
  				value=0.1,
  				width="100%",
  				round=-2,step=0.01
  			)
  		)
  	),

	fluidRow(
    	column(2,
    		textAreaInput(inputId="symbols",
				label="Gene symbols",
				value="AGRN\nTTLL10\nGATA1",
				width="150px",height="170px")
    	),
    	column(3,
			textAreaInput(inputId="coord",
				label="Enhancers (or genomic) regions",
				value="chr1:1003244-1003582\nchr1:1003777-1004190\nchr1:1004275-1004379\nchr1:1004480-1005034\nchr1:1005102-1006328\nchrX:48345389-48347146\nchrX:48356496-48357210\nchrX:48624027-48624211\nchrX:48652718-48654436",
				width="230px",height="170px"),
    	)	
  	),

  	fluidRow(
  		column(5,verbatimTextOutput(outputId="info.user"))
  	),

  	fluidRow(
  		column(5,verbatimTextOutput(outputId="info.filter"))
  	),

	

	#textAreaInput("input.coord","Enhancer coordinate (chr:start-end)", "chr10:52190910-52191939"),
	#fileInput("file1", "Bed file with enhancer coordinates"),

	DT::dataTableOutput("table",width = "80%")

)


#input<-list()
#input$coord<-"chr1:1003244-1003582\nchr1:1003777-1004190\nchr1:1004275-1004379\nchr1:1004480-1005034\nchr1:1005102-1006328\nchrX:48345389-48347146\nchrX:48356496-48357210\nchrX:48624027-48624211\nchrX:48652718-48654436"
#input$symbols<-"AGRN\nTTLL10\nGATA1"


server<-function(input,output){


	#output$info.filter<-renderText({"Nulla per ora"})

	output$info.user<-renderText({ ## to improve
		paste0("Summary user filters: \n",
			"\U2022 ",length(unlist(strsplit(input$symbols,split="\n")))," selected genes\n",
			"\U2022 ",length(unlist(strsplit(input$coord,split="\n")))," selected enhancers or regions\n",
			"\U2022 -log2(FDR) >= ", round(-log2(input$pval),2) 
		)
	})


	output$table<-DT::renderDataTable({
		
		# initialize
		ind.filter.1<-ind.filter.2<-1:nrow(Data)

		## -- Filter gene symbols -- ##
		symbols<-unlist(strsplit(input$symbols,split="\n"))
		if(length(symbols)>0){ 
			ind.filter.1<-which(Data$symbol %in% symbols) #improve as genes are merged by promters regions
		}

		## -- Filter enhancer regions -- ##
		e.chr<-data.table::data.table(do.call("rbind",strsplit(strsplit(input$coord,split="\n")[[1]],split="[:-]")))
		if(nrow(e.chr)>0){
			colnames(e.chr)<-c("chr","start","end")
			GR.input<-GenomicRanges::GRanges(
				seqnames=Rle(e.chr$chr),
				ranges=IRanges(start=as.numeric(e.chr$start),end=as.numeric(e.chr$end)),
				strand=Rle(rep("*",nrow(e.chr)))
			)

			# find overlaps
			ind.filter.2<-unique(GenomicRanges::findOverlaps(GR.E,GR.input,type="any",select="all")@from)
			output$info.filter<-renderText({ paste0("E:", length(unique(GR.E[ind.filter.2])) ) })

		}

		## -- Filter fdr -- ##
		ind.filter.3<- which(Data$log2.pval.adapt>= -log2(input$pval))



		DT::datatable(
			Data[Reduce(intersect, list(ind.filter.1,ind.filter.2,ind.filter.3)), ],
			options=list(pageLength=15)
		)
		
	})


}

shinyApp(ui,server)

#autoWidth=TRUE,scrollX = TRUE
#options=list(autoWidth=TRUE)

