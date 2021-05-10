require(shiny)
require(DT) # graphical render: A Wrapper of the JavaScript Library 'DataTables'
require(data.table)
require(R.utils)
require(GenomicRanges)


## -- load complete list of candifate ETG pairs -- ##
Data<-data.table::fread(file="Results/Formatted_Candidate_ETG_20200305_TAD10Kb.tsv.gz")

## Enhancer regions
GR.E<-GenomicRanges::GRanges(
		seqnames=Rle(Data$chr),
		ranges=IRanges(start=as.numeric(Data$start),end=as.numeric(Data$end)),
		strand=Rle(rep("*",nrow(Data))) 
	)

tab<-strsplit(Data$symbol,split=";")
n<-sapply(tab,length)

## Gene symbols (Note: genes with overlapping promoter regions are merged toghether)
G.splitted<-data.table::data.table(
	symbol=unlist(tab),
	index=rep(1:nrow(Data),times=n)
)

rm(tab); rm(n)





server<-function(input,output){

	message<-reactive({
		paste0("Summary user filters: \n",
			"\U2022 ",length(unlist(strsplit(input$symbols,split="\n")))," selected genes\n",
			"\U2022 ",length(unlist(strsplit(input$coord,split="\n")))," selected enhancers or regions\n",
			"\U2022 -log2(FDR) >= ", round(-log2(input$pval),2) 
		)
	})

	output$info.user<-renderText( message() )

	data<-reactive({

		# initialize
		ind.filter.1<-ind.filter.2<-1:nrow(Data)

		## -- 1. Filter gene symbols -- ##
		symbols<-unlist(strsplit(input$symbols,split="\n"))
		if(length(symbols)>0){ 
			ind.filter.1<-G.splitted$index[G.splitted$symbol %in% symbols]
		}

		## -- 2. Filter enhancer regions -- ##
		e.chr<-data.table::data.table(do.call("rbind",strsplit(strsplit(input$coord,split="\n")[[1]],split="[:-]")))
		if(nrow(e.chr)>0){
			colnames(e.chr)<-c("chr","start","end")
			GR.input<-GenomicRanges::GRanges(
				seqnames=Rle(e.chr$chr),
				ranges=IRanges(start=as.numeric(e.chr$start),end=as.numeric(e.chr$end)),
				strand=Rle(rep("*",nrow(e.chr)))
			)

			# find overlaps
			ind.over<-GenomicRanges::findOverlaps(GR.E,GR.input,type="any",select="all")
			ind.filter.2<-unique(ind.over@from)
			
			#output$info.filter<-renderText({ paste0("E:", length(unique(GR.E[ind.filter.2])) ) })

		}

		## -- 3. Filter fdr -- ##
		ind.filter.3<- which(Data$log2.pval.adapt>= -log2(input$pval))


		ind.select<-Reduce(intersect, list(ind.filter.1,ind.filter.2,ind.filter.3))

		## -- Summary filtered results -- ##
		output$info.filter<-renderText({ 
			paste0("Summary of the filters: \n",
			"\U2022 ",sum(symbols %in% G.splitted$symbol)," out of ",length(symbols)," provided genes identified\n",
			"\U2022 ",length(unique(ind.over@from))," enhancers overlapping ",nrow(e.chr)," provided regions\n",
			"\U2022 ",length(ind.filter.3), " pairs with -log2(FDR) >= ",round(-log2(input$pval),2),"\n\n",
			"Total (intersection): ",length(ind.select)," ETG pairs"
			)
		})


		Data[ind.select, ]
	})



	output$size<-renderText({ object.size() })

	output$table<-DT::renderDataTable({	
		DT::datatable(data(),options=list(pageLength=15))
	})


	output$downloadData <- downloadHandler(
	   filename = function() {
	     paste('Salviato_etal_3DETG_', Sys.Date(), '.tsv.gzip', sep='')
	   },
	   content = function(file) {
	     data.table::fwrite(data(), file, sep="\t",col.names=TRUE,compress="gzip",quote=FALSE)
	   }
	)


}


ui<-fluidPage(

	br(),
	h2("Leveraging three-dimensional chromatin architecture for effective reconstruction of enhancer-target gene regulatory interactions"),
	h5("Elisa Salviato, Vera Djordjilović, Judith M. Hariprakash, Ilario Tagliaferri, Koustav Pal, Francesco Ferrari"),
	h5(a("doi: https://doi.org/10.1101/2021.03.01.432687",href="https://doi.org/10.1101/2021.03.01.432687")),
	br(),

	sidebarLayout(position="right",
		sidebarPanel(width=4,
			helpText("Clean up all the preset filters to obtained the complete list of ETG interactions."),	
			fluidRow(
				column(11,
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
	   			column(4,
	    			textAreaInput(inputId="symbols",
					label="Gene symbols:",
					value="AGRN\nTTLL10\nGATA1\nABCB7",
					width="130px",height="140px")
	    		),
	    		column(8,
					textAreaInput(inputId="coord",
					label="Enhancer or genomic regions",
					value="chr1:1003244-1003582\nchr1:1003777-1004190\nchr1:1004275-1004379\nchr1:1004480-1005034\nchr1:1005102-1006328\nchrX:100000-150000000",
					width="230px",height="140px"),
	    		)
    		),
    		fluidRow(
  				column(12,verbatimTextOutput(outputId="info.filter"))
  			),
  			fluidRow(
  				column(width=12,downloadButton("downloadData","Download the ETG filtered list (.gzip)"))
  			)	
  		),
		mainPanel(width=7,
			br(),
			p(style="text-align: justify;",
				"The following table provide the complete list of candidate enhancer-target 
				pairs (ETG) inferred in ",a("Salviato et al. ",href="https://doi.org/10.1101/2021.03.01.432687"), 
				"The schematic workflow of our methodological framework is available at ",a("Figure 1",href="https://doi.org/10.1101/2021.03.01.432687")," of the paper.","The source code of the method is available in the ",a("3D-ETG",href="https://github.com/ElisaSalviato/3D-ETG")," github repository."),
			p(style="text-align: justify;",
				"Each line corresponds to an ETG pair, where:",
			HTML(
				"<ul style='padding-bottom: 5px;'>
					<li style='padding-bottom: 6px;'> <b>chr</b>, <b>start</b>, <b>end</b>: enhancer region (hg19). Promoter proximal elements are not considered (3.5 kb upstream and 1.5 kb downstream of coding and non-coding TSS). </li>
					<li style='padding-bottom: 6px;'> <b>symbol</b>: target-gene symbol of protein coding genes, based on RefSeq annotations in UCSC. Target genes with overlapping promoter regions are merged as a single pairs (separated by semicolon). </li>
					<li style='padding-bottom: 6px;'> <b>distance</b>: distance between mid-points of gene promoter and enhancer regions. Promoters are defined as 1.5kb upstream and 0.5kb downstream regions of TSS of coding genes.</li>
					<li style='padding-bottom: 6px;'> <b>cca</b>: first canonical correlation coefficient. It is based on the enrichment of DNase-seq and H3K27ac for enhancers and DNase-seq, H3K27ac and H3K4me3 for promoters. </li>
					<li style='padding-bottom: 6px;'> <b>HC</b>: Hierarchical Contact (HC) score. HC score is proportional to the likelihood of enhancer-promoter pairs co-localization within hierarchy of TADs across multiple Hi-C matrices.</li>
					<li style='padding-bottom: 6px;'> <b>Pval.raw</b>: -log10 p-value obtained by testing the canonical correlation coefficient. It quantifies the amount of evidence provided by the data for the presence of synchronized activity between the enhancer and promoter gene.</li>
					<li style='padding-bottom: 6px;'> <b>Pval.adapt</b>: -log10 adjusted p-values by considering the 3D co-localization information encoded in the HC score. </li>
					<li> <b>GTEx</b>, <b>PancanQTL</b>, <b>pcHiC</b>: 1 if supported by the Genotype-Tissue Expression project, the pan-cancer eQTL analysis, or capture Hi-C experiments, respectively. </li>
				</ul>")),
			br(),
			p("For more details, please refer to the ",em("Materials and methods")," section of the paper.")
		)
	),
	br(),br(),

	

	DT::dataTableOutput("table",width="90%")

)





shinyApp(ui,server)


