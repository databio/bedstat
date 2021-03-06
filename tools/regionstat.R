library(GenomicDistributions)
library(GenomicDistributionsData)
library(optparse)
library(tools)
data(TSS_hg38)

option_list = list(
    make_option(c("--bedfile"), type="character", default=NULL, 
              help="path to a BED file to process", metavar="character"),
	make_option(c("--fileId"), type="character", default=NULL,
              help="BED file ID to use for output files prefix", metavar="character"),
	make_option(c("--openSignalMatrix"), type="character",
			  help="path to the open signal matrix required for the tissue specificity plot", metavar="character"),
    make_option(c("--digest"), type="character", default=NULL,
              help="digest of the BED file", metavar="character"),
    make_option(c("--outputfolder"), type="character", default="output",
              help="base output folder for results", metavar="character"),
    make_option(c("--genome"), type="character", default="hg38",
              help="genome reference to calculate against", metavar="character"))
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);    

if (is.null(opt$bedfile)) {
    print_help(opt_parser)
    stop("Bed file input missing.")
}

if (is.null(opt$fileId)) {
    print_help(opt_parser)
    stop("fileId input missing.")
}

if (is.null(opt$digest)) {
    print_help(opt_parser)
    stop("digest input missing.")
}


plotBoth <- function(plotPth, g){
    print(paste0("Plotting: ", plotPth))
    ggplot2::ggsave(paste0(plotPth, ".png"), g, device="png", width=8, height=8, units="in")
    ggplot2::ggsave(paste0(plotPth, ".pdf"), g, device="pdf", width=8, height=8, units="in")
}

doItAall <- function(query, fname, fileId, genome, cellMatrix) {
    plots = data.frame(stringsAsFactors=F)
    bsGenomeAvail = ifelse((requireNamespace(BSg, quietly=TRUE) | requireNamespace(BSgm, quietly=TRUE)), TRUE, FALSE)
    ## continue on with calculations
	TSSdist = calcFeatureDistRefTSS(query, genome)
	plotId = "tssdist"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotFeatureDist(TSSdist, featureName="TSS"))
	newPlot = data.frame("name"=plotId, "caption"="Region-TSS distance distribution")
    plots = rbind(plots, newPlot)
    
    
    # Chromosomes region distribution plot
	x = calcChromBinsRef(query, genome)
    plotId = "chrombins"
    plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
             plotChromBins(x))
    newPlot = data.frame("name"=plotId, "caption"="Regions distribution over chromosomes")
    plots = rbind(plots, newPlot)
    
	# OPTIONAL: Plot GC content only if proper BSgenome package is installed. 
	if (bsGenomeAvail) {
		gcvec = calcGCContentRef(query, genome)
		plotId = "gccontent"
		plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
					plotGCContent(gcvec))
		newPlot = data.frame("name"=plotId, "caption"="GC content")
		plots = rbind(plots, newPlot)
	}
    # Partition Plots, default to percentages
	gp = calcPartitionsRef(query, genome)
	plotId = "partitions"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotPartitions(gp))
	newPlot = data.frame("name"=plotId, "caption"="Regions distribution over genomic partitions")
	plots = rbind(plots, newPlot)

	ep = calcExpectedPartitionsRef(query, genome)
	plotId = "expected_partitions"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotExpectedPartitions(ep))
	newPlot = data.frame("name"=plotId, "caption"="Expected distribution over genomic partitions")
	plots = rbind(plots, newPlot)

	cp = calcCumulativePartitionsRef(query, genome)
	plotId = "cumulative_partitions"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
			 plotCumulativePartitions(cp))
	newPlot = data.frame("name"=plotId, "caption"="Cumulative distribution over genomic partitions")
	plots = rbind(plots, newPlot)


	# flatten the result returned by the function above
	partiotionNames = as.vector(gp[,"partition"])
	partitionsList = list()
	for(i in seq_along(partiotionNames)){
	    partitionsList[[paste0(partiotionNames[i], "_frequency")]] = 
	        as.vector(gp[,"Freq"])[i]
	    partitionsList[[paste0(partiotionNames[i], "_percentage")]] = 
	        as.vector(gp[,"Freq"])[i]/length(query)	        
	}
	
	# QThist plot
	widths = calcWidth(query)
	plotId = "widths_histogram"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
		plotQTHist(widths))
	newPlot = data.frame("name"=plotId, "caption"="Quantile-Trimmed Histogram of Widths")
	plots = rbind(plots, newPlot)

	# Neighbor regions distance plots
	dist = calcNeighborDist(query)
	plotId = "neighbor_distances"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotNeighborDist(dist))
	newPlot = data.frame("name"=plotId, "caption"="Distance between neighbor regions")
	plots = rbind(plots, newPlot)

	# OPTIONAL: Add tissue specificity plot if open signal matrix is provided
	if (cellMatrix == "None") {
		message("open signal matrix not provided. Skipping tissue specificity plot ... ")
	} else {
		matrix = data.table::fread(cellMatrix)
		op = calcOpenSignal(query, matrix)
		plotId = "open_chromatin"
		plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
		         plotOpenSignal(op))
		newPlot = data.frame("name"=plotId, "caption"="Cell specific enrichment for open chromatin")
		plots = rbind(plots, newPlot)
	}

	# Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/const.py
	bedmeta = list(
	    id=fileId,
		gc_content=ifelse(bsGenomeAvail, mean(gcvec), NA),
		regions_no=length(query),
		mean_absolute_TSS_dist=mean(abs(TSSdist), na.rm=TRUE),
		mean_region_width=mean(widths),
		md5sum=opt$digest,
		plots=plots,
		bedfile_path=fname
	)
	write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), paste0(outfolder, "/", fileId, ".json"))
}

# define values and output folder for doitall()
fileId = opt$fileId
fn = opt$bedfile
outfolder = opt$outputfolder
genome = opt$genome
cellMatrix = opt$openSignalMatrix
orgName = "Mmusculus"

# build BSgenome package ID to check whether it's installed
if (startsWith(genome, "hg") | startsWith(genome, "grch")) orgName = "Hsapiens"

BSg = paste0("BSgenome.", orgName , ".UCSC.", genome)
BSgm = paste0(BSg, ".masked")

# read bed file and run doitall()
query = LOLA::readBed(fn)
doItAall(query, fn, fileId, genome, cellMatrix)


