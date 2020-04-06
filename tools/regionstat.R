library(GenomicDistributions)
library(optparse)
library(tools)

option_list = list(
    make_option(c("--bedfile"), type="character", default=NULL, 
              help="path to a BED file to process", metavar="character"),
	make_option(c("--fileId"), type="character", default=NULL,
              help="BED file ID to use for output files prefix", metavar="character"),
	make_option(c("--openSignalMatrix"), type="character", default=NULL,
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

doitall <- function(query, fname, fileId, genome, cellmatrix=NULL) {
    plots = data.frame(stringsAsFactors=F)

    ## continue on with calculations
	TSSdist = calcFeatureDistRefTSS(query, genome)
	plotId = "tssdist"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotFeatureDist(TSSdist, featureName="TSS"))
	newPlot = data.frame("name"=plotId, "caption"="Region-TSS distance distribution")
    plots = rbind(plots, newPlot)
    
    
	#x = calcChromBinsRef(query, genome)
    #plotId = "chrombins"
    #plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
             #plotChromBins(x))
    #newPlot = data.frame("name"=plotId, "caption"="Regions distribution over chromosomes")
    #plots = rbind(plots, newPlot)
    
	gcvec = calcGCContentRef(query, genome)
	plotId = "gccontent"
    plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
             plotGCContent(gcvec))
    newPlot = data.frame("name"=plotId, "caption"="GC content")
    plots = rbind(plots, newPlot)
    
	gp = calcPartitionsRef(query, genome)
	plotId = "partitions"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotPartitions(gp))
	newPlot = data.frame("name"=plotId, "caption"="Regions distribution over genomic partitions")
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
	
	# Add plots from hackaton
	# Add QThist plot
	widths = calcWidth(query)
	plotId = "widths_histogram"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
		plotQTHist(widths))
	newPlot = data.frame("name"=plotId, "caption"="Quantile-Trimmed Histogram of Widths")
	plots = rbind(plots, newPlot)

	# Add tissue specificity plot
	if (!is.null(cellmatrix)) {
		op = calcOpenSignal(query, cellmatrix)
		plotId = "open_chromatin"
		plotBoth(paste0(outfolder, "/", fileId, "_", plotId),
			plotOpenSignal(op))
		newPlot = data.frame("name"=plotId, "caption"="Cell specific enrichment for open chromatin")
		plots = rbind(plots, newPlot)
	}

	# Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/const.py
	bedmeta = list(
	    id=fileId,
		gc_content=mean(gcvec),
		regions_no=length(query),
		mean_absolute_TSS_dist=mean(abs(TSSdist), na.rm=TRUE),
		mean_region_width=mean(widths),
		md5sum=opt$digest,
		plots=plots,
		bedfile_path=fname
	)
	write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), paste0(outfolder, "/", fileId, ".json"))
}

# set query to bed file
fileId = opt$fileId
fn = opt$bedfile
outfolder = opt$outputfolder
genome = opt$genome

query = LOLA::readBed(fn)
if (!is.null(opt$openSignalMatrix)){
	osm = opt$openSignalMatrix
	cellMatrix = data.table::fread(osm)
	doitall(query, fn, fileId, genome, cellMatrix)
} else {
	doitall(query, fn, fileId, genome)
}

