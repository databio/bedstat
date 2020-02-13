library(GenomicDistributions)
library(optparse)
library(tools)

option_list = list(
    make_option(c("--bedfile"), type="character", default=NULL, 
              help="bed file to process", metavar="character"),
	make_option(c("--fileId"), type="character", default=NULL,
              help="fileId to use for output files prefix", metavar="character"),
    make_option(c("--outputfolder"), type="character", default="output",
              help="base output folder for results", metavar="character"),
    make_option(c("--genome"), type="character", default="hg38",
              help="genome to calculate against", metavar="character"))
 
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

plotBoth <- function(plotPth, g){
    print(paste0("Plotting: ", plotPth))
    ggplot2::ggsave(paste0(plotPth, ".png"), g, device="png", width=8, height=8, units="cm")
    ggplot2::ggsave(paste0(plotPth, ".pdf"), g, device="pdf", width=12, height=12, units="cm")
}

doitall <- function(query, fname, fileId, genome) {
    plots = data.frame(stringsAsFactors=F)
    ## calculate md5 sum of bedfile
    md5s <- as.vector(md5sum(fname))

    ## continue on with calculations
	TSSdist = calcFeatureDistRefTSS(query, genome)
	plotId = "tssdist"
	plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
	         plotFeatureDist(TSSdist, featureName="TSS"))
	newPlot = data.frame("name"=plotId, "caption"="Region-TSS distance distribution")
    plots = rbind(plots, newPlot)
    
    
	x = calcChromBinsRef(query, genome)
    plotId = "chrombins"
    plotBoth(paste0(outfolder, "/", fileId, "_", plotId), 
             plotChromBins(x))
    newPlot = data.frame("name"=plotId, "caption"="Regions distribution over chromosomes")
    plots = rbind(plots, newPlot)
    
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
	# Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/const.py
	bedmeta = list(
	    id=fileId,
		gc_content=mean(gcvec),
		regions_no=length(query),
		mean_absolute_TSS_dist=mean(abs(TSSdist), na.rm=TRUE),
		md5sum=md5s,
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
doitall(query, fn, fileId, genome)