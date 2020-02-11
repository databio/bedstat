library(GenomicDistributions)
library(optparse)
library(tools)

option_list = list(
    make_option(c("--bedfile"), type="character", default=NULL, 
              help="bed file to process", metavar="character"),
	make_option(c("--fileid"), type="character", default=NULL,
              help="fileID to use for output files prefix", metavar="character"),
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

if (is.null(opt$fileid)) {
    print_help(opt_parser)
    stop("fileID input missing.")
}


doitall <- function(query, fname, fileid, genome) {

    ## calculate md5 sum of bedfile
    md5s <- as.vector(md5sum(fname))

    ## continue on with calculations
	TSSdist = calcFeatureDistRefTSS(query, genome)
	g = plotFeatureDist(TSSdist, featureName="TSS")
    # different width/heights for presenting on screen in HTML and for PDFs
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_tssdist.png"), g, device="png", width=8, height=8, units="cm")
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_tssdist.pdf"), g, device="pdf", width=12, height=12, units="cm")
    
	x = calcChromBinsRef(query, genome)
	g = plotChromBins(x)
    # different width/heights for presenting on screen in HTML and for PDFs    
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_chrombins.png"), g, device="png", width=12, height=8, units="cm")
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_chrombins.pdf"), g, device="pdf", width=18, height=12, units="cm")
    
	gcvec = calcGCContentRef(query, genome)
	g = plotGCContent(gcvec)
    # different width/heights for presenting on screen in HTML and for PDFs
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_gccontent.png"), g, device="png", width=8, height=8, units="cm")
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_gccontent.pdf"), g, device="pdf", width=12, height=12, units="cm")
    
	gp = calcPartitionsRef(query, genome)
	# flatten the result returned by the function above
	partiotionNames = as.vector(gp[,"partition"])
	partitionsList = list()
	for(i in seq_along(partiotionNames)){
	    partitionsList[[paste0(partiotionNames[i], "_frequency")]] = 
	        as.vector(gp[,"Freq"])[i]
	    partitionsList[[paste0(partiotionNames[i], "_percentage")]] = 
	        as.vector(gp[,"Freq"])/length(query)[i]
	}
	g = plotPartitions(gp)
    # different width/heights for presenting on screen in HTML and for PDFs    
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_partitions.png"), g, device="png", width=8, height=8, units="cm")
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_partitions.pdf"), g, device="pdf", width=12, height=12, units="cm")
    
	# Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/const.py
	bedmeta = list(
	    id=fileid,
		GC_content=mean(gcvec),
		number_of_regions=length(query),
		mean_absolute_TSS_distance=mean(abs(TSSdist), na.rm=TRUE),
		md5sum=md5s,
		plots=data.frame(name=c("tssdist","chrombins","gccontent","partitions"), caption=c("Region-TSS distance distribution", "Regions distribution over chromosomes", "GC content", "Regions distribution over genomic partitions"))
	)
	write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), paste0(outfolder, "/", fileid, ".json"))
}

# set query to bed file
fileid = opt$fileid
fn = opt$bedfile
outfolder = opt$outputfolder
genome = opt$genome

query = LOLA::readBed(fn)
doitall(query, fn, fileid, genome)