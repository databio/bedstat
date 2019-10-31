library(GenomicDistributions)
library(optparse)
library(tools)

## we need to parameterize LOLACore location
## and file name for bed file

defaultLOLADB = "/ext/qumulo/resources/regions/LOLACore/hg38"

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

##rdb = LOLA::loadRegionDB("/home/maketo/dev/sheffield/ext/qumulo/LOLAweb/databases/LOLACore/hg38")

##rdb$regionGRL
##rdb$regionAnno[1,"filename"]

##length(rdb$regionGRL)
##lapply(seq(1, 3), 
##	function(x) {
##		doitall(rdb$regionGRL[[x]], rdb$regionAnno[x,]$filename)
##	}
##	)

##fileid <- "E125-DNase.macs2"
##fn <- "/home/maketo/dev/sheffield/ext/qumulo/LOLAweb/databases/LOLACore/hg38/cistrome_cistrome/regions/Human_HepG2-liver-cells_FoxA2_No-treatment_Wadelius.bed"

doitall <- function(query, fname, fileid, genome) {

    ## calculate md5 sum of bedfile
    md5s <- as.vector(md5sum(fname))

    ## continue on with calculations
	TSSdist = calcFeatureDistRefTSS(query, genome)
	g = plotFeatureDist(TSSdist, featureName="TSS")
    # different width/heights for presenting on screen in HTML and for PDFs
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_tssdist.png"), g, device="png", width=12, height=12, units="cm")
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_tssdist.pdf"), g, device="pdf", width=15, height=15, units="cm")
    
	x = calcChromBinsRef(query, genome)
	g = plotChromBins(x)
    # different width/heights for presenting on screen in HTML and for PDFs    
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_chrombins.png"), g, device="png", width=12, height=12, units="cm")
    ggplot2::ggsave(paste0(outfolder, "/", fileid, "_chrombins.pdf"), g, device="pdf", width=15, height=15, units="cm")
    
	gcvec = calcGCContentRef(query, genome)
	g = plotGCContent(gcvec)
    # different width/heights for presenting on screen in HTML and for PDFs
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_gccontent.png"), g, device="png", width=12, height=12, units="cm")
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_gccontent.pdf"), g, device="pdf", width=15, height=15, units="cm")
    
	gp = calcPartitionsRef(query, genome)
	gp$Perc = gp$Freq/length(query)
	g = plotPartitions(gp)
    # different width/heights for presenting on screen in HTML and for PDFs    
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_partitions.png"), g, device="png", width=12, height=12, units="cm")
	ggplot2::ggsave(paste0(outfolder, "/", fileid, "_partitions.pdf"), g, device="pdf", width=15, height=15, units="cm")
    
	##l = list()
	bedmeta = list(id=fileid,
		gc_content=mean(gcvec),
		num_regions=length(query),
		mean_abs_tss_dist=mean(abs(TSSdist), na.rm=TRUE),
		genomic_partitions=gp,
        md5sum=md5s)
	##l[[fileid]]=bedmeta
	##l

	write(jsonlite::toJSON(bedmeta, pretty=TRUE), paste0(outfolder,"/",fileid,".json"))

}

# set query to bed file
fileid = opt$fileid
fn = opt$bedfile
outfolder = opt$outputfolder
genome = opt$genome

#' query = rtracklayer::import(fn)
query = LOLA::readBed(fn)
doitall(query,fn,fileid,genome)
