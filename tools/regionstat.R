library(GenomicDistributions)
library(optparse)

## we need to parameterize LOLACore location
## and file name for bed file

defaultLOLADB = "/ext/qumulo/resources/regions/LOLACore/hg38"

option_list = list(
    make_option(c("--bedfile"), type="character", default=NULL, 
              help="bed file to process", metavar="character"),
    make_option(c("--lolaloc"), type="character", default=defaultLOLADB,
              help="location of LOLA Core database to process", metavar="character"),
	make_option(c("--fileid"), type="character", default=NULL,
              help="fileID to use for output files prefix", metavar="character"),
    make_option(c("--outputfolder"), type="character", default="output",
              help="base output folder for results", metavar="character"))
 
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

doitall <- function(query, fileid) {

	TSSdist = calcFeatureDistRefTSS(query, "hg38")
	g = plotFeatureDist(TSSdist, featureName="TSS")
	ggplot2::ggsave(paste0(outfolder, fileid, "_tssdist.png"), g)

	x = calcChromBinsRef(query, "hg38")
	g = plotChromBins(x)
	ggplot2::ggsave(paste0(outfolder, fileid, "_chrombins.png"), g)

	gcvec = calcGCContentRef(query, "hg38")
	g = plotGCContent(gcvec)
	ggplot2::ggsave(paste0(outfolder, fileid, "_gccontent.png"), g)

	gp = calcPartitionsRef(query, "hg38")
	gp$Perc = gp$Freq/length(query)
	g = plotPartitions(gp)
	ggplot2::ggsave(paste0(outfolder, fileid, "_partitions.png"), g)

	l = list()
	bedmeta = list(id=fileid,
		gc_content=mean(gcvec),
		num_regions=length(query),
		mean_abs_tss_dist=mean(abs(TSSdist), na.rm=TRUE),
		genomic_partitions=gp)
	l[[fileid]]=bedmeta
	l

	write(jsonlite::toJSON(l, pretty=TRUE), paste0(outfolder, fileid,".json"))

}

# set query to bed file
rdb = LOLA::loadRegionDB(opt$lolaloc)
fileid = opt$fileid
fn = opt$bedfile

outfolder = opt$outputfolder  # Set to '' for cwd

#' query = rtracklayer::import(fn)
query <- LOLA::readBed(fn)
doitall(query,fileid)
