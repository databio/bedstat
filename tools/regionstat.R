library(GenomicDistributions)
library(GenomicDistributionsData)
library(GenomeInfoDb)
library(ensembldb)
library(optparse)
library(tools)
library(R.utils)

trim <- IRanges::trim

option_list = list(
  make_option(c("--bedfilePath"), type="character", default=NULL, 
              help="full path to a BED file to process", metavar="character"),
  make_option(c("--fileId"), type="character", default=NULL,
              help="BED file ID to use for output files prefix", metavar="character"),
  make_option(c("--openSignalMatrix"), type="character",
              help="path to the open signal matrix required for the tissue specificity plot", metavar="character"),
  make_option(c("--digest"), type="character", default=NULL,
              help="digest of the BED file", metavar="character"),
  make_option(c("--outputFolder"), type="character", default="output",
              help="base output folder for results", metavar="character"),
  make_option(c("--genome"), type="character", default="hg38",
              help="genome reference to calculate against", metavar="character"),
  make_option(c("--ensDbGtf"), type="character",
              help="path to the Ensembl annotation gtf file", metavar="character")
)


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);    

if (is.null(opt$bedfilePath)) {
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

buildEnsDb <- function(genome){
  ## generate the SQLite database file
  DB <- ensDbFromGtf(gtf = ensDbGtf)
  ## load the DB file directly
  EDB <- EnsDb(DB)
  
  return(EDB)
}

myGeneModels <- function(genome) {
  geneModels = tryCatch({ 
    message("building geneModels using Ensembl")
    EnsDb = buildEnsDb(genome)
    codingFilter = AnnotationFilter::AnnotationFilter(~ gene_biotype == "protein_coding")
    geneFeats = ensembldb::genes(EnsDb, filter = codingFilter, columns=NULL)
    exonFeats = ensembldb::exons(EnsDb, filter = codingFilter, columns=NULL)
    UTR5Feats = ensembldb::fiveUTRsByTranscript(EnsDb, 
                                                filter = codingFilter,
                                                columns = NULL)
    UTR3Feats = ensembldb::threeUTRsByTranscript(EnsDb, 
                                                 filter = codingFilter, 
                                                 columns = NULL)
    UTR5Feats = unlist(UTR5Feats)
    UTR3Feats = unlist(UTR3Feats)
    # Smash 
    geneFeats = reduce(geneFeats)
    exonFeats = reduce(exonFeats)
    UTR5Feats = reduce(UTR5Feats)
    UTR3Feats = reduce(UTR3Feats)
    # Keep only standard chromosomes
    geneFeats= keepStandardChromosomes(geneFeats, pruning.mode = "coarse")
    exonFeats = keepStandardChromosomes(exonFeats, pruning.mode = "coarse")
    UTR5Feats = keepStandardChromosomes(UTR5Feats, pruning.mode = "coarse")
    UTR3Feats = keepStandardChromosomes(UTR3Feats, pruning.mode = "coarse")
    # Since we're storing this data, we want it to be small.
    elementMetadata(geneFeats) = NULL
    elementMetadata(exonFeats) = NULL
    elementMetadata(UTR5Feats) = NULL
    elementMetadata(UTR3Feats) = NULL
    # Change from ensembl-style chrom annotation to UCSC_style
    seqlevels(geneFeats) = paste0("chr", seqlevels(geneFeats))
    seqlevels(exonFeats) = paste0("chr", seqlevels(exonFeats))
    seqlevels(UTR5Feats) = paste0("chr", seqlevels(UTR5Feats))
    seqlevels(UTR3Feats) = paste0("chr", seqlevels(UTR3Feats))
    list(genesGR=geneFeats, exonsGR=exonFeats, threeUTRGR=UTR3Feats, 
         fiveUTRGR=UTR5Feats)
  }, error=function(err){
    # Try a TxDb instead
    message("Failed, trying a UCSC TxDb instead")
    message("Here's the original error message:")
    message(err)
    txdb = loadTxDb(genome)
    exonFeats = GenomicFeatures::exons(txdb)
    geneFeats = GenomicFeatures::transcripts(txdb)
    UTR3Feats = GenomicFeatures::threeUTRsByTranscript(txdb)
    UTR5Feats = GenomicFeatures::fiveUTRsByTranscript(txdb) 
    geneFeats = reduce(geneFeats)
    exonFeats = reduce(exonFeats)
    list(genesGR=geneFeats, exonsGR=exonFeats,  threeUTRGR=UTR3Feats, 
         fiveUTRGR=UTR5Feats)
  })
  return(geneModels)
}

myTSS <- function(genome) {
  feats = tryCatch({
    message("building TSS using Ensembl")
    EnsDb =  buildEnsDb(genome)
    codingFilter = AnnotationFilter::AnnotationFilter(
      ~ gene_biotype == "protein_coding")
    featsWide = ensembldb::genes(EnsDb, filter=codingFilter)
    # Grab just a single base pair at the TSS
    feats = promoters(featsWide, 1, 1)
    # Change from ensembl-style chrom annotation to UCSC_style
    seqlevels(feats) = paste0("chr", seqlevels(feats))
    feats
  }, error=function(err){
    message("failed, trying a UCSC TxDb instead")
    message("Here's the original error message:")
    message(err)
    # Using TxDb
    txdb = loadTxDb(genome)
    # Grab just a single base pair at the TSS
    GenomicFeatures::promoters(txdb, 1, 1)
  })
  return(feats)
}

myChromSizes <- function(genome){
  if (requireNamespace(BSgm, quietly=TRUE)){
    library (BSgm, character.only = TRUE)
    BSG = eval(as.name(BSgm))
  } else {
    library (BSg, character.only = TRUE)
    BSG = eval(as.name(BSg))
  }
  chromSizesGenome = seqlengths(BSG)
  return(chromSizesGenome)
}

plotBoth <- function(plotId, g){
  pth = paste0(opt$outputFolder, "/", fileId, "_", plotId)
  print(paste0("Plotting: ", pth))
  ggplot2::ggsave(paste0(pth, ".png"), g, device="png", width=8, height=8, units="in")
  ggplot2::ggsave(paste0(pth, ".pdf"), g, device="pdf", width=8, height=8, units="in")
}

getPlotReportDF <- function(plotId, title){
  pth = paste0(opt$outputFolder, "/", fileId, "_", plotId)
  rel_pth = getRelativePath(pth, paste0(opt$outputFolder, "/../../../"))
  print(paste0("Writing plot json: ", rel_pth))
  newPlot = data.frame(
    "name"=plotId, 
    "title"=title, 
    "thumbnail_path"=paste0(rel_pth, ".png"), 
    "path"=paste0(rel_pth, ".pdf")
  )
  return(newPlot)
}


doItAall <- function(query, fileId, genome, cellMatrix) {
  plots = data.frame(stringsAsFactors=F)
  bsGenomeAvail = ifelse((requireNamespace(BSg, quietly=TRUE) | requireNamespace(BSgm, quietly=TRUE)), TRUE, FALSE)
  # TSS distance plot
  tryCatch(
    expr = {
      if (genome %in% c("hg19", "hg38", "mm10", "mm9")){
        plotBoth("tssdist", plotFeatureDist( calcFeatureDistRefTSS(query, genome), featureName="TSS"))
      } else{
        if (ensDbGtf == "None") {
          message("Ensembl annotation gtf file not provided. Skipping TSS distance plot ... ")
        } else {
        tss = myTSS(genome)
        plotBoth("tssdist", plotFeatureDist( calcFeatureDist(query, tss), featureName="TSS"))
        }
      }
      plots = rbind(plots, getPlotReportDF("tssdist", "Region-TSS distance distribution"))
      message("Successfully calculated and plot TSS distance.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
  
  
  # Chromosomes region distribution plot
  tryCatch(
    expr = {
      if (genome %in% c("mm39", "dm3", "dm6", "ce10", "ce11", "danRer10", "danRer10", "T2T")){
        chromSizes = myChromSizes(genome)
        genomeBins  = getGenomeBins(chromSizes)
        plotBoth("chrombins", plotChromBins(calcChromBins(query, genomeBins)))
      } else{
        plotBoth("chrombins", plotChromBins(calcChromBinsRef(query, genome)))
      }
      
      plots = rbind(plots, getPlotReportDF("chrombins", "Regions distribution over chromosomes"))
      message("Successfully calculated and plot chromosomes region distribution.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
  
  
  # OPTIONAL: Plot GC content only if proper BSgenome package is installed. 
  if (bsGenomeAvail) {
    tryCatch(
      expr = {
        if (requireNamespace(BSgm, quietly=TRUE)){
          library (BSgm, character.only = TRUE)
          bsg = eval(as.name(BSgm))
          gcvec = calcGCContent(query, bsg)
        } else {
          library (BSg, character.only = TRUE)
          bsg = eval(as.name(BSg))
          gcvec = calcGCContent(query, bsg)
        }
        plotBoth("gccontent", plotGCContent(gcvec))
        plots = rbind(plots, getPlotReportDF("gccontent", "GC content"))
        message("Successfully calculated and plot GC content.")
      },
      error = function(e){
        message('Caught an error!')
        print(e, gcvec)
      }
    ) 
  }
  
  # Partition plots, default to percentages
  tryCatch(
    expr = {
      gp = calcPartitionsRef(query, genome)
      plotBoth("paritions", plotPartitions(gp))
      plots = rbind(plots, getPlotReportDF("paritions", "Regions distribution over genomic partitions"))
      # flatten the result returned by the function above
      partiotionNames = as.vector(gp[,"partition"])
      partitionsList = list()
      for(i in seq_along(partiotionNames)){
        partitionsList[[paste0(partiotionNames[i], "_frequency")]] = 
          as.vector(gp[,"Freq"])[i]
        partitionsList[[paste0(partiotionNames[i], "_percentage")]] = 
          as.vector(gp[,"Freq"])[i]/length(query)	        
      }
      message("Successfully calculated and plot regions distribution over genomic partitions.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
  
  # Expected partition plots
  tryCatch(
    expr = {
      plotBoth("expected_partitions", plotExpectedPartitions(calcExpectedPartitionsRef(query, genome)))
      plots = rbind(plots, getPlotReportDF("expected_partitions", "Expected distribution over genomic partitions"))
      message("Successfully calculated and plot expected distribution over genomic partitions.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
 
  # Cumulative partition plots
  tryCatch(
    expr = {
      plotBoth("cumulative_partitions", plotCumulativePartitions(calcCumulativePartitionsRef(query, genome)))
      plots = rbind(plots, getPlotReportDF("cumulative_partitions", "Cumulative distribution over genomic partitions"))
      message("Successfully calculated and plot cumulative distribution over genomic partitions.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
  
  # QThist plot
  tryCatch(
    expr = {
      widths = calcWidth(query)
      plotBoth("widths_histogram", plotQTHist(widths))
      plots = rbind(plots, getPlotReportDF("widths_histogram", "Quantile-trimmed histogram of widths"))
      message("Successfully calculated and plot quantile-trimmed histogram of widths.")
    },
    error = function(e){
      message('Caught an error!')
      print(e, widths)
    }
  ) 
  
  # Neighbor regions distance plots
  tryCatch(
    expr = {
      plotBoth("neighbor_distances", plotNeighborDist(calcNeighborDist(query)))
      plots = rbind(plots, getPlotReportDF("neighbor_distances", "Distance between neighbor regions"))
      message("Successfully calculated and plot distance between neighbor regions.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    }
  ) 
  
  # Tissue specificity plot if open signal matrix is provided
  if (cellMatrix == "None") {
    message("open signal matrix not provided. Skipping tissue specificity plot ... ")
  } else {
    tryCatch(
      expr = {
        plotBoth("open_chromatin", plotSummarySignal(calcSummarySignal(query, data.table::fread(cellMatrix))))
        plots = rbind(plots, getPlotReportDF("open_chromatin", "Cell specific enrichment for open chromatin"))
        message("Successfully calculated and plot cell specific enrichment for open chromatin.")
      },
      error = function(e){
        message('Caught an error!')
        print(e)
      }
    ) 
  }
  
  # Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/schemas/bedfiles_schema.yaml
  bedmeta = list(
    name=fileId,
    regions_no=length(query),
    mean_region_width=ifelse(exists('widths'), signif(mean(widths), digits = 4), NA),
    md5sum=opt$digest
  )
  if (exists('gcvec') && !isEmpty(gcvec)){
    gc_content <- list(gc_content = signif(mean(gcvec), digits = 4))
    bedmeta = append(bedmeta, gc_content)
  }
  if (exists('TSSdist') && !all(is.na(TSSdist))){
    tss <- list(median_TSS_dist = signif(median(abs(TSSdist), na.rm=TRUE), digits = 4))
    bedmeta = append(bedmeta, tss)
  }
  if (exists('partitionsList')){
    write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), paste0(outfolder, "/", fileId, ".json"))
  } else {
     write(jsonlite::toJSON(c(bedmeta), pretty=TRUE), paste0(outfolder, "/", fileId, ".json"))
    }
  
  if (exists('plots')){
    write(jsonlite::toJSON(plots, pretty=TRUE), paste0(outfolder, "/", fileId, "_plots.json"))
  }
}

# define values and output folder for doitall()
fileId = opt$fileId
bedPath = opt$bedfilePath
outfolder = opt$outputFolder
genome = opt$genome
cellMatrix = opt$openSignalMatrix
ensDbGtf = opt$ensDbGtf


# build BSgenome package ID to check whether it's installed
if (genome == "T2T"){
  BSg = "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
} else {
  if (startsWith(genome, "hg") | startsWith(genome, "grch")) {
  orgName = "Hsapiens"
  } else if (startsWith(genome, "mm") | startsWith(genome, "grcm")){
  orgName = "Mmusculus"
  } else if (startsWith(genome, "dm")){
  orgName = "Dmelanogaster"
  } else if (startsWith(genome, "ce")){
  orgName = "Celegans"
  } else if (startsWith(genome, "danRer")){
  orgName = "Drerio"
  } 
  BSg = paste0("BSgenome.", orgName , ".UCSC.", genome)
}

BSgm = paste0(BSg, ".masked")

# read bed file and run doitall()
query = LOLA::readBed(bedPath)
doItAall(query, fileId, genome, cellMatrix)
