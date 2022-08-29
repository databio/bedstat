library(GenomicDistributions)
library(GenomicDistributionsData)
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
              help="genome reference to calculate against", metavar="character"))

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
      TSSdist = calcFeatureDistRefTSS(query, genome)
      plotBoth("tssdist", plotFeatureDist(TSSdist, featureName="TSS"))
      plots = rbind(plots, getPlotReportDF("tssdist", "Region-TSS distance distribution"))
      message("Successfully calculated and plot TSS distance.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    }
  ) 
  
  
  # Chromosomes region distribution plot
  tryCatch(
    expr = {
      plotBoth("chrombins", plotChromBins(calcChromBinsRef(query, genome)))
      plots = rbind(plots, getPlotReportDF("chrombins", "Regions distribution over chromosomes"))
      message("Successfully calculated and plot chromosomes region distribution.")
    },
    error = function(e){
      message('Caught an error!')
      print(e)
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
    }
  ) 
  
  
  # OPTIONAL: Plot GC content only if proper BSgenome package is installed. 
  if (bsGenomeAvail) {
    tryCatch(
      expr = {
        gcvec = calcGCContentRef(query, genome)
        plotBoth("gccontent", plotGCContent(gcvec))
        plots = rbind(plots, getPlotReportDF("gccontent", "GC content"))
        message("Successfully calculated and plot GC content.")
      },
      error = function(e){
        message('Caught an error!')
        print(e, gcvec)
      },
      warning = function(w){
        message('Caught an warning!')
        print(w)
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
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
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
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
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
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
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
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
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
    },
    warning = function(w){
      message('Caught an warning!')
      print(w)
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
      },
      warning = function(w){
        message('Caught an warning!')
        print(w)
      }
    ) 
  }
  
  # Note: names of the list elements MUST match what's defined in: https://github.com/databio/bbconf/blob/master/bbconf/schemas/bedfiles_schema.yaml
  bedmeta = list(
    name=fileId,
    gc_content=ifelse(bsGenomeAvail, signif(mean(gcvec), digits = 4), NA),
    regions_no=length(query),
    median_TSS_dist=ifelse(exists('TSSdist'), signif(median(abs(TSSdist), na.rm=TRUE), digits = 4), NA),
    mean_region_width=ifelse(exists('widths'), signif(mean(widths), digits = 4), NA),
    md5sum=opt$digest
  )
  if (exists('partitionsList')){
    write(jsonlite::toJSON(c(bedmeta, partitionsList), pretty=TRUE), paste0(outfolder, "/", fileId, ".json"))
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
orgName = "Mmusculus"

# build BSgenome package ID to check whether it's installed
if (startsWith(genome, "hg") | startsWith(genome, "grch")) orgName = "Hsapiens"

BSg = paste0("BSgenome.", orgName , ".UCSC.", genome)
BSgm = paste0(BSg, ".masked")

# read bed file and run doitall()
query = LOLA::readBed(bedPath)
doItAall(query, fileId, genome, cellMatrix)