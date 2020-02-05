install.packages("BiocManager"); library(BiocManager)
install.packages("optparse")
install.packages("devtools")
BiocManager::install("GenomicRanges")
devtools::install_github("databio/GenomicDistributions")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19") # depending on the genome used
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19.masked")
BiocManager::install("LOLA")