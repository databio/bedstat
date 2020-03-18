.install_pkg = function(p, bioc=FALSE) {
    if(!require(package = p, character.only=TRUE)) {
        if(bioc) {
            BiocManager::install(pkgs = p)
        } else {
            install.packages(pkgs = p)   
        }
    }
}

.install_pkg("BiocManager")
.install_pkg("optparse")
.install_pkg("devtools")
devtools::install_github("databio/GenomicDistributions", ref="dev")
genomes = list(Hsapiens = c("hg18","hg19","hg38"), 
               Mmusculus = c("mm10","mm9"))
for(name in names(genomes)) {
    for(genome in genomes[[name]]) {
        # should install non-masked too
        .install_pkg(p=paste0("BSgenome.", name, 
                                ".UCSC.", genome,".masked"), 
                     bioc=TRUE)
    }
}
.install_pkg("GenomicRanges", bioc=TRUE)
.install_pkg("LOLA", bioc=TRUE)
