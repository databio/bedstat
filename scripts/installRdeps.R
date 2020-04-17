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
devtools::install_github("databio/GenomicDistributions")
.install_pkg("GenomicRanges", bioc=TRUE)
.install_pkg("LOLA", bioc=TRUE)
