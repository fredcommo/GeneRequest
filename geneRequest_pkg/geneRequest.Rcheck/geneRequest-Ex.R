pkgname <- "geneRequest"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('geneRequest')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("geneRequest-package")
### * geneRequest-package

flush(stderr()); flush(stdout())

### Name: geneRequest-package
### Title: NCBI Query for Gene Annotation.
### Aliases: geneRequest-package

### ** Examples

geneRequest("ercc1")
geneRequest(2067)



cleanEx()
nameEx("geneRequest")
### * geneRequest

flush(stderr()); flush(stdout())

### Name: geneRequest
### Title: Gene annotations via Web Queries
### Aliases: geneRequest

### ** Examples

# Simple query using symbol or entrezgeneId
geneRequest("ercc1")
geneRequest(2067, bySymb=FALSE)

# Multiple queries
ids <- c("egfr", "erbb2", "fgfr1")
annots <- lapply(ids, function(id) geneRequest(id) )
annots <- do.call(rbind, annots)
annots



cleanEx()
nameEx("hg19")
### * hg19

flush(stderr()); flush(stdout())

### Name: hg19
### Title: HG19 chromosomes length
### Aliases: hg19
### Keywords: datasets

### ** Examples

  load(system.file("extdata", "hg19.rda", package="geneRequest"))



### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
