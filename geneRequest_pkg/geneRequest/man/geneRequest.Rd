\name{geneRequest}
\alias{geneRequest}

\title{
  Gene annotations via Web Queries
}
\description{
  This function simply constructs an appropriate link to query the NCBI databases.
  Information is extracted from The xml file returned.
}
\usage{
  geneRequest(geneId, DB = "gene", bySymb = TRUE, verbose = TRUE)
}

\arguments{
  \item{geneId}{
: character or numeric. \code{geneId} can be either a valid HUGO symbol (not case sensitive), or an entrezgene id. See \code{note}.
}
  \item{DB}{
: string. The NCBI database to query. This must be one of the NCBI databases supporting web queries. 
}
  \item{bySymb}{
: logical. Must be set to \code{TRUE} when a HUGO symbols is used.
}
  \item{verbose}{
: logical. If \code{TRUE}, progress in prompted.
}
}
\value{
  A vector of values.  
}
\author{
  Frederic Commo
}
\note{
See \href{http://www.ncbi.nlm.nih.gov/books/NBK25501/}{Introduction to the E-utilities} for more details.
}

\examples{
# Simple query using symbol or entrezgeneId
geneRequest("ercc1")
geneRequest(2067, bySymb=FALSE)

# Multiple queries
ids <- c("egfr", "erbb2", "fgfr1")
annots <- lapply(ids, function(id) geneRequest(id) )
annots <- do.call(rbind, annots)
annots
}