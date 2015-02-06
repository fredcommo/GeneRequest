\name{geneRequest-package}
\alias{geneRequest-package}

\docType{package}
\title{
  NCBI Query for Gene Annotation.
}
\description{
  Retrieve gene annotations from NCBI, given a HUGO symbol or an entrezgene id, using web requests.
}
\details{
  \tabular{ll}{
    Package: \tab geneRequest\cr
    Type: \tab Package\cr
    Version: \tab 1.0\cr
    Date: \tab 2015-02-06\cr
    License: \tab GPL\cr
  }
  This function simply constructs an appropriate link to query the NCBI databases.
  Information is extracted from The xml file returned.
}
\author{
Frederic Commo <frederic.commo@gustaveroussy.fr>
}

\examples{
geneRequest("ercc1")
geneRequest(2067)
}
