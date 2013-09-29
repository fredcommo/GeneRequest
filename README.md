```
# If not installed yet
install.packages("devtools")
install_github("rGithubClient", "brian_bot")
source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")
```

```
require("devtools")
require("rGithubClient")
require(synapseClient)

# Use the synapseLogin() function: synapseLogin("me@myemail.com", "mypwd")
# or visit https://www.synapse.org/#!Synapse:syn1834618/WIKI/56017 to configure your .Rprofile

getFilesList <- function(git, tag = ''){
  flist <- git@tree$path
  return(flist[grep(tag, flist)])
}

git <- getRepo('fredcommo/GeneRequest')
Rlist <- getFilesList(git, '.R')

# source GeneRequest.v6.4
sourceRepoFile(git, Rlist[5])

geneList <- c("egfr", "fgfr1", "IwantRedHair")
myList <- lapply(geneList, function(gene) GeneRequest.v6.4(gene))
as.data.frame(do.call(rbind, myList))
```

**- version v6.3 takes a vector of genes (symbols or entrezids) as argument and use a for loop.**

**- version v6.4 takes only one gene at a time - combined with lapply(), this version works faster.**

**- version 7 searches in a table built by one of the versions above. I have to rename this guy!**
