[synapse - getting started]:https://www.synapse.org/#!Synapse:syn1834618/WIKI/56017
```
# If not installed yet
install.packages("devtools")
source('http://depot.sagebase.org/CRAN.R')
pkgInstall("synapseClient")
```

```
require("devtools")
require(synapseClient)
install_github("rGithubClient", "brian-bot")
require("rGithubClient")
```

Use the synapseLogin() function to log into synapse:
```synapseLogin("me@myemail.com", "mypwd")```  
or visit  [synapse - getting started] to configure your .Rprofile for an automatic synapse login.


```
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
