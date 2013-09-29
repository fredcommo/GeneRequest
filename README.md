```
# If not installed yet
install.packages("devtools")
install_github("rGithubClient", "brian_bot")
```

```
require("devtools")
require("rGithubClient")
require(RCurl)
require(foreign)
url <- "https://raw.github.com/fredcommo/GeneRequest/master/human.chrom.info.hg19.FC.txt"
hg19 <- read.csv(textConnection(getURL(url)), header = T, sep = '\t')

getFilesList <- function(git, tag = ''){
  flist <- git@tree$path
  return(flist[grep(tag, flist)])
}

git <- getRepo('fredcommo/GeneRequest')
Rlist <- getFilesList(git, '.R')

# source GeneRequest.v6.4
sourceRepoFile(git, Rlist[5])

geneList <- c("egfr", "fgfr1", "greenBanana")
myList <- lapply(geneList, function(gene) GeneRequest.v6.4(gene))
as.data.frame(do.call(rbind, myList))
```

** version v6.3 takes a list of genes as argument and use a for loop.**

** version v6.4 takes only one gene at a time - combined with lapply(), this version works faster.**

** version 7 searches in a table built by one of the versions above. I have to rename this guy!**
