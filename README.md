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

git <- getRepo('fredcommo/IC50')
Rlist <- getFilesList(git, '.R')

# source GeneRequest.v6.4
sourceRepoFile(git, Rlist[5])

geneList <- c("egfr", "fgfr1", "greenBanana")
GeneRequest.v6.4(geneList)
```
