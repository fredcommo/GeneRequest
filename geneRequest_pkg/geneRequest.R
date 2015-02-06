#################################################################################################
# if(!exists("HG")){
#   require(synapseClient)
#   cat("Loading hg19 from synapse...")
#   # use hg19 chromosome size here.
#   ent <- synGet("syn2141399")
#   HG <- read.csv(ent@filePath, header = TRUE, sep = "\t")
#   cat("Done.\n")
# }

# require(XML)

geneRequest <- function(geneId, DB = "gene", bySymb = TRUE, verbose = TRUE){
  ow <- options("warn")
  options(warn = -1)
  cumLen <- cumsum(as.numeric(HG$length))
  cumLen <- c(0, cumLen[1:23])
	Sys.sleep(0.01)

  entrezId <- NA
  ids <- geneId

  if(bySymb){
    query <- toupper(geneId)
    ids <- unlist(.gsearch(paste0(query, "%5Bsymbol%5D%20homo%20sapiens"), DB = DB, kTries = 10))
		if(is.null(ids)){
      if(verbose){
        txt <- sprintf("\n*** Can't find \"%s\". Stop kidding! ***\n", query)
        cat(txt)
        }
      return(NULL)
		  }
    }

  if(!is.null(ids)){
    if(verbose)
      cat(query, "found:", length(ids), "id(s)...\t")
  	j = 1
	  id <- ids[j]
	  geneSummary <- .gsummary(id, DB = DB, kTries = 10)
	  if(length(ids)>1)
		  while (!.checkSummary(query, geneSummary) & j < length(ids)){
        j = j + 1
				id <- ids[j]
				geneSummary <- .gsummary(id, DB = DB, kTries = 10)
				}
    if(.checkSummary(query, geneSummary)){
      out <- .getSummary(geneSummary, cumLen)
      entrezId <- as.numeric(id)
      if(verbose)
        cat('Done.\n')
    	return(
        c("query" = query, out, "entrezgeneId" = as.numeric(entrezId))
        )
      }
  	}
    return(NULL)
    options(ow)
}

##################################
# Helper Functions

.gsearch <- function (geneSymb, DB, kTries){
  
  # geneSymb : use official symbols. Multiple requests are accepted, e.g. "EGFR, Homo sapiens"
  # database : have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
  # 		on available databases and other e-tools as well.
  # ! This function can return more than one Id !
   
  gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  gsrch.mode <- paste("db=", DB, "&retmode=xml","&term=", sep = "")
  URL <- paste(gsrch.stem, gsrch.mode, geneSymb, sep = "")
  k = 1
  doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
  while(class(doc)[1] == 'try-error' & k <= kTries){
  	doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
	k = k + 1
	cat("Connexion error during gsearch - I'm trying again...", k, "\n")
  }
  if(class(doc)[1] == 'try-error') stop ('Connection temporarily unavailable. Check it manually:\n', URL)
  sapply(c("//Id"), xpathApply, doc = doc, fun = xmlValue)
}

.gsummary <- function (id, DB, kTries){

  # id is provided by gsearch()
  # database : let's have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
  # 		on available databases and other e-tools as well.  

  sum.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
  sum.mode <- paste("db=", DB, "&id=", sep = "")
  k = 1
  doc <- try(xmlTreeParse(paste(sum.stem, sum.mode, id, sep = ""), isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
  while(class(doc)[1] == 'try-error' & k < kTries){
    cat("Connexion error during gsummary - I'm trying again...", k, "\n")
    k = k + 1
    doc <- try(xmlTreeParse(paste(sum.stem, sum.mode, id, sep = ""), isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
  }
  items <- unlist(sapply(c("//Item"), xpathApply, doc = doc, fun = xmlGetAttr, name="Name")[1:30])
  values <- unlist(sapply(c("//Item"), xpathApply, doc = doc, fun = xmlValue)[1:30])
  names(values) <- items
  return( values )
}
.getItem <- function(dat, item=NULL){
  if(item %in% names(dat))
    return( as.character(dat[item]) )
  return( NA )
}
# .inList <- function(query, geneTable){
#   if(is.null(nrow(geneTable)))
#     return(FALSE)
#   query <- sprintf("^%s$", toupper(query))
#   return(any(grepl(query, geneTable$symbol)))
# }
.getChr <- function(geneSummary){
  chr <- .getItem(geneSummary, "Chromosome")
  if(chr == "") chr <- NA
  if(!is.na(chr) & chr == "X") chr <- 23
  if(!is.na(chr) & chr == "Y") chr <- 24
  if(!is.na(chr) & chr == "X, Y") chr <- 23
  return(as.numeric(chr))
}
.getChrLoc <- function(geneSummary){
    chrAccVer <- .getItem(geneSummary, "ChrAccVer")
    chrStart <- as.numeric(.getItem(geneSummary, "ChrStart"))
    chrEnd <-  as.numeric(.getItem(geneSummary, "ChrStop"))
  return(cbind.data.frame(chrAccVer = chrAccVer, chrStart = chrStart, chrEnd = chrEnd))
}
.getSummary <- function(geneSummary, cumLen){
  chr <- .getChr(geneSummary)
  geneLoc <- .getChrLoc(geneSummary)
  
  return(c(symbol = .getItem(geneSummary, "Name"),
              fullName = .getItem(geneSummary, "NomenclatureName"),
              alias = .getItem(geneSummary, "OtherAliases"),
              organism = .getItem(geneSummary, "ScientificName"),
              verifId = .getItem(geneSummary, "Status"),
              status = .getItem(geneSummary, "NomenclatureStatus"),
              chr = chr,
              cytoband = .getItem(geneSummary, "MapLocation"),
              exonCount = .getItem(geneSummary, "ExonCount"),
              accVersion = as.character(geneLoc$chrAccVer),
              chrStart = geneLoc$chrStart,
              chrEnd = geneLoc$chrEnd,
              genomStart = geneLoc$chrStart + cumLen[chr],
              genomEnd = geneLoc$chrEnd + cumLen[chr]
           )
         )
}
.checkSummary <- function(query, geneSummary){
  symbol <- .getItem(geneSummary, "Name")
  alias <- .getItem(geneSummary, "OtherAliases")
  alias <- unlist(strsplit(alias, ', '))
  correctSymbol <- toupper(query) %in% c(symbol, alias)
  official <- .getItem(geneSummary, "NomenclatureStatus") == "Official"
  return(correctSymbol && official)
}
#   organism <- .getItem(geneSummary, "ScientificName")
#   return(organism == "Homo sapiens" &&
#            .getItem(geneSummary, "NomenclatureStatus") == "Official" #|
#          )
#}
# .default <- function(){
#   return(list(symbol = NA, fullName = NA, org = NA, verifId = NA,
#               chr = NA, cytoband = NA, alias = NA, description = NA,
#               chrStart = NA, chrEnd = NA, rangeId = NA, genomStart = NA, genomEnd = NA)
#     )
# }
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# End functions
