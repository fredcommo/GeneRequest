#################################################################################################

GeneRequest.v6.4 <- function(geneId, DB = "gene", HG = hg19, bySymb = TRUE, verbose = TRUE){
  ow <- options("warn")
  options(warn = -1)
  require(XML)
  cumLen <- cumsum(as.numeric(HG$length))
  cumLen <- c(0, cumLen[1:23])
	Sys.sleep(0.01)

  entrezId <- NA
  Symbol <- "Not found"
  out <- .default()
  ids <- geneId
    
  if(bySymb){
    query <- geneId
    ids <- unlist(gsearch(paste0(toupper(query), "%5Bsymbol%5D%20homo%20sapiens"), DB = DB, kTries = 10)) #homo sapiens
		if(is.null(ids) & verbose) cat("\n ***", query, "... Can't find this guy: Stop kidding *** !\n\n")
    }
  if(!is.null(ids)){
    if(verbose) cat(query, "found:", length(ids), "ids...\t")
  	j = 1
	  Id <- ids[j]
	  geneSummary <- unlist(gsummary(Id, DB = DB, kTries = 10))					#homo sapiens
	  if(length(ids)>1)
		  while (!.checkSummary(query, geneSummary, cumLen) & j < length(ids)){
        j = j + 1
				Id <- ids[j]
				geneSummary <- unlist(gsummary(Id, DB = DB, kTries = 10))			#homo sapiens
				}
    if(.checkSummary(query, geneSummary, cumLen)){
      out <- .getSummary(geneSummary, cumLen)
      entrezId <- Id
      }
    if(verbose) cat('Done.\n')
  	}
    options(ow)
	return(c("query" = query, out, "entrezgeneId" = as.numeric(entrezId)))
}

##################################
# Helper Functions

gsearch <- function (geneSymb, DB, kTries){
  
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

gsummary <- function (id, DB, kTries){

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
  sapply(c("//Item"), xpathApply, doc = doc, fun = xmlValue)
}
.inList <- function(query, geneTable){
  if(is.null(nrow(geneTable)))
    return(FALSE)
  query <- sprintf("^%s$", toupper(query))
  return(any(grepl(query, geneTable$symbol)))
}
.getChr <- function(geneSummary){
  chr <- geneSummary[6]
  if(chr == "") chr <- NA
  if(!is.na(chr) & chr == "X") chr <- 23
  if(!is.na(chr) & chr == "Y") chr <- 24
  if(!is.na(chr) & chr == "X, Y") chr <- 23
  return(as.numeric(chr))
}
.getChrLoc <- function(geneSummary){
  if(length(geneSummary)!=27){
    rangeId <- geneSummary[19]
    chrStart <- geneSummary[20]
    chrEnd <- geneSummary[21]
  }
  else{
    rangeId <- geneSummary[20]
    chrStart <- geneSummary[21]
    chrEnd <- geneSummary[22]
  }
  return(cbind.data.frame(rangeId = rangeId, chrStart = as.numeric(chrStart), chrEnd = as.numeric(chrEnd)))
}
.genomeLoc <- function(Chr, loc, chrLen){
  return(as.numeric(loc) + chrLen[Chr])
}
.getSummary <- function(geneSummary, cumLen){
  symbol <- geneSummary[1]
  chr <- .getChr(geneSummary)
  geneLoc <- .getChrLoc(geneSummary)
#  if(!bySymb) cat(as.character(Id), "found as", Symbol, "\n")
  return(list(symbol = symbol,
              fullName = geneSummary[2],
              org = geneSummary[3],
              verifId = as.numeric(geneSummary[5]),
              chr = chr,
              cytoband = geneSummary[8],
              alias = geneSummary[9],
              description = geneSummary[10],
              chrStart = geneLoc$chrStart,
              chrEnd = geneLoc$chrEnd,
              rangeId = geneLoc$rangeId,
              genomStart = .genomeLoc(chr, geneLoc$chrStart, cumLen),
              genomEnd = .genomeLoc(chr, geneLoc$chrEnd, cumLen)
              #status = geneSummary[2],
              )
  )
}
.checkSummary <- function(query, geneSummary, cumLen){
  
  tmp <- .getSummary(geneSummary, cumLen)
  alias <- unlist(strsplit(tmp$alias, ', '))
  correctSymbol <- toupper(query) %in% c(tmp$symbol, alias)
  return(# toupper(query) == tmp$Symbol |
           tmp$org == "Homo sapiens" |
           tmp$fullName != "reserved" |
           tmp$verifId == 0 |
           !tmp$chr %in% c("", NA))
}
.default <- function(){
  return(list(symbol = NA, fullName = NA, org = NA, verifId = NA,
              chr = NA, cytoband = NA, alias = NA, description = NA,
              chrStart = NA, chrEnd = NA, rangeId = NA, genomStart = NA, genomEnd = NA)
    )
}
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# End functions
