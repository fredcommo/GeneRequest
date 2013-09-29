#################################################################################################

# Exemple liste SAFIR

GeneRequest.v6.3 <- function(geneList, DB = "gene", bySymb = TRUE, verbose = TRUE, connectTries = 15){
  require(synapseClient)
  require(XML)

  ent <- synGet('syn2141399')
  hg19 <- read.csv(ent@filePath, header = TRUE, sep = '\t')
  cumLen <- cumsum(as.numeric(hg19$length))
  cumLen <- c(0, cumLen[1:23])

  geneList <- as.character(geneList)
  Total = length(geneList)
	# pb <- tkProgressBar(title = "The R'atWork BaBar", min = 0, max = Total, width = 500)

	# cum.len <- cumsum(hg19.info$length)

  geneTable <- c()
  
	for(i in 1:Total){
		Sys.sleep(0.1)

		# launch & increment the pBar
	  # setTkProgressBar(pb, i, label = paste("Take it easy... I'm workin' for U... :p ", round(i/Total*100, 2), "% done"))

		entrezId <- NA
    Symbol <- "Not found"
    out <- .default()
    ids <- geneList[i]
    
    if(bySymb){
      query <- geneList[i]
      ids <- unlist(gsearch(paste0(toupper(query), "%5Bsymbol%5D%20homo%20sapiens"), DB = DB, kTries = connectTries))					#homo sapiens
		  if(is.null(ids)) cat("\n ***", query, "... Can't find this guy: Stop kidding *** !\n\n")
      }

    if(.isAlias(query, geneTable[,c(2,8,15)]))
      cat(query, 'already in table: probable aliases.\n')
    else{
      if(!is.null(ids)){
        if(verbose) cat(query, "found:", length(ids), "id(s).\n")
    		j = 1
	    	Id <- ids[j]
	    	geneSummary <- unlist(gsummary(Id, DB = DB, kTries = connectTries))					#homo sapiens
	  	  if(length(ids)>1)
		  	  while (!.checkSummary(query, geneSummary, cumLen) & j < length(ids)){
			  	  j = j + 1
				    Id <- ids[j]
				    geneSummary <- unlist(gsummary(Id, DB = DB, kTries = connectTries))			#homo sapiens
				    }
        if(.checkSummary(query, geneSummary, cumLen)){
          out <- .getSummary(geneSummary, cumLen)
          entrezId <- Id
          }
  		  }
		  geneTable <- rbind(geneTable, cbind(query = query, as.data.frame(out), entrezgeneId = entrezId))
      }
	  }
  rownames(geneTable) <- seq(1, nrow(geneTable))
	return(geneTable)
}

##################################
# Helper Functions

gsearch <- function (geneSymb, DB, kTries){
  
  # geneSymb : use official symbols. Multiple requests are accepted, e.g. "EGFR, Homo sapiens"
  # database : have a look at 'http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html' for details
  #             on available databases and other e-tools as well.
  # ! This function can return more than one Id !
  
  gsrch.stem <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
  gsrch.mode <- paste("db=", DB, "&retmode=xml","&term=", sep = "")
  URL <- paste0(gsrch.stem, gsrch.mode, geneSymb)
  k = 1
  doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
  while(class(doc)[1] == 'try-error' & k <= kTries){
    cat("Connexion error during gsearch() - I'm trying again...", k, "\n")
    doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
    k = k + 1
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
  URL <- paste0(sum.stem, sum.mode, id)
  k = 1
  doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
  while(class(doc)[1] == 'try-error' & k <= kTries){
    cat("Connexion error during gsummary() - I'm trying again...", k, "\n")
    doc <- try(xmlTreeParse(URL, isURL = TRUE, useInternalNodes = TRUE), silent = TRUE)
    k = k + 1
  }
  if(class(doc)[1] == 'try-error') stop ('Connection temporarily unavailable. Check it manually:\n', URL)
  sapply(c("//Item"), xpathApply, doc = doc, fun = xmlValue)
}
.isAlias <- function(query, Table){
  if(is.null(ncol(Table)))
    return(FALSE)
  out <- sapply(1:ncol(Table), function(i) any(grepl(paste0('^',query,'$'), Table[,i])))
  return(any(out == TRUE))
}
.getChr <- function(geneSummary){
  Chr <- geneSummary[6]
  if(Chr == "") Chr <- NA
  if(!is.na(Chr) & Chr == "X") Chr <- 23
  if(!is.na(Chr) & Chr == "Y") Chr <- 24
  if(!is.na(Chr) & Chr == "X, Y") Chr <- 23
  return(as.numeric(Chr))
}
.getChrLoc <- function(geneSummary){
  if(length(geneSummary)!=27){
    RangeId <- geneSummary[19]
    ChrStart <- geneSummary[20]
    ChrStop <- geneSummary[21]
  }
  else{
    RangeId <- geneSummary[20]
    ChrStart <- geneSummary[21]
    ChrStop <- geneSummary[22]
  }
  return(cbind.data.frame(RangeId = RangeId, ChrStart = as.numeric(ChrStart), ChrStop = as.numeric(ChrStop)))
}
.genomeLoc <- function(Chr, loc, chrLen){
  return(as.numeric(loc) + chrLen[Chr])
}
.getSummary <- function(geneSummary, cumLen){
  Symbol <- geneSummary[1]
  Chr <- .getChr(geneSummary)
  geneLoc <- .getChrLoc(geneSummary)
#  if(!bySymb) cat(as.character(Id), "found as", Symbol, "\n")
  return(list(Symbol = Symbol,
              FullName = geneSummary[2],
              Org = geneSummary[3],
              verifId = as.numeric(geneSummary[5]),
              Chr = Chr,
              Cytoband = geneSummary[8],
              Alias = geneSummary[9],
              Description = geneSummary[10],
              ChrStart = geneLoc$ChrStart,
              ChrStop = geneLoc$ChrStop,
              RangeId = geneLoc$RangeId,
              genomicStart = .genomeLoc(Chr, geneLoc$ChrStart, cumLen),
              genomicStop = .genomeLoc(Chr, geneLoc$ChrStop, cumLen)
              #status = geneSummary[2],
              )
  )
}
.checkSummary <- function(query, geneSummary, cumLen){
  tmp <- .getSummary(geneSummary, cumLen)
  Alias <- unlist(strsplit(tmp$Alias, ', '))
  correctSymbol <- toupper(query) %in% c(tmp$Symbol, Alias)
  return( toupper(query) == tmp$Symbol |
           tmp$Org == "Homo sapiens" |
           tmp$FullName != "reserved" |
           tmp$verifId == 0 |
           !tmp$Chr %in% c("", NA))
}
.default <- function(){
  return(list(Symbol = NA, FullName = NA, Org = NA, verifId = NA,
              Chr = NA, Cytoband = NA, Alias = NA, Description = NA,
              ChrStart = NA, ChrStop = NA, RangeId = NA, genomicStart = NA, genomicStop = NA)
    )
}
# http://eutils.ncbi.nlm.nih.gov/entrez/eutils/egquery.fcgi

# End functions