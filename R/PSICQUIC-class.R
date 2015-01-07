.PSICQUIC <- setClass("PSICQUIC",
                      slots=c(df="DataFrame",
                              providers="character"))

setGeneric("interactions", signature="object",
             function(object,
                      id=NA,
                      species=NA,
                      speciesExclusive=TRUE,
                      type=NA,
                      provider=NA,
                      detectionMethod=NA,
                      publicationID=NA,
                      quiet=TRUE)
           standardGeneric ("interactions"))

setGeneric("count", signature="object", function(object)
           standardGeneric ("count"))

setGeneric("providers", signature="object", function(object)
           standardGeneric ("providers"))

setGeneric("providerUrl", signature="object", function(object, provider)
           standardGeneric ("providerUrl"))

setGeneric("rawQuery", signature="object", function(object, provider, rawArgs)
           standardGeneric ("rawQuery"))
#-------------------------------------------------------------------------------
.printf <- function(...) print(noquote(sprintf(...)))
#-------------------------------------------------------------------------------
setValidity("PSICQUIC", function(object) {
    msg = NULL

    if (length (object) ==  0)
       return (TRUE)

    df <- object@df
    
    msg2 <- validObject(df)
    if(!msg2)
        msg <- c(msg, msg2)
    
    if(!all(colnames(df) %in% c("url")))
        msg = c(msg, "must have column named 'url'")

    if (is.null(msg)) TRUE else msg
    
})
#-------------------------------------------------------------------------------
PSICQUIC <- function(test=FALSE)
{
   object <- .PSICQUIC()
   registry.tbl <- .loadRegistry()

   if(all(is.na(registry.tbl)))
       return(NA)
   
   object@df <- .loadRegistry()
   providers <- rownames(object@df)
   
   deleters <- grep("GeneMANIA", providers, ignore.case=TRUE)
       # GeneMANIA too voluminous, too slow, often simply does not work
   if(length(deleters) > 0)
      providers <- providers[-deleters]

   if(test){ # to test robustness against missing providers:
      count <- length(providers)
      selected <- sample(1:count, 3)
      providers <- providers[selected]
      .printf("test providers: %s", paste(providers, collapse=","))
      }

   object@providers <- providers

   object

} # PSICQUIC ctor
#-------------------------------------------------------------------------------
setMethod ("count", signature=c(object="PSICQUIC"),

   function (object) {
      providers[providers!="GeneMANIA"]   # too voluminous, too slow
      })

#-------------------------------------------------------------------------------
setMethod ("providers", signature=c(object="PSICQUIC"),

   function (object) {
      object@providers
      })

#-------------------------------------------------------------------------------
setMethod ("providerUrl", signature=c(object="PSICQUIC"),

   function (object, provider) {
      df <- object@df
      stopifnot(provider %in% rownames(df))
      df[provider, "url"]
      })

#-------------------------------------------------------------------------------
.loadRegistry <- function()
{
    url <- "http://www.ebi.ac.uk/Tools/webservices/psicquic/registry/registry?action=ACTIVE&format=txt"

    if(!url.exists(url)){
       warning(sprintf("failed to get response from PSICQUIC registry at %s", url))
       return(NA)
      }

    
    emptyResult <- DataFrame(url=vector(mode="character", length=0),
                             row.names=character(0))

    tryCatch({
        txt <- getURL(url)
        if(nchar(txt) == 0)
           return (emptyResult)
        lines <- strsplit(txt, "\n")[[1]]
        tokenized.lines <- strsplit(lines, "=")
        providers <- sapply(tokenized.lines, "[", 1)
        urls <- sapply(tokenized.lines, "[", 2)
        #printf(".loadRegistry, providers: %s", paste(providers, collapse=", "))
        return(DataFrame(row.names=providers, url=urls))
        }, error=function(err) {
           print(sprintf("%s, %s", str(err), "server not responding"))
           print(sprintf("failed url: %s", url))
           return(emptyResult)
       })

} # .loadRegistry
#-------------------------------------------------------------------------------
setMethod('show', 'PSICQUIC',

    function(object) {
      providers <- providers(object)
      count <- length(providers)
      msg = sprintf("PSICQUIC object with %d currently responding providers", count)
      if(count == 0)
          return()
      
      cat(msg, ":\n\n", sep="")
      for(i in 1:count){
          cat(sprintf("%15s: %s\n", providers[i], object@df[i,]))
          } # for i
    })
#-------------------------------------------------------------------------------
setMethod ("interactions", "PSICQUIC",

   function (object, id, species, speciesExclusive, type,
             provider, detectionMethod, publicationID, quiet) {

       if(any(is.na(species)))   # can't be exclusive if none defined.
           speciesExclusive <- FALSE
       
       requested.providers = provider
       if(length(requested.providers) == 1 && is.na(requested.providers))
            requested.providers <- providers(object)

       actual.providers <- intersect(providers(object), requested.providers)

          # we discourage GeneMANIA as a default provider, but
          # permit for those asking for it explicilty

       if("GeneMANIA" %in% requested.providers &&
          "GeneMANIA" %in% rownames(object@df))
          actual.providers <- c(actual.providers, "GeneMANIA")

       unrecognized.providers <- setdiff(requested.providers, actual.providers)
       if(length(unrecognized.providers) > 0){
          stop(sprintf("Error! these providers not recognized: %s",
                          paste(unrecognized.providers, collapse=",")))
          }

       result <- data.frame(stringsAsFactors=FALSE)
       pairs <- .enumerateSearchPairs(id)
       pair.count <- length(pairs$a)
       
       if(pair.count > 1) { # warn user that possibly very many separate queries are needed
         s <- sprintf("About to execute %d PSICQUIC queries (%d gene pairings x %d provider/s)",
                      pair.count * length(actual.providers), pair.count,
                      length(actual.providers))
         message(s)
         }
       
       for(provider in actual.providers){
           if(!quiet) .printf("retrieving from %s", provider)
           df <- object@df
           base.url <- df[provider, "url"]
           for(p in 1:pair.count){
               a <- pairs$a[p]
               b <- pairs$b[p]
               tbl <- .runQuery(base.url, a, b, species, speciesExclusive,
                                type, detectionMethod,
                                publicationID, quiet)
               if(!ncol(tbl) %in%  c(0,15)){
                  if(!quiet) .printf("unexpected colnumber: %d (%s)",ncol(tbl), provider)
                   next;
                 }
           
                   # unless speciesExclusive is FALSE -- in which case
                   # each row needs to have just one of the specified 
                   # species -- eliminate all rows with any other
                   # species involved
                   # columns 10 and 11 are taxonA and taxonB respectively
               
               if(nrow(tbl) > 0){
                 tbl <- cbind(tbl,provider=rep(provider, nrow(tbl)),
                              stringsAsFactors=FALSE)
                 if(speciesExclusive)
                     tbl <- .restrictBySpecies(tbl, species)
                 result <- rbind(result, tbl)
                 } # if nrow > 0
              }# for pair
           } # for provider

       official.column.names <-  c("A", "B", "altA", "altB", "aliasA", "aliasB",
                                   "detectionMethod",
                                   "firstAuthor",
                                   "publicationID",
                                   "taxonA",
                                   "taxonB",
                                   "type",
                                   "sourceDatabases",
                                   "interactionID",
                                   "confidenceScore",
                                   "provider")

       if(nrow(result) > 0){
           if(ncol(result) != 16){
              if(!quiet)
                .printf("%d cols for %s, %s-%s", ncol(result), provider, a, b)
              }
            colnames(result) <- official.column.names
            # convert "-" to NA
            scores <- result$confidenceScore
            scores[scores=="-"] <- NA_character_
            result$confidenceScore <- scores
            } # if nrow(result)

             # create an empty data.frame with the proper columns
       if(nrow(result) == 0)
           result <- as.data.frame(sapply(official.column.names,
                                          function(cn) character(0)),
                                   stringsAsFactors=FALSE)
       result
       }) # interactions

#-------------------------------------------------------------------------------
# query language described here:
#   http://code.google.com/p/psicquic/wiki/MiqlReference27
# 
.runQuery <- function(siteSpecificUrl, a, b, species, speciesExclusive,
                      type, detectionMethod, pubID, quiet=TRUE)
{
    a.specified <- !is.na(a)
    b.specified <- !is.na(b)
    idString <- ""

    if(a.specified & !b.specified)
        idString <- sprintf("identifier:%s", a)

    if(!a.specified & b.specified)
        idString <- sprintf("identifier:%s", b)

    if(a.specified & b.specified)
        idString <- sprintf("identifier:(%s AND %s)", a, b)

    speciesString <- ""
    speciesConjunction <- "AND"
    if(!speciesExclusive)
       speciesConjunction <- "OR"

    if(!any(is.na(species))){
        speciesString <- sprintf("species:%s", species)
        tmp <- "species:"
        if(nchar(idString) > 0)
           tmp <- " AND species:"
        if(length(species) == 1)
           tmp <- sprintf("%s%s", tmp, species[1])
        else{
           tmp <- sprintf("%s(%s", tmp, species[1])
           for(aSpecies in species[2:length(species)])
               tmp <- sprintf("%s %s %s", tmp, speciesConjunction,
                              aSpecies)
           tmp <- sprintf("%s)", tmp)
           }
        speciesString <- tmp
        } # !na(species)


    detectionMethodString <- ""
    if(!any(is.na(detectionMethod))){
        if(nchar(paste(idString,
                       speciesString,
                       sep="")) > 0)
           tmp <- " AND detmethod:"
        else
           tmp <- "detmethod:"
        if(length(detectionMethod) == 1)
           tmp <- sprintf("%s%s", tmp, detectionMethod[1])
        else{
           tmp <- sprintf("%s(%s", tmp, detectionMethod[1])
           for(detMeth in detectionMethod[2:length(detectionMethod)])
               tmp <- sprintf("%s OR %s", tmp, detMeth)
           tmp <- sprintf("%s)", tmp)
           }
        detectionMethodString <- tmp
        #.printf("detectionMethodString: %s", detectionMethodString)
        } # detectionMethod
    
    pubIDString <- ""

    if(!any(is.na(pubID))){
        if(nchar(paste(idString,
                       speciesString,
                       detectionMethodString,
                       sep="")) > 0)
           tmp <- " AND pubid:"
        else
           tmp <- "pubid:"
        if(length(pubID) == 1)
           tmp <- sprintf("%s%s", tmp, pubID[1])
        else{
           tmp <- sprintf("%s(%s", tmp, pubID[1])
           for(pmid in pubID[2:length(pubID)])
               tmp <- sprintf("%s OR %s", tmp, pmid)
           tmp <- sprintf("%s)", tmp)
           }
        pubIDString <- tmp
        } # pubID

    typeString <- ""
    if(!any(is.na(type))){
        if(nchar(paste(idString,
                       speciesString,
                       detectionMethodString,
                       typeString,
                       sep="")) > 0)
           tmp <- " AND type:"
        else
           tmp <- "type:"
        if(length(type) == 1)
           tmp <- sprintf("%s%s", tmp, type[1])
        else{
           tmp <- sprintf("%s(%s", tmp, type[1])
           for(pmid in type[2:length(type)])
               tmp <- sprintf("%s OR %s", tmp, pmid)
           tmp <- sprintf("%s)", tmp)
           }
        typeString <- tmp
        } # detectionMethod

    
    fixed.site.url <- sub("psicquic$", "current/search/query/",
                          siteSpecificUrl)
    query.url <- sprintf("%s%s%s%s%s%s", fixed.site.url, idString,
                         speciesString, detectionMethodString,
                         pubIDString, typeString)

    if(!quiet) .printf("PSICQUIC:::.runQuery: %s", query.url)
    result <- .retrieveData(query.url, quiet)
    return(result)

} # .runQuery
#-------------------------------------------------------------------------------
setMethod ("rawQuery", "PSICQUIC",

   function(object, provider, rawArgs) {
       if(length(provider) > 1){
           msg <- sprintf("rawQuery permits only one provider, you supplied %d: %s",
                          length(provider), paste(provider, collapse=","))
           stop(msg)
           }
       df <- object@df
       base.url <- df[provider, "url"]
       fixed.site.url <- sub("psicquic$", "current/search/query/", base.url)
       complete.url <- sprintf("%s%s", fixed.site.url, rawArgs)
       result <- .retrieveData(complete.url)
       return(result)
       })

#-------------------------------------------------------------------------------
.retrieveData <- function(query.url, quiet=TRUE)
{
       # replace spaces with url encoding equivalent
  
    query.url <- gsub(" ", "%20", query.url, fixed=TRUE)
    query.url <- gsub("(", "%28", query.url, fixed=TRUE)
    query.url <- gsub(")", "%29", query.url, fixed=TRUE)
    
    if(!quiet) .printf("final query.url: %s", query.url)
    
    tryCatch({
        txt <- getURL(query.url)
        if(nchar(txt) == 0)
           return (data.frame(stringsAsFactors=FALSE))
        ftmp <- tempfile()
        write(txt, file=ftmp)
        result <-read.table(file=ftmp, , sep="\t", header=FALSE,
                            fill=TRUE, quote='"',
                            stringsAsFactors=FALSE)
        if(!quiet)
           .printf("--- %s result: %d %d", query.url,
                   nrow(result), ncol(result))
        return(result)
        }, error=function(err) {
           print(sprintf("%s, %s", str(err), "server not responding"))
           print(sprintf("failed url: %s", query.url))
           return(data.frame(stringsAsFactors=FALSE))
       })

} # .retrieveData
#-------------------------------------------------------------------------------
# psicquic semantics are, at root:  query(A, B)
# but with this class we want to support a broader scheme:
#          query(listOfIds)
# this function translates a list of ids into a list of a/b pairs
# including the degenerate case in which there is no B id

.enumerateSearchPairs <- function(ids)
{
    a <- NA
    b <- NA
    
    if(length(ids) == 1){
       a <- ids[1]
    } else if (length(ids) == 2){
       a <- ids[1]
       b <- ids[2]
    } else if(length(ids) > 2){
       tmp <- outer(ids, ids, paste, sep=",")
       pair.strings <- tmp[upper.tri(tmp)]
       gene.pairs <- strsplit(pair.strings, ",")
       a <- sapply(gene.pairs, "[", 1)
       b <- sapply(gene.pairs, "[", 2)
       }

    list(a=a, b=b)

} # .enumerateSearchPairs
#-------------------------------------------------------------------------------
.restrictBySpecies <- function(tbl, species)
{
    found.in.a <- c()
    found.in.b <- c()

    for(aSpecies in species){
        a.matches <- grep(aSpecies, tbl[,10])
        b.matches <- grep(aSpecies, tbl[,11])
        found.in.a <- c(found.in.a, a.matches)
        found.in.b <- c(found.in.b, b.matches)
        }

    keepers <- intersect(unique(found.in.a), unique(found.in.b))

    tbl[keepers,]

} # .restrictBySpecies
#-------------------------------------------------------------------------------
interactionTypes <- function()
{
    browseURL("http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MI&termId=MI:0190&termName=interaction%20type")

} # interactionTypes
#-------------------------------------------------------------------------------
detectionMethods <- function()
{
    browseURL("http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MI&termId=MI:0001&termName=interaction%20detection%20method")

} # detectionMethods
#-------------------------------------------------------------------------------
speciesIDs <- function()
{
    browseURL("http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html")

} # speciesIDs
#-------------------------------------------------------------------------------
