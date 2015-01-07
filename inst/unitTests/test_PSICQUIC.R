# test_PSICQUIC.R
#-------------------------------------------------------------------------------
# options(stringsAsFactors=FALSE)
#-------------------------------------------------------------------------------
# move these to paulsTests when development slows
library(PSICQUIC)
library(RUnit)
library(RCurl)
#-------------------------------------------------------------------------------
kras <- "3845"
tp53 <- "7157"
myc  <- "4609"

hepC  <- "11103"
human <- "9606"
rat   <- "10116"
stickleback <- "69293"
#-------------------------------------------------------------------------------
if(!exists("psicquic"))
    psicquic <- PSICQUIC(test=TRUE)
#-------------------------------------------------------------------------------
printf <- function(...) print(noquote(sprintf(...)))
#-------------------------------------------------------------------------------
paulsTests <- function()
{
    test_initialConditions()
    test_.loadRegistry()
    test_ctor()
    test_.retrieveData()
    test_rawQuery()
    test_.enumerateSearchPairs()
    test_.restrictBySpecies()
    test_interactions()
    test_interactionsTwoGenes()
    test_interactionsFourGenes()
    test_handleEmbeddedSingleQuotes()
    test_retrieveByPubmedID()
    test_retrieveByPubmedIdOnly()
    test_retrieveByOmimId()
    test_retrieveByInteractionType()
    test_retrieveByDetectionMethod()
    test_retrieveBySpeciesId()
    test_smallMoleculeWithoutSpeciesDesignation()
    
} # paulsTests
#-------------------------------------------------------------------------------
test_initialConditions <- function()
{
    print("--- test_initialConditions")
        # the PSICQUIC class avoids factors, using character strings always
        # make sure that nothing in the environment will hide any failure
        # to accomplish that
    x <- data.frame(a=letters, b=LETTERS)
    checkTrue(all(sapply(x, is.factor)))
    
} # test_initialConditions
#-------------------------------------------------------------------------------
test_.loadRegistry <- function()
{
    print("--- test_.loadRegistry")
    tbl <- PSICQUIC:::.loadRegistry()

        # do some simple reasonableness tests
    checkEquals(colnames(tbl), c("url"))
    checkEquals(grep("^http://", tbl$url), 1:nrow(tbl))
    checkTrue(nrow(tbl) > 1)   # ridiculously conservative number of rows
    
} # test_.loadRegistry 
#-------------------------------------------------------------------------------
test_ctor <- function()
{
    print("--- test_ctor")
    psicquic <- PSICQUIC(test=TRUE)
    checkTrue(is(providers(psicquic), "character"))
    checkTrue(length(providers(psicquic)) > 0)
    
} # test_ctor
#-------------------------------------------------------------------------------
# this utility method is the ultimate endpoint of all queries for the
# interactions, and rawQuery methods.
# make sure it works when called directly
test_.retrieveData <- function()
{
    print("--- test_.retrieveData")
    psicquic <- PSICQUIC(test=TRUE)
    available.providers <- providers(psicquic)
    source <- "BioGrid"

    if(!source %in% available.providers)
        source <- "BIND"

    if(!source %in% available.providers)
        source <- "mentha"

    if(!source %in% available.providers){
        print("expected psicquic sources not available, skipping test_.retrieveData");
        return(TRUE)
        }
    
    base.url <- PSICQUIC:::providerUrl(psicquic, source)
    fixed.site.url <- sub("psicquic$", "current/search/query/", base.url)
    args <- "identifier:TP53%20AND%20species:9606"
    full.url <- sprintf("%s%s", fixed.site.url, args)
    tbl <- PSICQUIC:::.retrieveData(full.url)
       # no factors!
    checkTrue(all(sapply(tbl, is.character)))

    checkEquals(ncol(tbl), 15)
    checkTrue(nrow(tbl) > 0)     # 62 x 15 on (11 sep 2013)

    suppressMessages({tbl <- PSICQUIC:::.retrieveData("http://intentionally.bogus")})
    checkEquals(dim(tbl), c(0,0))

} # test_.retrieveData
#-------------------------------------------------------------------------------
test_rawQuery <- function()
{
    print("--- test_rawQuery")

    psicquic <- PSICQUIC(test=TRUE)

        # the minimal query: id and species

    available.providers <- providers(psicquic)
    source <- "BioGrid"

    if(!source %in% available.providers)
        source <- "BIND"

    if(!source %in% available.providers)
        source <- "mentha"

    if(!source %in% available.providers){
        print("expected psicquic sources not available, skipping test_.retrieveData");
        return(TRUE)
        }


    rawArgs.0 <- "identifier:TP53 AND species:9606"
    tbl <- rawQuery(psicquic, source, rawArgs.0)
    checkEquals(ncol(tbl), 15)
    checkTrue(nrow(tbl) > 20)
       # no column names added
    checkEquals(colnames(tbl), paste("V", 1:15, sep=""))


       # specify a AND b.  IntAct and iRefIndex both
       # report interactions between these two
    
    if("iRefIndex" %in% available.providers){
        rawArgs.1 <- "identifier:(ALK AND MAP3K3) AND species:9606"
        tbl.2 <- rawQuery(psicquic, "iRefIndex", rawArgs.1)
        checkEquals(nrow(tbl.2), 4)
    
           # two references support interactions between these two at iRefIndex:
           #  15657099: 1 row
           #  14743216: 3 rows
    
        rawArgs <- sprintf("%s AND pubid:(15657099 OR 14743216)", rawArgs.1)
        tbl.3 <- rawQuery(psicquic, "iRefIndex", rawArgs)
        checkEquals(nrow(tbl.3), 4)
    
        rawArgs <- sprintf("%s AND pubid:15657099", rawArgs.1)
        tbl.4 <- rawQuery(psicquic, "iRefIndex", rawArgs)
        checkEquals(nrow(tbl.4), 1)
    
        rawArgs <- sprintf("%s AND pubid:14743216", rawArgs.1)
        tbl.5 <- rawQuery(psicquic, "iRefIndex", rawArgs)
        checkEquals(nrow(tbl.5), 3)
        } # if iRefIndexl

       # BioGrid ALK yields these interaction types
       #   psi-mi:MI:0407(direct interaction) 7
       #   psi-mi:MI:0915(physical association) 55

    if("BioGrid" %in% available.providers){
        rawArgs <- sprintf("%s AND type:physical association", rawArgs.0)
        tbl.6 <- rawQuery(psicquic, "BioGrid", rawArgs)
        checkTrue(nrow(tbl.6) > 50)  # 1331 (aug 2014)  but just 55 (dec 2013)
        checkEquals(ncol(tbl.6), 15)
        
        rawArgs <- sprintf("%s AND type:direct interaction", rawArgs.0)
        tbl.7 <- rawQuery(psicquic, "BioGrid", rawArgs)
        checkTrue(nrow(tbl.7) >= 15) # 15 (dec 2013) 801 (aug 2014)
        checkEquals(ncol(tbl.7), 15)
    
        rawArgs <- sprintf("%s AND type:(direct interaction OR physical association)",
                           rawArgs.0)
        tbl.8 <- rawQuery(psicquic, "BioGrid", rawArgs)
        checkTrue(nrow(tbl.8) >= 62) # 62 (dec 2013)  2184 (aug 2014)
        } # if BioGrid


} # test_rawQuery
#-------------------------------------------------------------------------------
test_.enumerateSearchPairs <- function()
{
    print("--- test_.enumerateSearchPairs")

    checkEquals(PSICQUIC:::.enumerateSearchPairs(c()), list(a=NA, b=NA))
    checkEquals(PSICQUIC:::.enumerateSearchPairs(c("a")), list(a="a", b=NA))
    checkEquals(PSICQUIC:::.enumerateSearchPairs(c("a", "b")), list(a="a", b="b"))

    pairs <- PSICQUIC:::.enumerateSearchPairs(c("a", "b", "c"))
    checkEquals(pairs, list(a=c("a", "a", "b"),
                            b=c("b", "c", "c")))

} # test_.enumerateSearchPairs
#-------------------------------------------------------------------------------
test_interactions <- function()
{
    print("--- test_interactions")
    psicquic <- PSICQUIC(test=FALSE)
    providers <- c("mentha", "BioGrid")

    if(!all(providers %in% providers(psicquic))){
       print("expected psicquic sources not available, test_interactionTwoGenes");
       return(TRUE)
       }
       # make sure that a failed request returns an empty data.frame
       # with all of the expected column names

    tbl.empty <- interactions(psicquic, provider=providers,
                              id="fubarCompletelyUnrecognizableGeneName",
                              species=stickleback)

    checkEquals(dim(tbl.empty), c(0,16))
                
       # ALK: anaplastic lymphoma kinase, entrez geneID 238
    
    tbl <- interactions(psicquic, id="ALK",
                        species="9606", speciesExclusive=TRUE,
                        provider=providers, quiet=TRUE)

       # no factors.  confidence scores are (as of oct 2013) unparsed strings
    checkTrue(all(sapply(tbl, is.character)))
    
    #browser()
    checkTrue(nrow(tbl) > 1)    # 411 on (10 sep 2013), 165 on (6 jan 2015)

       # except for "confidenceScore", which is numeric,
       # all columns are character
  
              
       # cross-tabulation gives a useful summary of the
       # interactions: how determined, of what nature:
    
    xtab <- as.data.frame(with(tbl, table(detectionMethod,type)),
                          stringsAsFactors=FALSE)
    

    xtab <- xtab[xtab$Freq > 0,]
    xtab <- xtab[order(xtab$Freq, decreasing=TRUE),]

       # do a very simple single test:
       # "physical association" is a common interaction type
   checkTrue(length(grep("physical association", xtab$type)) > 1)
    
       # test now for a nonsensical call
   checkException(interactions(psicquic, provider="bogus provider"),
                   silent=TRUE)
   
} # test_interactions
#-------------------------------------------------------------------------------
test_interactionsTwoGenes <- function()
{
    print("--- test_interactionsTwoGenes")
    psicquic <- PSICQUIC()
    providers <- c("mentha", "BioGrid")
    if(!all(providers %in% providers(psicquic))){
       print("expected psicquic sources not available, test_interactionTwoGenes");
       return(TRUE)
       }

       # ALK: anaplastic lymphoma kinase, entrez geneID 238
       # SHC3: entrez geneID 53358
       # JAK3: entrez geneID 3718
       # HSPD1: entrez geneID 3329

    tbl <- interactions(psicquic, species="9606", id=c("ALK", "JAK3"),
                        provider=providers, quiet=TRUE)

    checkTrue(nrow(tbl) > 1)  # (found nrow=5 on 6 jan 2015)
     
} # test_interactionsTwoGenes
#-------------------------------------------------------------------------------
test_interactionsFourGenes <- function()
{
    print("--- test_interactionsFourGenes")
    psicquic <- PSICQUIC(test=FALSE) 
    
    providers <- c("mentha", "BioGrid", "BIND")
    if(!all(providers %in% providers(psicquic))){
       print("expected psicquic sources not available, test_interactionFourGenes");
       return(TRUE)
       }

       # ALK: anaplastic lymphoma kinase, entrez geneID 238
       # SHC3: entrez geneID 53358
       # JAK3: entrez geneID 3718
       # HSPD1: entrez geneID 3329

    tbl <- interactions(psicquic, species="9606", id=c("ALK", "JAK3", "SHC3", "HSPD1"),
                        provider=providers, quiet=TRUE)
    checkTrue(nrow(tbl) > 1)

    mapper <- IDMapper("9606")
    tbl.2 <- addGeneInfo(mapper, tbl)
    checkEquals(sort(unique(c(tbl.2$A.name, tbl.2$B.name))),
                c("ALK", "HSPD1", "JAK3", "SHC3"))
     
} # test_interactionsFourGenes
#-------------------------------------------------------------------------------
test_retrieveByPubmedID <- function()
{
    print("--- test_retrieveByPubmedID")
    psicquic <- PSICQUIC(test=TRUE)

    provider <- "iRefIndex"
    
    if(provider %in% providers(psicquic)){
       genes <- c("ALK", "MAP3K3")
       tbl <- interactions(psicquic,
                           id=genes,
                           species="9606",
                           provider=provider, quiet=TRUE)
      
       checkEquals(dim(tbl), c(4, 16))
   
       tbl.2 <-interactions(psicquic, species="9606",
                            id=genes,
                            provider=provider, publicationID="15657099",
                            quiet=TRUE)
       checkEquals(dim(tbl.2), c(1, 16))
   
       tbl.3 <-interactions(psicquic, id=genes, species="9606",
                            provider=provider, publicationID="14743216",
                            quiet=TRUE)
       checkEquals(dim(tbl.3), c(3, 16))
       
       tbl.4 <-interactions(psicquic, id=genes, species="9606",
                            provider=provider,
                            publicationID=c("14743216", "15657099"))
   
       checkEquals(dim(tbl.4), c(4, 16))
       } # if provider

} # test_retrieveByPubmedID
#-------------------------------------------------------------------------------
test_retrieveByPubmedIdOnly <- function()
{
    print("--- test_retrieveByPubmedIdOnly")

       # first ensure that this works as a rawQuery
    psicquic <- PSICQUIC(test=TRUE)
    provider <- "IntAct";
    if(!provider %in% providers(psicquic)){
       printf("%s not available, skipping test_retrieveByPubmedIdOnly", provider);
       return(TRUE)
       }
      
       # pmid: 20936779
       # Nature Methods, 2010
       # A human MAP kinase interactome, Sourav Bandyopadhyay et al
    tbl <- rawQuery(psicquic, provider, "species:9606 AND pubid:20936779")
    checkEquals(ncol(tbl), 15)
    checkTrue(nrow(tbl) > 500)

       #  same results through the recommended api?
    tbl.2 <- interactions(psicquic, species="9606",
                          provider=provider, publicationID="20936779",
                          quiet=TRUE)
    checkEquals(ncol(tbl.2), 16)
    checkTrue(nrow(tbl.2) > 500)

} # test_retrieveByPubmedIdOnly
#-------------------------------------------------------------------------------
test_retrieveByOmimId <- function()
{
    print("--- test_retrieveByOmimId")

    psicquic <- PSICQUIC(test=TRUE)

    omim.1 <- "00109135"
    omim.2 <- "00137800"
    
       # http://www.ncbi.nlm.nih.gov/omim/?term=00109135
       # describes the AXL receptor tryosine kinase, geneID 558,
       # oddly, this gene appears in only 7/184 of the interactions
       # returned, all of which come from STRING
       # thus, a buggy retrieval, which we support nonetheless

    tbl.1 <- interactions(psicquic, species="9606", publicationID=omim.1)
    checkEquals(length(grep(omim.1, tbl.1$publicationID)), nrow(tbl.1))

    tbl.2 <- interactions(psicquic, species="9606", publicationID=omim.2)
    checkEquals(length(grep(omim.2, tbl.2$publicationID)), nrow(tbl.2))

    if("STRING" %in% providers(psicquic)){
       tbl.3 <- interactions(psicquic, provider="STRING",
                             publicationID=c(omim.1, omim.2),
                             quiet=TRUE)
       # a weak test, appropriately so:  the retrieval of interactions
       # by STRING, from omim ids, is rather hard to fathom

      checkTrue(nrow(tbl.3) > nrow(tbl.1))
      checkTrue(nrow(tbl.3) > nrow(tbl.2))
      } # if STRING is available

} # test_retrieveByPubmedIdOnly
#-------------------------------------------------------------------------------
test_retrieveByInteractionType <- function()
{
    print("--- test_retrieveByInteractionType")
    psicquic <- PSICQUIC(test=TRUE)

    provider <- "BioGrid"
    if(!provider %in% providers(psicquic))
        return(sprintf("%s not available", provider))
    
    tbl <- interactions(psicquic, id="ALK", species="9606", provider="BioGrid")
    checkTrue(nrow(tbl) > 50)
    checkEquals(ncol(tbl), 16)

       # BioGrid ALK yields these interaction types
       #   psi-mi:MI:0407(direct interaction) 7
       #   psi-mi:MI:0915(physical association) 55

    tbl.1 <- interactions(psicquic, id="ALK",
                          species="9606", provider="BioGrid",
                          type="physical association", quiet=TRUE)

    checkTrue(nrow(tbl.1) > 50)
    checkEquals(ncol(tbl.1), 16)

    tbl.2 <- interactions(psicquic, id="ALK",
                          species="9606", provider="BioGrid",
                          type="direct interaction", quiet=TRUE)

    checkTrue(nrow(tbl.2) > 5)
    checkEquals(ncol(tbl.2), 16)

    tbl.3 <- interactions(psicquic, id="ALK",
                          species="9606", provider="BioGrid",
                          type=c("physical association",
                                            "direct interaction"),
                          quiet=TRUE)

    checkTrue(nrow(tbl.3) > 50)
    checkEquals(ncol(tbl.3), 16)

} # test_retrieveByInteractionType
#-------------------------------------------------------------------------------
# ALK interactions learned from pull down studies ("psi-mi:MI:0096(pull down)")
# are reported (as of 19 sep 2013) by
# BioGrid,  InnateDB, IntAct, MINT, Reactome-FIs, STRING, UniProt, and iRefIndex
# the related laboratory method, "psi-mi:MI:0006(anti bait coip)"
# is reported too, but just from iRefIndex.
# here we try a few combinations of provider & detection method

test_retrieveByDetectionMethod <- function()
{
    
    psicquic <- PSICQUIC(test=FALSE)
    providers <- c("BioGrid","InnateDB","IntAct","MINT","Reactome-FIs","STRING",
                   "UniProt","iRefIndex")
    if(!all(providers %in% providers(psicquic))){
       print("expected psicquic sources not available, test_retrieveByDetectionMethod");
       return(TRUE)
       }

    tbl.0 <- interactions(psicquic, id="ALK", species="9606",
                          detectionMethod="pull down")
    checkEquals(unique(tbl.0$detectionMethod), "psi-mi:MI:0096(pull down)")
    pullDown.count <- nrow (tbl.0)

    tbl.1 <- interactions(psicquic, id="ALK", species="9606",
                          detectionMethod="anti bait coip")
        # some surprising results here: both "anti bait coip" and
        # "anti tag coip are returned
    
    antiBaitCoip.count <- nrow(tbl.1)
    methods <- sort(unique(tbl.1$detectionMethod))

    checkTrue(length(grep("MI:0006", methods)) > 0)

    both.methods <- c("pull down", "anti bait coip")
    tbl.2 <- interactions(psicquic, id="ALK", species="9606",
                          detectionMethod=both.methods)

    bothMethods.count <- nrow (tbl.2)
    checkTrue(bothMethods.count >= (pullDown.count + antiBaitCoip.count))
              
} # test_retrieveByDetectionMethod
#-------------------------------------------------------------------------------
# most interactions are for human, mouse, yeast
# some interactions are between an infecting agent and the human host
# explore that capability here
#
#  hepatitus C virus:  taxid:11103
#  speciesIds() drives your browser to the NCBI taxonomy browser
#
test_retrieveBySpeciesId <- function()
{
    print("--- test_retrieveBySpeciesId")
    psicquic <- PSICQUIC(test=TRUE)

    provider <- "BioGrid"
    if(!provider %in% providers(psicquic))
        return(sprintf("%s not available", provider))
    
    tbl <- interactions(psicquic,
                        species=hepC, speciesExclusive=FALSE,
                        provider="BioGrid", quiet=TRUE)

    checkTrue(nrow(tbl) > 100)  # 148 on (1 oct 2013)

    counts <- as.list(with(tbl, table(c(taxonA, taxonB))))

      # hepC ids should be present in every interaction
    checkEquals(counts$`taxid:11103`, nrow(tbl))

      # make sure that there is at least one hepC/human interaction
    checkTrue("taxid:9606" %in% names(counts))


      # as of (01 october 2013) BioGrid has no interactions
      # restricted to hepC: every hepC interactions is with
      # either human or rat proteins.
      # therefore a speciesExclusive species query, limited
      # to hepC, will at present return an empty data.frame
      # but in time this may change, and a non-empty data.frame
      # will be returned, but with only hepC proteins mentioned.
    
    tbl.1 <- interactions(psicquic,
                          species=hepC, speciesExclusive=TRUE,
                          provider="BioGrid", quiet=TRUE)
    if(nrow(tbl.1) > 0)
      checkEquals(unique(c(tbl.1$taxonA, tbl.1$taxonB)), "taxid:11103")


      # now request interactions between hepC AND human from BioGrid
      # the total count should be equal to counts$`taxid:9606`

    tbl.2 <- interactions(psicquic, species=c(hepC, human), provider="BioGrid",
                          quiet=TRUE, speciesExclusive=TRUE)
    checkEquals(nrow(tbl.2), counts$`taxid:9606`)
    

    tbl.3 <- interactions(psicquic, species=c(stickleback, human), provider="BioGrid",
                          quiet=TRUE, speciesExclusive=TRUE)

   
} # test_retrieveBySpeciesId
#-------------------------------------------------------------------------------
test_.restrictBySpecies <- function()
{
    print("--- test_.restrictBySpecies")
    psicquic <- PSICQUIC()
    provider <- "mentha"   # reports human, mouise, fly and cattle interactions for ALK
    if(!provider %in% providers(psicquic)){
        print("did not find %s in providers, skipping test_.restrictBySpecies")
        return(TRUE)
        }

    tbl <- interactions(psicquic, id="ALK", provider=provider)
    if(length(unique(c(tbl$taxonA, tbl$taxonB))) == 1){
        print("did not find multiple species in available providers, skipping test_.restrictBySpecies")
        return(TRUE)
        }

    checkTrue(length(unique(c(tbl$taxonA, tbl$taxonB))) > 1)

       # ensure that there are some non-9606 taxa reported
       # in this result.  for a detailed look:
    
    xtab.taxa <- as.data.frame(table(c(tbl$taxonA, tbl$taxonB)))
    taxa.variants <- xtab.taxa$Var1
    checkTrue(length(grep("9606", taxa.variants)) < length(taxa.variants)) 
       
    tbl.2 <- PSICQUIC:::.restrictBySpecies(tbl, "9606")

        # ensure that every surving row has a the "9606" string
        # in both taxon entries
    
    all.taxa.restricted <- unique(c(tbl.2$taxonA, tbl.2$taxonB))
    checkEquals(length(all.taxa.restricted), 1)

} # test_.restrictBySpecies
#-------------------------------------------------------------------------------
# 403 lines from MINT for TP53
# something awry on line 223?
# In scan(file, what, nmax, sep, dec, quote, skip, nlines, na.strings,  :
#   EOF within quoted string
# dim(tbl) [1] 222  16
test_handleEmbeddedSingleQuotes <- function()
{
    print("--- test_handleEmbeddedSingleQuote")
    psicquic <- PSICQUIC(test=TRUE)
    provider <- "MINT"
    if(!provider %in% providers(psicquic))
       return(sprintf("%s not available", provider))

        
    mint.tp53.url <-  paste("http://www.ebi.ac.uk/Tools/webservices/psicquic/mint",
                            "/webservices/current/search/query",
                            "/identifier:TP53%20AND%20species:9606",
                            sep="")
    raw.text <- getURL(mint.tp53.url)
      # first authors O'Connor and O'Neill bolix things up if
      # single quotes are not disabled as quoting characters
      # make sure that this query is still a good one.
    checkEquals(grep("'", raw.text), 1)
 
    raw.row.count <- length(strsplit(raw.text, "\n")[[1]])  # 445 on (10 oct 2013)
 
      # make sure we get the same number of rows
    
    tbl <- interactions(psicquic, species="9606", provider=provider, id="TP53",
                        quiet=TRUE, speciesExclusive=FALSE)

    checkEquals(nrow(tbl), raw.row.count)

} # test_handleEmbeddedSingleQuotes
#-------------------------------------------------------------------------------
# if species for a or b is unspecified, probably by the designation "-",
# we want to make sure that we accept the accompanying interaction
# one case where this occurs:  ChEMBL, where a drug belongs to
# to no species.  the current solution is to issue the query without
# species restriction.  

test_smallMoleculeWithoutSpeciesDesignation <- function()
{
    print("--- test_smallMoleculeWithoutSpeciesDesignation")
    psicquic <- PSICQUIC(test=FALSE)
    provider <- "ChEMBL"

    if(!provider %in% providers(psicquic))
       return(sprintf("%s not available", provider))
    
    tbl.chembl <- interactions(psicquic, id="imatinib",
                               provider=provider,
                               speciesExclusive=FALSE,
                               quiet=FALSE)
    hit.count <- nrow(tbl.chembl)
    checkTrue(hit.count > 100)   # 141 on (23 dec 2013)
    freq <- with(tbl.chembl, as.list(table(c(taxonA, taxonB))))

       # 3 species:  "-" indicates "not assigned" or "unspecified"
    checkEquals(sort(names(freq)), c("-", "taxid:10090(mouse)", "taxid:9606(human)"))

       # since the query is only for the drug imatinib, every returned row
       # will include it, and thus the "unassigned" species token should
       # be found once in every row -- though perhaps sometimes in
       # taxonA and sometimes in taxonB
    
    checkEquals(hit.count, freq[["-"]])
    
    
} # test_smallMoleculeWithoutSpeciesDesignation
#-------------------------------------------------------------------------------
# multi-species queries are, in general, less interesting than single-species
# queries, with two important exceptions:
# 
#   1) infection
#   2) small-molecule protein interactions
#
# this exploratory function identifies the ins and outs of case 2, using the
# ChEMBL data source
#    ChEMBL or ChEMBLdb is a manually curated chemical database of
#     bioactive molecules with drug-like properties.[1] It is maintained
#     by the European Bioinformatics Institute (EBI), based on the
#     Wellcome Trust Genome Campus, Hinxton, UK.
#
explore_multiSpeciesQuery <- function()
{
    print("--- explore_multiSpeciesQuery")

    provider <- "ChEMBL"
    psicquic <- PSICQUIC(test=TRUE)
    
    if(!provider %in% providers(psicquic))
       return(sprintf("%s not available", provider))

    base.url <- PSICQUIC:::providerUrl(psicquic, provider)
    fixed.site.url <- sub("psicquic$", "current/search/query/", base.url)
    #args <- "identifier:imatinib"   
    args <- sprintf("identifier:(%s AND %s)", "imatinib", "ABL1")
    full.url <- sprintf("%s%s", fixed.site.url, args)
    tbl <- PSICQUIC:::.retrieveData(full.url)
       # no factors!
    colnames(tbl) <- c("A","B","altA","altB","aliasA","aliasB","detectionMethod",
                       "firstAuthor","publicationID","taxonA","taxonB","type",
                       "sourceDatabases","interactionID","confidenceScore")#,"provider")
    coi <- c("A", "B", "aliasA", "aliasB", "taxonA", "taxonB")
    coi <- c("A", "B", "taxonA", "taxonB")

    tbl <- unique(tbl[, coi])
    checkEquals(tbl$taxonA, rep("-", 3))
    checkEquals(sort(tbl$taxonB), c("taxid:10090(mouse)","taxid:9606(human)",
                                    "taxid:9606(human)"))

    tbl.1 <- interactions(psicquic, id=c("imatinib","ABL1"), provider=provider)
    freq <- as.list(table(c(tbl.1$taxonA, tbl.1$taxonB)))
    checkEquals(sum(as.integer(freq)), 2 * nrow(tbl.1))

       # identify a test query, and ensure that it is actually
       # a meaninful test.  the species-neutral query produces
       # at least one human, one mouse, and many "-" (unspecified)
       # for imatinib
    
    checkTrue(length(grep("mouse", names(freq))) > 0)   # just 1
    checkTrue(length(grep("human", names(freq))) > 0)   # 35
    checkTrue(freq[["-"]] > 30)  # (36 23 dec 2013)

       # should get 35 interactions, with the one mouse interaction
       # left out.   currently returns zero...
    
    tbl.human <- interactions(psicquic, id=c("imatinib","ABL1"),
                              species=c("9606", "-"),
                              provider=provider)
    

} # explore_multiSpeciesQuery
#-------------------------------------------------------------------------------
# skipped. GENEMANIA
doNot_test_GeneMANIA <- function()
{
   psicquic <- PSICQUIC(test=TRUE)
   tbl <- interactions(psicquic, species="9606", id="TP53", provider="GeneMANIA")

   # this test fails now (30 oct 2014) due to GeneMANIA problems
   # checkTrue(nrow(tbl) > 0)  

    
} # doNot_test_genemania
#-------------------------------------------------------------------------------
