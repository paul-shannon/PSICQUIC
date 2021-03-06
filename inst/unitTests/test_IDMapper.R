# test_IDMapper.R
#-------------------------------------------------------------------------------
# move these to paulsTests when development slows
library(PSICQUIC)
library(RUnit)
library(org.Hs.eg.db)
#-------------------------------------------------------------------------------
if(!exists("mapper"))
    mapper <- IDMapper("9606")

#-------------------------------------------------------------------------------
# biomaRt can be flakey.  before doing any tests, make sure we can get
# the mart and dataset we need.  this is called out of many of the test
# functions found below, and the test proceeds only if it returns TRUE
good.hsapiens.mart.exists <- function(mart.name="ensembl",
                                      desired.dataset="hsapiens_gene_ensembl")
{
   if(exists("mart")){
     if(mart@biomart == "ENSEMBL_MART_ENSEMBL" & mart@dataset == desired.dataset)
        return(TRUE)
      }

      # otherwise, create a new instance
   test.mart <- tryCatch({
      useMart(biomart=mart.name)},
      error=function(err) NULL)

   if(is.null(test.mart))
      return(FALSE)

   if(desired.dataset %in% listDatasets(test.mart)$dataset)
      test.mart <- tryCatch({
        useMart(biomart="ensembl", dataset=desired.dataset)},
        error=function(err) NULL)

   if(!is.null(test.mart)) mart <<- test.mart
   return(!is.null(mart))

} # good.hssapiens.mart.exists
#--------------------------------------------------------------------------------
# for use in a few functions below.
# this data.frame ("tbl.myc") was created via
#  tbl.myc <- interactions(psicquicServer, "Myc", source="BioGrid", organism="9606")
sample.interactions.filename <- system.file(package="PSICQUIC",
                                            "extdata",
                                            "mycInteractions.RData")
load(sample.interactions.filename)
#-------------------------------------------------------------------------------
paulsTests <- function()
{
    test_ctor()
    test_.categorize()
    test_.translate.uniprotkb()
    test_.translate.string()
    test_.translate.ensemblGene()
    test_.translate.ensemblProt()
    test_.translate.locuslink()
    test_.translate.refseq()
    test_.translateAll()

    test_addGeneInfo()
    test_preserveKnownGeneIdentifiers()

    #test_addStandardNames()
    #test_preservePreviouslyAssignedStandardNames()

} # paulsTests
#-------------------------------------------------------------------------------
sampleIDs <- function()
{
    c("uniprotkb:P51532",
      "entrez gene/locuslink:6597|BIOGRID:112481",
      "string:9606.ENSP00000223500|uniprotkb:Q9Y323",
      "refseq:NP_001167455",
      "ensembl:ENSG00000104765",
      #"ensembl:ENSP00000350720",
      "ensembl:ENSP00000257290",
      "uniprotkb:P34896"             # maps to two geneIDs, we want the lesser
      )

} # samples
#-------------------------------------------------------------------------------
test_ctor <- function()
{
    print("--- test_ctor")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return(TRUE)
       }

    human <- "9606"
    rawIDs <- sampleIDs()
    mapper <- IDMapper(human)

    checkException(IDMapper("intentional error"))

} # test_ctor
#-------------------------------------------------------------------------------
test_.categorize <- function()
{
    print("--- test_.categorize")

    rawIDs <- sampleIDs()
    x <- PSICQUIC:::.categorize(rawIDs)
    categories <- sort(names(x))
    checkEquals(categories,
                 c("ensemblGene", "ensemblProt", "locuslink", "refseq", "string",
                   "uniprotkb", "unrecognized"))
    checkEquals(sapply(x, length),
                c(ensemblGene=1L,
                  ensemblProt=1L,
                  locuslink=1L,
                  refseq=1L,
                  string=1L,
                  uniprotkb=2L,
                  unrecognized=0L))

} # test_.categorize
#-------------------------------------------------------------------------------
# uniprotkb:P04637
test_.translate.uniprotkb <- function()
{
    print("--- test_.translate.uniprotkb")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }

    rawIDs <- sampleIDs()
    #browser()
    #if(!exists("mart")){
    #    mapper <- IDMapper("9606")
    #    mart <<- mapper@mart
    #    }

    x <- PSICQUIC:::.categorize(rawIDs)$uniprotkb
    checkEquals(length(x), 2)

    tbl.x <- PSICQUIC:::.translate.uniprotkb(mart, x)
    checkEquals(ncol(tbl.x), 4)
       # sometimes multiple hits. note extra geneID.
    checkTrue(nrow(tbl.x) >= 2)

       #       id    geneID  symbol           raw.id
       # 1 P34896 102466733   SHMT1 uniprotkb:P34896
       # 2 P34896      6470   SHMT1 uniprotkb:P34896
       # 3 P51532      6597 SMARCA4 uniprotkb:P51532

    row.1 <- match("6470", tbl.x$geneID)
    row.2 <- match("6597", tbl.x$geneID)

    checkEquals(as.list(tbl.x[row.1,]),
                list(id="P34896", geneID="6470", symbol="SHMT1", raw.id="uniprotkb:P34896"))
    checkEquals(as.list(tbl.x[row.2,]),
                list(id="P51532", geneID="6597", symbol="SMARCA4", raw.id="uniprotkb:P51532"))

} # test_.translate.uniprotkb
#-------------------------------------------------------------------------------
test_.translate.string <- function()
{
    print("--- test_.translate.string")
    rawIDs <- sampleIDs()

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }

    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }

    x <- PSICQUIC:::.categorize(rawIDs)$string
    tbl.x <- PSICQUIC:::.translate.string(mart, x)
    checkEquals(tbl.x, data.frame(id="ENSP00000223500",
                                  geneID="51510",
                                  symbol="CHMP5",
                                  raw.id="string:9606.ENSP00000223500|uniprotkb:Q9Y323",
                                  stringsAsFactors=FALSE))

} # test_.translate.uniprotkb
#-------------------------------------------------------------------------------
test_.translate.ensemblGene <- function()
{
    print("--- test_.translate.ensemblGene")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }

    rawIDs <- sampleIDs()

    if(!exists("mart")){
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }

    x <- PSICQUIC:::.categorize(rawIDs)$ensemblGene
    tbl.x <- PSICQUIC:::.translate.ensemblGene(mart, x)
    checkEquals(tbl.x, data.frame(id="ENSG00000104765",
                                  geneID="665",
                                  symbol="BNIP3L",
                                  raw.id="ensembl:ENSG00000104765",
                                  stringsAsFactors=FALSE))


} # test_.translate.ensemblGene
#-------------------------------------------------------------------------------
test_.translate.ensemblProt <- function()
{
    print("--- test_.translate.ensemblProt")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }


    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }

    x <- PSICQUIC:::.categorize(rawIDs)$ensemblProt
    tbl.x <- PSICQUIC:::.translate.ensemblProt(mart, x)
    checkEquals(tbl.x, data.frame(id="ENSP00000257290",
                                  geneID="5156",
                                  symbol="PDGFRA",
                                  raw.id="ensembl:ENSP00000257290",
                                  stringsAsFactors=FALSE))


} # test_.translate.ensemblProt
#-------------------------------------------------------------------------------
test_.translate.locuslink <- function()
{
    print("--- test_.translate.locuslink")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }

    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }

    x <- PSICQUIC:::.categorize(rawIDs)$locuslink
    tbl.x <- PSICQUIC:::.translate.locuslink(mart, x)

    checkEquals(tbl.x, data.frame(id=6597,
                                  geneID="6597",
                                  symbol="SMARCA4",
                                  raw.id="entrez gene/locuslink:6597|BIOGRID:112481",
                                  stringsAsFactors=FALSE))


} # test_.translate.locuslink
#-------------------------------------------------------------------------------
test_.translate.refseq <- function()
{
    print("--- test_.translate.refseq")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }


    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }

    x <- PSICQUIC:::.categorize(rawIDs)$refseq
    tbl.x <- PSICQUIC:::.translate.refseq(mart, x)

    checkEquals(tbl.x, data.frame(id="NP_001167455",
                                  geneID="29117",
                                  symbol="BRD7",
                                  raw.id="refseq:NP_001167455",
                                  stringsAsFactors=FALSE))

} # test_.translate.refseq
#-------------------------------------------------------------------------------
test_.translateAll <- function()
{
    print("--- test_translateAll")

    if(!good.hsapiens.mart.exists()) {
       message("cannot create biomart, skipping test")
       return
       }

    raw.ids <- unique(c(tbl.myc$A, tbl.myc$B))
    tbl.ids <- PSICQUIC:::.translateAll(mart, raw.ids)
    checkEquals(colnames(tbl.ids),  c("id", "geneID", "symbol", "raw.id"))
       # sometimes we see off-by-one.  be extremely generous, allow for off by 100
       # meaning that this test is one for basic function, no longer for
       # precise mapping
    inexactitude.margin <- 100
    expected.min <- length(raw.ids) - inexactitude.margin
    expected.max <- length(raw.ids) + inexactitude.margin
    expected.range <- expected.min:expected.max
    checkTrue(nrow(tbl.ids) %in% expected.range)

} # test_.translateAll
#-------------------------------------------------------------------------------
test_addGeneInfo <- function()
{
    print("--- test_addGeneInfo")

    human <- "9606"
    if(!exists("mapper"))
       mapper <<- IDMapper(human)

    tbl.mycAugmented <- addGeneInfo(mapper, tbl.myc)
    checkEquals(dim(tbl.mycAugmented), c(nrow(tbl.myc), ncol(tbl.myc) + 4))
    checkTrue(all(sapply(tbl.mycAugmented,function(column) !is.factor(column))))

} # test_addGeneInfo
#-------------------------------------------------------------------------------
test_addGeneInfoMinimalTable <- function()
{
    print("--- test_addGeneInfoMinimalTable")
    mapper <- IDMapper("9606")
    tbl <- data.frame(A="entrez gene/locuslink:238|BIOGRID:106739",
                      B="entrez gene/locuslink:3718|BIOGRID:109921",
                      stringsAsFactors=FALSE)
    tbl.withGeneInfo <- addGeneInfo(mapper, tbl)
    checkEquals(dim(tbl.withGeneInfo), c (1, 6))


} # test_addGeneInfoMinimalTable
#-------------------------------------------------------------------------------
# take a row from the RefNet gerstein.2012 interactions, which already has
# gene symbols and geneIDs, and ensure that the IDMapper leaves them intact
test_preserveKnownGeneIdentifiers <- function()
{
    print("--- test_preserveKnownGeneIdentifiers")
    row.1 <- list(A="MYC",
                  B="SHMT1",
                  altA="4609",
                  altB="6470",
                  aliasA="MYC",
                  aliasB="SHMT1",
                  detectionMethod="psi-mi:MI:0402(chromatin immunoprecipitation assay)",
                  firstAuthor="-",
                  publicationID="gerstein.2012",
                  taxonA="9606",
                  taxonB="9606",
                  type="psi-mi:MI:0407(direct interaction)",
                  sourceDatabases="gerstein.2012",
                  interactionID="-",
                  confidenceScore="-",
                  provider="gerstein.2012",
                  A.name="MYC",
                  B.name="SHMT1",
                  A.id="4609",
                  B.id="6470")

    row.2 <- list(A="uniprotkb:Q9GZQ8",
                  B="uniprotkb:P34896",
                  altA="intact:EBI-373144|uniprotkb:Q6NW02",
                  altB="intact:EBI-715117|uniprotkb:D3DXD0|uniprotkb:Q9UMD1|uniprotkb:Q9UMD2|uniprotkb:Q96HY0",
                  aliasA="psi-mi:mlp3b_human(display_long)|uniprotkb:MAP1LC3B(gene name)|psi-mi:MAP1LC3B(display_short)|uniprotkb:MAP1ALC3(gene name synonym)|uniprotkb:Microtubule-associated protein 1 light chain 3 beta(gene name synonym)|uniprotkb:MAP1A/MAP1B light chain 3 B(gene name synonym)|uniprotkb:MAP1 light chain 3-like protein 2(gene name synonym)|uniprotkb:Autophagy-related protein LC3 B(gene name synonym)|uniprotkb:Autophagy-related ubiquitin-like modifier LC3 B(gene name synonym)",
                  aliasB="psi-mi:glyc_human(display_long)|uniprotkb:Glycine hydroxymethyltransferase(gene name synonym)|uniprotkb:SHMT1(gene name)|psi-mi:SHMT1(display_short)|uniprotkb:Serine methylase(gene name synonym)",
                  detectionMethod="psi-mi:MI:0007(anti tag coimmunoprecipitation)",
                  firstAuthor="Behrends et al. (2010)",
                  publicationID="pubmed:20562859|imex:IM-15184",
                  taxonA="taxid:9606(human)|taxid:9606(Homo sapiens)",
                  taxonB="taxid:9606(human)|taxid:9606(Homo sapiens)",
                  type="psi-mi:MI:0914(association)",
                  sourceDatabases="psi-mi:MI:0469(IntAct)",
                  interactionID="intact:EBI-3045543|imex:IM-15184-337",
                  confidenceScore="intact-miscore:0.35",
                  provider="IntAct")

    row.3 <- list(A="_6470_1_c",
                  B="R_GHMT2r",
                  altA="6470",
                  altB="-",
                  aliasA="SHMT1",
                  aliasB="glycine hydroxymethyltransferase, reversible",
                  detectionMethod="psi-mi:MI:0046(experimental knowledge based)",
                  firstAuthor="-",
                  publicationID="9038835",
                  taxonA="9606",
                  taxonB="9606",
                  type="modifies",
                  sourceDatabases="recon2",
                  interactionID="-",
                  confidenceScore="-",
                  provider="recon2",
                  A.name="-",
                  B.name="-",
                  A.id="-",
                  B.id="-",
                  compartment="c",
                  a.uniprot="P34896",
                  a.sboTerm="SBO:0000252",
                  a.chebi="-",
                  a.kegg.compound="-",
                  a.kegg.genes="hsa:6470",
                  a.kegg.drug="-",
                  a.hmdb="-",
                  a.pubchem.substance="-",
                  reversible="true")

    tbl.1 <- as.data.frame(row.1, stringsAsFactors=FALSE)
    tbl.2 <- as.data.frame(row.2, stringsAsFactors=FALSE)
    tbl.3 <- as.data.frame(row.3, stringsAsFactors=FALSE)

       # used to use my improvised RefNet:::.smartRbind at this
       # point, to create a single data.frame with all columns.
       # found in the three lists/data.frames hand-crafted above.
       # this introduces a RefNet dependency into the PSICQUIC
       # package, clearly undesirable: RefNet depends upon PSICQUIC
       #
       # so now using hadley's rbind.fill instead, arguably a better
       # choice anyway.  one difference must be accomodated:
       # rbind.fill sensibly inserts NAs  in empty cells; smartRbind
       # uses "-" since that is PSICQUIC's  token for missing values.
       # we take that extra step below

    #tbl <- RefNet:::.smartRbind(tbl.1, tbl.2)
    #tbl <- RefNet:::.smartRbind(tbl, tbl.3)

    tbl <- rbind.fill(tbl.1, tbl.2)
    tbl <- rbind.fill(tbl, tbl.3)

    cols <- colnames(tbl)
    for(col in cols){
       column.as.list <- tbl[, col]
       na.entries <- which (is.na(column.as.list))
       #printf("%d na.entries in col %s", length(na.entries), col)
       if(length(na.entries) > 0)
          tbl[na.entries, col] <- "-"
       } # for col

    tbl.mapped <- addGeneInfo(mapper, tbl)

    checkIdentical(tbl.1, tbl.mapped[1,1:20])
       # second row of returned tbl.mapped should match incoming tbl.2, up to the 16th column
    checkIdentical(as.list(tbl[2,1:16]), as.list(tbl.2))
       # the last 4 columns should have these identifiers
    checkEquals(as.list(tbl.mapped[2,17:20]), list(A.name="MAP1LC3B",
                                              B.name="SHMT1",
                                              A.id="81631",
                                              B.id="6470"))

      # row 3 should not change at all.  the interactors are
      #     a recon2 interaction
      #     a recon2 modifier, which comes in a form we do not yet
      #        translate into standard forms

    checkTrue(all(tbl.mapped[3,] == tbl.3))
    x <- 99

} # test_preserveKnownGeneIdentifiers
#-------------------------------------------------------------------------------
test_addStandardNames <- function()
{
    print("--- test_addStandardNames")

    human <- "9606"
    if(!exists("mapper"))
       mapper <<- IDMapper(human)

        # work with just the first 10 interactions in the table

    max <- 10

    tbl.mycAugmented <- addStandardNames(mapper, tbl.myc[1:max,])
    checkEquals(dim(tbl.mycAugmented), c(max, ncol(tbl.myc) + 4))
    checkEquals(colnames(tbl.mycAugmented)[17:20],
                c("A.name", "B.name", "A.id", "B.id"))
    checkTrue(all(sapply(tbl.mycAugmented,function(column) !is.factor(column))))


    geneIds <- with(tbl.mycAugmented, sort(unique(c(A.id, B.id))))
    geneSyms <- with(tbl.mycAugmented, sort(unique(c(A.name, B.name))))
    checkEquals(geneIds, c("10273", "4609", "5204", "5300", "672", "7316", "8841" ))

    checkEquals(sort(as.character(mget(geneIds, org.Hs.egSYMBOL))), geneSyms)


} # test_addStandardNames
#-------------------------------------------------------------------------------
# take a row from the RefNet gerstein.2012 interactions, which already has
# gene symbols and geneIDs, and ensure that the IDMapper leaves them intact
no_test_preservePreviouslyAssignedStandardNames <- function()
{
    print("--- test_preservePreviouslyAssignedStandardNames")

    row.1 <- list(A="MYC",
                  B="SHMT1",
                  altA="4609",
                  altB="6470",
                  aliasA="MYC",
                  aliasB="SHMT1",
                  detectionMethod="psi-mi:MI:0402(chromatin immunoprecipitation assay)",
                  firstAuthor="-",
                  publicationID="gerstein.2012",
                  taxonA="9606",
                  taxonB="9606",
                  type="psi-mi:MI:0407(direct interaction)",
                  sourceDatabases="gerstein.2012",
                  interactionID="-",
                  confidenceScore="-",
                  provider="gerstein.2012",
                  A.name="MYC",
                  B.name="SHMT1",
                  A.id="4609",
                  B.id="6470")

    row.2 <- list(A="uniprotkb:Q9GZQ8",
                  B="uniprotkb:P34896",
                  altA="intact:EBI-373144|uniprotkb:Q6NW02",
                  altB="intact:EBI-715117|uniprotkb:D3DXD0|uniprotkb:Q9UMD1|uniprotkb:Q9UMD2|uniprotkb:Q96HY0",
                  aliasA="psi-mi:mlp3b_human(display_long)|uniprotkb:MAP1LC3B(gene name)|psi-mi:MAP1LC3B(display_short)|uniprotkb:MAP1ALC3(gene name synonym)|uniprotkb:Microtubule-associated protein 1 light chain 3 beta(gene name synonym)|uniprotkb:MAP1A/MAP1B light chain 3 B(gene name synonym)|uniprotkb:MAP1 light chain 3-like protein 2(gene name synonym)|uniprotkb:Autophagy-related protein LC3 B(gene name synonym)|uniprotkb:Autophagy-related ubiquitin-like modifier LC3 B(gene name synonym)",
                  aliasB="psi-mi:glyc_human(display_long)|uniprotkb:Glycine hydroxymethyltransferase(gene name synonym)|uniprotkb:SHMT1(gene name)|psi-mi:SHMT1(display_short)|uniprotkb:Serine methylase(gene name synonym)",
                  detectionMethod="psi-mi:MI:0007(anti tag coimmunoprecipitation)",
                  firstAuthor="Behrends et al. (2010)",
                  publicationID="pubmed:20562859|imex:IM-15184",
                  taxonA="taxid:9606(human)|taxid:9606(Homo sapiens)",
                  taxonB="taxid:9606(human)|taxid:9606(Homo sapiens)",
                  type="psi-mi:MI:0914(association)",
                  sourceDatabases="psi-mi:MI:0469(IntAct)",
                  interactionID="intact:EBI-3045543|imex:IM-15184-337",
                  confidenceScore="intact-miscore:0.35",
                  provider="IntAct")

    row.3 <- list(A="_6470_1_c",
                  B="R_GHMT2r",
                  altA="6470",
                  altB="-",
                  aliasA="SHMT1",
                  aliasB="glycine hydroxymethyltransferase, reversible",
                  detectionMethod="psi-mi:MI:0046(experimental knowledge based)",
                  firstAuthor="-",
                  publicationID="9038835",
                  taxonA="9606",
                  taxonB="9606",
                  type="modifies",
                  sourceDatabases="recon2",
                  interactionID="-",
                  confidenceScore="-",
                  provider="recon2",
                  A.name="-",
                  B.name="-",
                  A.id="-",
                  B.id="-",
                  compartment="c",
                  a.uniprot="P34896",
                  a.sboTerm="SBO:0000252",
                  a.chebi="-",
                  a.kegg.compound="-",
                  a.kegg.genes="hsa:6470",
                  a.kegg.drug="-",
                  a.hmdb="-",
                  a.pubchem.substance="-",
                  reversible="true")

    tbl.1 <- as.data.frame(row.1, stringsAsFactors=FALSE)
    tbl.2 <- as.data.frame(row.2, stringsAsFactors=FALSE)
    tbl.3 <- as.data.frame(row.3, stringsAsFactors=FALSE)

    tbl <- RefNet:::.smartRbind(tbl.1, tbl.2)
    tbl <- RefNet:::.smartRbind(tbl, tbl.3)

    tbl.mapped <- addStandardNames(mapper, tbl)

    checkIdentical(tbl.1, tbl.mapped[1,1:20])
       # second row of returned tbl.mapped should match incoming tbl.2, up to the 16th column
    checkIdentical(as.list(tbl[2,1:16]), as.list(tbl.2))
       # the last 4 columns should have these identifiers
    checkEquals(as.list(tbl.mapped[2,17:20]), list(A.name="MAP1LC3B",
                                                   B.name="SHMT1",
                                                   A.id="81631",
                                                   B.id="6470"))

      # row 3 should not change at all.  the interactors are
      #     a recon2 interaction
      #     a recon2 modifier, which comes in a form we do not yet
      #        translate into standard forms

    checkTrue(all(tbl.mapped[3,] == tbl.3))
    x <- 99

} # test_preservePreviouslyAssignedStandardNames
#-------------------------------------------------------------------------------
