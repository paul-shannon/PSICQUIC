# test_IDMapper.R
#-------------------------------------------------------------------------------
# move these to paulsTests when development slows
library(PSICQUIC)
library(RUnit)
#-------------------------------------------------------------------------------
if(!exists("mapper"))
    mapper <- IDMapper("9606")

#-------------------------------------------------------------------------------
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
    
} # paulsTests
#-------------------------------------------------------------------------------
sampleIDs <- function()
{
    c("uniprotkb:P51532",
      "entrez gene/locuslink:6597|BIOGRID:112481",
      "string:9606.ENSP00000223500|uniprotkb:Q9Y323",
      "refseq:NP_001167455",
      "ensembl:ENSG00000104765",
      "ensembl:ENSP00000350720"
      )

 } # samples
#-------------------------------------------------------------------------------
test_ctor <- function()
{
    print("--- test_ctor")

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
                  uniprotkb=1L,
                  unrecognized=0L))

} # test_.categorize
#-------------------------------------------------------------------------------
test_.translate.uniprotkb <- function()
{
    print("--- test_.translate.uniprotkb")
    
    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }
   
    x <- PSICQUIC:::.categorize(rawIDs)$uniprotkb
    tbl.x <- PSICQUIC:::.translate.uniprotkb(mart, x)
    checkEquals(tbl.x, data.frame(id="P51532", geneID="6597", symbol="SMARCA4",
                                 raw.id="uniprotkb:P51532", stringsAsFactors=FALSE))

} # test_.translate.uniprotkb
#-------------------------------------------------------------------------------
test_.translate.string <- function()
{
    print("--- test_.translate.string")
    rawIDs <- sampleIDs()
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
    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
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
    rawIDs <- sampleIDs()
    if(!exists("mart")){
        #mapper <- IDMapper("9606", rawIDs)
        mapper <- IDMapper("9606")
        mart <<- mapper@mart
        }
   
    x <- PSICQUIC:::.categorize(rawIDs)$ensemblProt
    tbl.x <- PSICQUIC:::.translate.ensemblProt(mart, x)
    checkEquals(tbl.x, data.frame(id="ENSP00000350720",
                                  geneID="6597",
                                  symbol="SMARCA4",
                                  raw.id="ensembl:ENSP00000350720",
                                  stringsAsFactors=FALSE))


} # test_.translate.ensemblProt
#-------------------------------------------------------------------------------
test_.translate.locuslink <- function()
{
    print("--- test_.translate.locuslink")

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

    if(!exists("mart")){
        mart <<- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
        }
    
    raw.ids <- unique(c(tbl.myc$A, tbl.myc$B))
    tbl.ids <- PSICQUIC:::.translateAll(mart, raw.ids)
    checkEquals(colnames(tbl.ids),  c("id", "geneID", "symbol", "raw.id"))
    checkTrue(nrow(tbl.ids) >= length(raw.ids))

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
                  A.sym="MYC",
                  B.sym="SHMT1",
                  A.geneID="4609",
                  B.geneID="6470")
  
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

    tbl.1 <- as.data.frame(row.1, stringsAsFactors=FALSE)
    tbl.2 <- as.data.frame(row.2, stringsAsFactors=FALSE)
    tbl <- RefNet:::.smartRbind(tbl.1, tbl.2)    
    tbl.3 <- addGeneInfo(mapper, tbl)
    checkIdentical(tbl.1, tbl.3[1,])
       # second row of returned tbl.3 should match incoming tbl.2, up to the 16th column
    checkIdentical(as.list(tbl[2,1:16]), as.list(tbl.2))
       # the last 4 columns should have these identifiers
    checkEquals(as.list(tbl.3[2,17:20]), list(A.sym="MAP1LC3B",
                                              B.sym="SHMT1",
                                              A.geneID="81631",
                                              B.geneID="6470"))

    
} # test_preserveKnownGeneIdentifiers
#-------------------------------------------------------------------------------
