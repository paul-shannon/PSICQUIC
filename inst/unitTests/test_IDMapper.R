# test_IDMapper.R
#-------------------------------------------------------------------------------
# move these to paulsTests when development slows
library(PSICQUIC)
library(RUnit)
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
