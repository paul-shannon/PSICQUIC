setClass("IDMapper",
         slots=c(species="character",
                 mart="Mart"))

#-------------------------------------------------------------------------------
IDMapper <- function(species)
{
    if(!species %in% "9606") {
        stop("IDMapper only supports human ID mapping for now")
        }


    self <- new("IDMapper")
    if(species == "9606"){
       self@species <- "9606"
       message("connecting to biomart...")
       self@mart <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
       }

    self

}
#-------------------------------------------------------------------------------
setGeneric("addGeneInfo", signature="object",
               function(object, tbl.interactions)
           standardGeneric("addGeneInfo"))
                            
#-------------------------------------------------------------------------------
.categorize <- function(rawIDs)
{

    x <- rawIDs
    
    string.entries <- grep("string:", x, value=TRUE)
    uniprotkb.entries <- grep("uniprotkb:", x, value=TRUE)
  
        # string entries can be misconstrued as uniprot entries. prevent that
        # here

    if(length(string.entries) > 0)
        uniprotkb.entries <- setdiff(uniprotkb.entries, string.entries)
  
    locuslink.entries <- grep("locuslink:", x, value=TRUE)
    refseq.entries <- grep("refseq:", x, value=TRUE)
    ensembl.gene.entries <- grep("ensembl:ENSG", x, value=TRUE)
    ensembl.prot.entries <- grep("ensembl:ENSP", x, value=TRUE)
    recognized.entries <- c(uniprotkb.entries, locuslink.entries,
                            refseq.entries, ensembl.gene.entries,
                            ensembl.prot.entries, string.entries)
    unrecognized.entries <- setdiff(x, recognized.entries)
  
    list(ensemblGene=ensembl.gene.entries,
         ensemblProt=ensembl.prot.entries,
         locuslink=locuslink.entries,
         refseq=refseq.entries,
         string=string.entries,
         uniprotkb=uniprotkb.entries,
         unrecognized=unrecognized.entries)
         

} # .categorize
#-------------------------------------------------------------------------------
.translate.uniprotkb <- function(mart, entries)
{
    uniprots <- gsub(".*uniprotkb:([A-Z0-9]*).*", "\\1", entries)
    names(uniprots) <- entries
    filter <- "uniprot_sptrembl"
    columns = c (filter, "entrezgene", "hgnc_symbol");
    tbl.uniprot.trembl <- getBM(filters=filter, values=uniprots,
                                attributes=columns, mart=mart)
    colnames(tbl.uniprot.trembl) <- c("id", "geneID", "symbol")
    raw.ids <- names(uniprots)[match(tbl.uniprot.trembl$id,
                                     as.character(uniprots))]
    tbl.uniprot.trembl$raw.id <- raw.ids
    
    filter <- "uniprot_swissprot_accession"
    columns = c (filter, "entrezgene", "hgnc_symbol");
    tbl.uniprot <- getBM(filters=filter, values=uniprots, attributes=columns,
                         mart=mart)
    colnames(tbl.uniprot) <- c("id", "geneID", "symbol")
    raw.ids <- names(uniprots)[match(tbl.uniprot$id, as.character(uniprots))]
    tbl.uniprot$raw.id <- raw.ids

    tbl <- rbind(tbl.uniprot.trembl, tbl.uniprot)
    tbl$geneID <- as.character(tbl$geneID)

    tbl

} # translate.uniprotkb
#-------------------------------------------------------------------------------
.translate.string <- function(mart, entries)
{
    tbl.string <- data.frame(id=character(0), geneID=character(0),
                             symbol=character(0), raw.id=character(0))
    if(length(entries) == 0)
        return(tbl.string)
    
    string.ensps <- gsub(".*\\.(ENSP[0-9]*)\\|.*", "\\1", entries)
    names(string.ensps) <- entries
    tbl.string <- getBM(filters="ensembl_peptide_id", values=string.ensps,
                        attributes=c("ensembl_peptide_id", "entrezgene",
                                     "hgnc_symbol"),
                        mart=mart)
    colnames(tbl.string) <- c("id", "geneID", "symbol")
    raw.ids <- names(string.ensps)[match(tbl.string$id, as.character(string.ensps))]
    tbl.string$raw.id <- raw.ids
    tbl.string$geneID <- as.character(tbl.string$geneID)

    tbl.string

} # .translate.string
#------------------------------------------------------------------------------- 
.translate.ensemblGene <- function(mart, entries)
{
    tbl.ensg <- data.frame(id=character(0), geneID=character(0),
                             symbol=character(0), raw.id=character(0))
    if(length(entries) == 0)
        return(tbl.ensg)
    
    ensembl.genes <- gsub(".*ensembl:([A-Z0-9_]*).*", "\\1", entries)
    names(ensembl.genes) <- entries
    tbl.ensg <- getBM(filters="ensembl_gene_id", values=ensembl.genes,
                     attributes=c("ensembl_gene_id", "entrezgene", "hgnc_symbol"),
                     mart=mart)
    colnames(tbl.ensg) <- c("id", "geneID", "symbol")
    raw.ids  <- names(ensembl.genes)[match(tbl.ensg$id,
                                           as.character(ensembl.genes))]
    tbl.ensg$raw.id <- raw.ids

    tbl.ensg$geneID <- as.character(tbl.ensg$geneID)

    tbl.ensg


} # .translate.ensemblGene
#-------------------------------------------------------------------------------
.translate.ensemblProt <- function(mart, entries)
{
    tbl.ensp <- data.frame(id=character(0), geneID=character(0),
                           symbol=character(0), raw.id=character(0))

    if(length(entries) == 0)
        return(tbl.ensp)

    ensembl.prots <- gsub(".*ensembl:([A-Z0-9_]*).*", "\\1", entries)
    names(ensembl.prots) <- entries
    tbl.ensp <- getBM(filters="ensembl_peptide_id", values=ensembl.prots,
                      attributes=c("ensembl_peptide_id", "entrezgene",
                                   "hgnc_symbol"),
                      mart=mart)
    colnames(tbl.ensp) <- c("id", "geneID", "symbol")
    raw.ids <- names(ensembl.prots)[match(tbl.ensp$id,
                                          as.character(ensembl.prots))]
    tbl.ensp$raw.id <- raw.ids

    tbl.ensp$geneID <- as.character(tbl.ensp$geneID)

    tbl.ensp


} # .translate.ensemblProt
#-------------------------------------------------------------------------------
.translate.locuslink <- function(mart, entries)
{

   tbl.entrezs <- data.frame(id=character(0), geneID=character(0),
                             symbol=character(0), raw.id=character(0))

    if(length(entries) == 0)
        return(tbl.entrezs)

    entrezs <- gsub(".*locuslink:([A-Z0-9]*).*", "\\1", entries)
    names(entrezs) <- entries
    tbl.entrezs <- getBM(filters="entrezgene", values=entrezs,
                         attributes=c("entrezgene", "entrezgene", "hgnc_symbol"),
                         mart=mart)
    colnames(tbl.entrezs) <- c("id", "geneID", "symbol")
    raw.ids <- names(entrezs)[match(tbl.entrezs$id, as.character(entrezs))]
    tbl.entrezs$raw.id <- raw.ids

    tbl.entrezs$geneID <- as.character(tbl.entrezs$geneID)

    tbl.entrezs


} # .translate.locuslink
#-------------------------------------------------------------------------------
.translate.refseq <- function(mart, entries)
{
    tbl.refseq <- data.frame(id=character(0), geneID=character(0),
                             symbol=character(0), raw.id=character(0))
    if(length(entries) == 0)
        return(tbl.refseq)

    refseqs <- gsub(".*refseq:([A-Z0-9_]*).*", "\\1", entries)
    names(refseqs) <- entries
    tbl.refseqs <- getBM(filters="refseq_peptide", values=refseqs,
                          attributes=c("refseq_peptide", "entrezgene",
                                       "hgnc_symbol"),
                          mart=mart)
    colnames(tbl.refseqs) <- c("id", "geneID", "symbol")
    raw.ids <- names(refseqs)[match(tbl.refseqs$id, as.character(refseqs))]
    tbl.refseqs$raw.id <- raw.ids

    tbl.refseqs$geneID <- as.character(tbl.refseqs$geneID)

    tbl.refseqs

} # .translate.refseq
#-------------------------------------------------------------------------------
.translate.geneSymbol <- function(mart, entries)
{
    tbl.geneSymbol <- data.frame(id=character(0), geneID=character(0),
                                 symbol=character(0), raw.id=character(0))
    if(length(entries) == 0)
        return(tbl.refseq)

    refseqs <- gsub(".*refseq:([A-Z0-9_]*).*", "\\1", entries)
    names(refseqs) <- entries
    tbl.refseqs <- getBM(filters="refseq_peptide", values=refseqs,
                          attributes=c("refseq_peptide", "entrezgene",
                                       "hgnc_symbol"),
                          mart=mart)
    colnames(tbl.refseqs) <- c("id", "geneID", "symbol")
    raw.ids <- names(refseqs)[match(tbl.refseqs$id, as.character(refseqs))]
    tbl.refseqs$raw.id <- raw.ids

    tbl.refseqs$geneID <- as.character(tbl.refseqs$geneID)

    tbl.refseqs

} # .translate.geneSymbol
#-------------------------------------------------------------------------------
.translateAll <- function(mart, raw.ids)
{
    categories <- .categorize(raw.ids)
    result <- data.frame()
    for(category in names(categories)){
        ids <- categories[[category]]
        if(length(ids) > 0){
            if(category == "ensemblGene")
                result <- rbind(result, .translate.ensemblGene(mart, ids))
             if(category == "ensemblProt")
                result <- rbind(result, .translate.ensemblProt(mart, ids))
            if(category == "locuslink")
                result <- rbind(result, .translate.locuslink(mart, ids))
            if(category == "refseq")
                result <- rbind(result, .translate.refseq(mart, ids))
            if(category == "string")
                result <- rbind(result, .translate.string(mart, ids))
            if(category == "uniprotkb")
                result <- rbind(result, .translate.uniprotkb(mart, ids))
          } # if length
       } # for category

    result

} # .translateAll
#-------------------------------------------------------------------------------
setMethod ("addGeneInfo", signature=c(object="IDMapper"),

   function(object, tbl.interactions) {

     #stopifnot(
     raw.ids <- unique(c(tbl.interactions$A, tbl.interactions$B))
     
     tbl.xref <- .translateAll(object@mart, raw.ids)
        # create two named lists, for fast lookup of tbl.interactions$A and $B
     syms <- tbl.xref$symbol
     names(syms) <- tbl.xref$raw.id

     geneIDs <- tbl.xref$geneID
     names(geneIDs) <- tbl.xref$raw.id


     max <- nrow(tbl.interactions)
     a.syms <- vector(length=max)
     b.syms <- vector(length=max)
     a.geneIDs <- vector(length=max)
     b.geneIDs <- vector(length=max)

     A <- tbl.interactions$A
     B <- tbl.interactions$B

     A.sym <- as.character(syms[A])
     B.sym <- as.character(syms[B])
     A.geneID <- as.character(geneIDs[A])
     B.geneID <- as.character(geneIDs[B])
     cbind(tbl.interactions, A.sym, B.sym, A.geneID, B.geneID,
           stringsAsFactors=FALSE)

     }) # addGeneInfo
#-------------------------------------------------------------------------------

