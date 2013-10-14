\name{PSICQUIC-class}
\docType{class}

\alias{class:PSICQUIC}
\alias{PSICQUIC-class}
\alias{PSICQUIC}
\alias{show,PSICQUIC-method}

\title{PSICQUIC}

\description{
PSICQUIC is an effort from the HUPO Proteomics Standard Initiative (HUPO-PSI) to
standardise programmatic access to molecular interaction databases.  The
Bioconductor PSICQUIC package provides a traditional R interface layered on top
of the PSICQUIC REST interface. Gene symbols are most commonly used in
queries; interactions are returned in a data.frame, characterized by
interaction type, detection method, and publication references.
Confidence scores are sometimes avaialable. Queries may be constrained
by many of these same attributes, i.e., interaction type, detection
method, species, publication identifier, and source database.

}
\details{
  There are two operational differences between the native PSIQUIC REST
  interface and that offered here via the interactions method:
  \describe{
    \item{species exclusivity:}{The REST interface requires only that one
      participant in an interaction be from the named species.  By
      default, we require that both participants are from the named
      species.  This can be controlled by the \code{speciesExclusive}
      logical argument to the interactions method.}
    \item{number of molecules (or identifiers):}{the REST interface permits only zero, one or two
      gene or protein identifiers per query.  We allow any number, zero or
      more and, when the number is greater than or equal to two, the
      interactions returned are only those in which any two of those identifiers
      participate. If you want all interactions which include any of a
      list of identifiers, and you don't care to control for their
      partners, you can accomplish this by issuing successive
      single-identifer interaction queries.}
  }
}



\section{Constructor}{
  \describe{
    \item{}{
      \code{PSICQUIC}: contacts the central PSICQUIC web server,
      discovers currently functioning servers, returns an object
      used in the methods below.}
   }
  }

\section{Methods}{
  \describe{
    \item{}{
      \code{providers(x)}: lists the short names of the data providers
     }
     \item{}{
       \code{interactions(x,id,
	 species, speciesExclusive,
	 type, provider,detectionMethod,
         publicationID, quiet)}: retrieves all interactions matching the specified pattern.
       }
     \item{}{
       \code{rawQuery(x, provider, rawArgs)}: query terms in native  PSICQUIC REST style
       }
   \item{}{
      \code{show(x)}: displays current providers and related data
       }
    }
   }


\section{Functions}{
  \describe{
    \item{}{
      \code{detectionMethods()}: your web browser will display the
                                 PSI-MI ontology for detection methods
      }
      \item{}{
	 \code{interactionTypes()}: your web browser will display the
                                    PSI-MI ontology for molecular
                                    interaction types
        }
      }
    }
   

\author{Paul Shannon}

\examples{
   psicquic <- PSICQUIC()
     # obtain the list of two dozen (or so) currently live
     # PSICQUIC-compliant data providers
   providers(psicquic)
     # a minimal call: get all interactions with MAP3K3, of all types,
     # from all providers.  a data.frame is returned
   tbl.0 <- interactions(psicquic, "MAP3K3", species="9606")

     # build a contingency table, sort it, and see
     # what kinds of interactions were returned, obtained
     # by what detection methods.
     # "-" is used when the provider does not specify a value.
     # you will see a wide range of specificity, in detection method,
     # interaction type, and number of interactions found.

   xtab <- with(tbl.0, as.data.frame(
                table(type, detectionMethod,  provider)))
   xtab <- subset(xtab, Freq > 0) # [order(xtab$Freq, decreasing=TRUE),]
   xtab <- xtab[order(xtab$Freq, decreasing=TRUE),]

     # what interactors were returned?  the IDMapper class in this
     # package converts many PSICQUIC providers' protein identifiers to
     # entrez geneIDs and HUGO gene symbols, via remote calls to
     # biomaRt:

   idMapper <- IDMapper("9606")
   tbl.0g <- addGeneInfo(idMapper, tbl.0)
   with(tbl.0g, head(unique(c(A.sym, B.sym))))

      # we see that MAP2K5 is the most frequently mentioned interacator:

   xtab.sym <- with(tbl.0g, table(c(A.sym, B.sym)))
   head(sort(xtab.sym, decreasing=TRUE))

      # PSIQUIC uses well-devloped ontologies -- controlled vocabularies --
      # which are currently best viewed in a web browser.
      # we provide two convenience functions which will display these
      # hierarchically defined vocabularies:

   # interactionTypes()
   # detectionMethods()

      # NCBI curates taxonomy codes, such as "9606" for Homo sapiens.
      # you can find these codes by using this method, which will
      # drive your browser to the appropriate NCBI web page.

   # speciesIds()

      # use terms from these vocabularies to retrieve interaction
      # information for these two proteins.  note that both of
      # these terms are mid-level in their respective hierarchies
      # and will likely retrieve more specific nested terms

   tbl.2 <- interactions(psicquic, id=c("MAP3K3", "MAP2K5"),
                         species="9606",
                         type="physical association",
                         detectionMethod="affinity chromatography technology")

      # add gene IDs and symbols

   tbl.2g <- addGeneInfo(idMapper, tbl.2)

      # how many publications lie behind these interactions?
   tbl.2g[, c("A.sym", "B.sym", "detectionMethod", "firstAuthor")]

      # the package also provides a convenience method for submitting
      # queries in native MIQL (Molecular Interaction Query Language).
      # the language is defined here:
      # http://code.google.com/p/psicquic/wiki/MiqlReference27

   tbl.3 <- rawQuery(psicquic, "BioGrid", "identifier:ALK AND species:9606")
      # what publications?
   table(tbl.3$V8)
}

\seealso{
  providers,
  interactions,
  rawQuery,
  addGeneInfo,
  IDMapper,
  interactionTypes,
  detectionMethods,
  speciesIDs
}



\keyword{methods}
\keyword{classes}
