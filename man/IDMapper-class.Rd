\name{IDMapper-class}
\docType{class}

\alias{class:IDMapper}
\alias{IDMapper-class}
\alias{IDMapper}

\alias{addGeneInfo}
\alias{addGeneInfo,IDMapper-method}
\alias{\S4method{addGeneInfo}{IDMapper}}


\title{IDMapper}

\description{
  PSICQUIC-compliant services typically return verbose protein identifiers
  for the interactions they report.   However, we often wish to know the
  entrez geneID or HUGO symbol of the gene associated with the
  interacting protein.  The IDMapper class, using biomaRt, provides that
  mapping.   The \code{addGeneInfo} method adds four columns to the
  data.frame it is passed: gene symbol and geneID for both interacting proteins.
}

\section{Constructor}{
   \describe{
    \item{}{
      \code{IDMapper(species)}: only human ("9606") is
         currently supported.
     }
   }
   }

\section{Methods}{
  \describe{
     \item{}{
       \code{addGeneInfo(x,tbl)}: add geneID and symbol for the
       A and B columns in the data.frame "tbl", which is expected
       to be a result returned by the interactions method.
       }
    }
  }

 
\usage{
IDMapper(species)
\S4method{addGeneInfo}{IDMapper}(object, tbl.interactions)
}

\arguments{
  \item{species}{
    a \code{character} string, as specified by the NCBI.  Only 
      "9606" (Homo sapiens) is currently supported.
  }
  
  \item{object}{
    a \code{PSICQUIC} object.
  }

  \item{tbl.interactions}{
    A \code{data.frame} such as one returned by the \code{interactions}
    method.  Columns named \code{A} and \code{B} are required, and
    are presumed to contain PSICQUIC provider identifers in a few
    formats.  Most, but not all identifiers returned by all PSICQUIC providers
    are recognized.
  }
}

\value{
  A copy of the data.frame passed in, with 4 columns added: A.sym,
  B.sym, A.geneID, B.geneID.

}

\author{Paul Shannon}

\examples{
   mapper <- IDMapper("9606")
   tbl <- data.frame(A="entrez gene/locuslink:238|BIOGRID:106739",
                     B="entrez gene/locuslink:3718|BIOGRID:109921",
                     stringsAsFactors=FALSE)
   tbl.withGeneInfo <- addGeneInfo(mapper, tbl)
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
