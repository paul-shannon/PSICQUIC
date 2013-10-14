\name{rawQuery}
\alias{rawQuery}
\alias{rawQuery,PSICQUIC-method}
\alias{\S4method{rawQuery}{PSICQUIC}}

\title{rawQuery}

\description{
  For exploratory situations, or by personal preference, one can bypass
  the interactions API (which we expect will be used by most people most
  of the time).  The rawQuery method allows you to craft your own MIQL
  (molecular interaction query language) query but without needing to
  obtain per-provider URLs.   The sixteen columns of data are returned
  in a data.frame without column titles.
  
}
\usage{
\S4method{rawQuery}{PSICQUIC}(object, provider, rawArgs)
}

\arguments{
  \item{object}{
    a \code{PSICQUIC} object.
  }

  \item{provider}{
    a \code{character} object, one of the approximately two dozen
    PSCIQUIC-compliant services.
  }

  \item{rawArgs}{

    A character string which follows MIQL syntax, which is defined here:

      http://code.google.com/p/psicquic/wiki/MiqlReference27
    
    The encoding of special characters (spaces and parentheses in
    particular) is handled by this method, so there is no need to
    replace, for instance, spaces or parens with their url-safe
    hex codes.
  }
  
}

\value{
  A data.frame with 15 columns, and one row for every interaction,
  described as an annotated A/B relationship, with self-describing
  column names.
}

\author{Paul Shannon}

\examples{
  psicquic <- PSICQUIC()
  providers(psicquic)
      # query all providers for all known Myc interactions
  tbl <- rawQuery(psicquic, "BioGrid", "identifier:ALK AND species:9606")

      # what publications?
  table(tbl$V8)

  rawArgs.1 <- "identifier:(ALK AND MAP3K3) AND species:9606"
  tbl.2 <- rawQuery(psicquic, "iRefIndex", rawArgs.1)

  rawArgs.2 <- paste(rawArgs.1, " AND pubid:(15657099 OR 14743216)", sep="")
  tbl.3 <- rawQuery(psicquic, "iRefIndex", rawArgs.2)

}
\seealso{
  providers,
  interactions,
  addGeneInfo,
  IDMapper,
  interactionTypes,
  detectionMethods,
  speciesIDs
}


\keyword{utilities}
