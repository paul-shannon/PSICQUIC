\name{providers}

\alias{providers}
\alias{providers,PSICQUIC-method}
\alias{\S4method{providers}{PSICQUIC}}

\title{providers}

\description{

  The central PSICQUIC server has only one purpose:  to return a list of, 
  and metadata about, the currently "live" PSICQUIC-compliant MI-QL
  websites.  A list of these providers by the PSICQUIC object you
  create, and is returned by this method.
}

\usage{
\S4method{providers}{PSICQUIC}(object)
}

\arguments{
  \item{object}{a \code{PSICQUIC} object.}
}

\value{
A list of character strings, short readable names for the providers,
e.g., "BioGrid", "IntAct".
}
\author{Paul Shannon}

\examples{
  psicquic <- PSICQUIC()
  providers(psicquic)

}
\seealso{
  interactions,
  rawQuery,
  addGeneInfo,
  IDMapper,
  interactionTypes,
  detectionMethods,
  speciesIDs
}


\keyword{utilities}
