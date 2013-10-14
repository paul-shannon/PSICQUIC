\name{PSICQUIC-package}
\alias{PSICQUIC-package}
\docType{package}
\title{
Protemics Standard Initiative Common QUery InterfaCe
}
\description{

This package provides a simple functional interface to the rich
resources of the \emph{PSICQUIC}, the HUPO Proteomics Standard
Initiative (HUPO-PSI) project to offer programmatic access to molecular
interaction databases.

 \url{http://code.google.com/p/psicquic/}

 There are two classes in the package:
\code{\linkS4class{PSICQUIC}}
provides the query
interface via the \emph{interactions} method, whose results are
returned in a data.frame; and \code{\linkS4class{IDMapper}}
adds gene symbol and entrez geneIDs to the PSICQUIC results data.frame,
for easier interpretation and subsequent use.
}


\author{
Paul Shannon<pshannon@fhcrc.org>
}
\references{
  http://code.google.com/p/psicquic/
  Aranda, Bruno, et al. "PSICQUIC and PSISCORE: accessing and scoring molecular interactions." Nature methods 8.7 (2011): 528-529.
}
\keyword{datasets}

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

