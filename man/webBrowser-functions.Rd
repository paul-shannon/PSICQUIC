\name{speciesIDs}

\alias{interactionTypes}
\alias{detectionMethods}
\alias{speciesIDs}

\title{Web Browser Functions}

\description{

PSICQUIC uses standard codes to identify species, detection methods, and
interaction types.   In the absence of easy-to-use data files which we
could store with, and provide from this package, we instead offer methods which
drive your web browser to pages which offer the relevant information.
}

\usage{
interactionTypes()
detectionMethods()
speciesIDs()
}


\value{
  None
}

\details{
interactionTypes will drive your web browser to a EMBL-EBI Ontology
Lookup Service devoted to protein-protein interaction types.


detectionMethods will drive your web browser to a EMBL-EBI Ontology
Lookup Service devoted to protein interaction detection methods.


speciesIDs helps in finding the NCBI taxon ids which are used to identify species (e.g., "9606" for Homo
sapiens). This method drives your web browser to the NCBI Taxonomy
Browser (the NCBI Taxonomy Homepage), where you can use the search
facilities to obtain the taxon id for your organism/s of interest.
}

\author{Paul Shannon}

\examples{
  \dontrun{
    psicquic <- PSICQUIC()
    speciesIDs()
    interactionTypes()
    detectionMethods()
    }
}
\seealso{
  providers,
  interactions,
  rawQuery,
  addGeneInfo,
  IDMapper,
}



\keyword{utilities}
