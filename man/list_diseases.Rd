% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_diseases.R
\name{list_diseases}
\alias{list_diseases}
\title{List diseases available in the Cell Marker Accordion disease database}
\usage{
list_diseases(
  species = c("Human", "Mouse"),
  tissue = NULL,
  aberrant_celltype = NULL
)
}
\arguments{
\item{species}{Character string or character string vector specifying the
species for which to extract the associate list of available diseases
Currently, either “Human” and/or “Mouse” are supported. Default is
c("Mouse",“Human”).}

\item{tissue}{Character string or character string vector specifying the tissue
for which to extract the associate list of available diseases.
If NULL, information from all tissues are retrieved.}

\item{aberrant_celltype}{Character string or character string vector
specifying the aberrant celltype for which to extract the associate list of
available tissues. If NULL, information from all aberrant cell types are
retrieved.
@return List of diseases available in the Cell Marker Accordion disease
database}
}
\description{
List diseases available in the Cell Marker Accordion disease database
}
