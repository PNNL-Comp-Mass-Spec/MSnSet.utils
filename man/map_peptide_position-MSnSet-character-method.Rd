% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/map_peptide_position.R
\docType{methods}
\name{map_peptide_position,MSnSet,character-method}
\alias{map_peptide_position,MSnSet,character-method}
\alias{map_peptide_position,MSnSet-character-method}
\title{Map positions of peptides on protein sequences}
\usage{
\S4method{map_peptide_position}{MSnSet,character}(
  object,
  fasta,
  accession_col = "accession",
  peptide_col = "peptide",
  ...
)
}
\arguments{
\item{object}{\code{MSnSet} object with \code{fData}.}

\item{fasta}{character; path to FASTA file.}

\item{accession_col}{character; name of column in \code{fData(object)}
containing accessions (usually proteins).}

\item{peptide_col}{character; name of column in \code{fData(object)}
containing peptides.}

\item{...}{additional arguments}
}
\description{
Map positions of peptides on protein sequences
}
