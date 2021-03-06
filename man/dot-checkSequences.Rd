% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createTaxedSeq.R
\name{.checkSequences}
\alias{.checkSequences}
\title{.checkSequences}
\usage{
.checkSequences(DNAdata)
}
\arguments{
\item{DNAdata}{taxedSeq object containing a DNAStringSet or only
the latter}
}
\value{
A DNAStringSet
}
\description{
Check for a DNAstringSet potentially within a taxedSeq object
}
\details{
This internal function checks for a DNAstringSet or taxedSeq
object. If hte latter, it etextracts this DNAStringSet
}
\author{
Emanuel Heitlinger
}
