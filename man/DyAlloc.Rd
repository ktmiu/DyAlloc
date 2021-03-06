% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DyAlloc.R
\name{DyAlloc}
\alias{DyAlloc}
\title{Provides the final treatment group allocation for a new subject entering the study.
 Incorporates all the previous functions and a practitioner can use this function solely
 to obtain treatment assignments for the study subjects. The site default is FALSE.}
\usage{
DyAlloc(df, N, covars, weight, newsub, site = FALSE)
}
\arguments{
\item{df}{data frame. The profiles of existing subjects in the study.}

\item{N}{numeric variable. The number of treatment arms}

\item{covars}{a vector. All the names of covariates in the study.}

\item{weight}{a vector. The weights for overall study, with-in stratum, (site),
and factors/covariates, in this particular order}

\item{newsub}{a vector. The factor profile of the new subject entering the study.}

\item{site}{a boolean value. Whether or not to account for site imbalance (default:FALSE)}
}
\value{
A treatment group allocation
}
\description{
Provides the final treatment group allocation for a new subject entering the study.
 Incorporates all the previous functions and a practitioner can use this function solely
 to obtain treatment assignments for the study subjects. The site default is FALSE.
}
