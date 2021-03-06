% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick.helpers.R
\name{quick.SSP.matrix}
\alias{quick.SSP.matrix}
\title{Put together SS from car::Anova package}
\usage{
quick.SSP.matrix(my.summary, marginality)
}
\arguments{
\item{my.model}{Model to be used}
}
\value{
List of lists with 3 pieces: null model, pre models, and full models
}
\description{
Helper function to make list of null model, pre models,
and full models for calculating type II SS (and eventually
type III SS).
}
\keyword{Explore}
