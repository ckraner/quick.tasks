% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick.plots.R
\name{quick.heat}
\alias{quick.heat}
\title{Heat Map of Correlations}
\usage{
quick.heat(myDF, use = "pairwise.complet.obs")
}
\arguments{
\item{myDF}{Dataframe}

\item{use}{use option passed to cor()}
}
\value{
GGplot of heat map
}
\description{
Helper function to quickly get list of column names
based on partial piece.
}
\examples{
quick.heat(myDF)

}
\keyword{Explore}
