% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick.tasks.R
\name{quick.multi.anode}
\alias{quick.multi.anode}
\title{Make ANODE table of pooled values}
\usage{
quick.multi.anode(my.mult, task, type = 3)
}
\arguments{
\item{my.mult}{multinomial mira dataframe}

\item{task}{ANODE for ANODE table, coef for z-table, R2 for pR2}

\item{type}{Type for car::Anova (default=3)}
}
\value{
Pooled values for task
}
\description{
Helper for pooled multinomial analyses analyses
}
