% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quick.tasks.R
\name{quick.contrast}
\alias{quick.contrast}
\title{Contrast Tables in Pixiedust}
\usage{
quick.contrast(my.model, SS.type = 3, adjustment = "bonferroni",
  test.stat = "Wilks", abbrev.length = 15, pix.int = T,
  pix.method = "html", my.factors = my.contrasts, my.type = my.reg.type,
  skip.me = F)
}
\arguments{
\item{my.model}{Model.}

\item{SS.type}{Type of sums of squares to report, either 2 or 3. (default = 3)}

\item{adjustment}{P-value adjustment sent to p.adjust (default = "bonferroni")}

\item{test.stat}{Character string giving the type of test statistic to be used (in MANOVA). (default="Wilks")}

\item{abbrev.length}{Integer telling how long of a label is too long. Longer than this and will be condensed (default=15)}

\item{pix.int}{Should this be viewable or is this for a document/dashboard? (default=T)}

\item{pix.method}{Print method. (default="html")}

\item{my.factors}{If you only want some of the factors, use this. Otherwise, factors are pulled from the regression model}

\item{my.type}{If you have problems, you can specify the regression type. This is pulled from the model}

\item{manova}{Is this a MAN(C)OVA?}
}
\value{
Either pixiedust object or code (in HTML or LaTeX) for table
}
\description{
Beautiful tables using PHIA and PixieDust
for lm, glm, and mancova.
}
\examples{
quick.contrast(my.model)
}
\keyword{Explore}
