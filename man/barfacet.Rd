% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barfacets_function.R
\name{barfacet}
\alias{barfacet}
\title{Graph: Bar graph for one factor with facets}
\usage{
barfacet(
  model,
  facet = NULL,
  theme = theme_bw(),
  geom = "bar",
  fill = "lightblue",
  pointsize = 4.5,
  width.bar = 0.15
)
}
\arguments{
\item{model}{DIC, DBC or DQL object}

\item{facet}{vector with facets}

\item{theme}{ggplot2 theme}

\item{geom}{graph type (columns or segments)}

\item{fill}{fill bars}

\item{pointsize}{Point size}

\item{width.bar}{width of the error bars of a regression graph.}
}
\value{
Returns a bar chart for one factor
}
\description{
This is a function of the bar graph for one factor with facets
}
\examples{
library(AgroR)
data("laranja")
a=with(laranja, DBC(trat, bloco, resp,
     mcomp = "sk",angle=45,
     ylab = "Number of fruits/plants"))
barfacet(a,c("S1","S1","S1","S1","S1",
             "S2","S2","S3","S3"))
}
\author{
Gabriel Danilo Shimizu, \email{shimizu@uel.br}

Leandro Simoes Azeredo Goncalves

Rodrigo Yudi Palhaci Marubayashi
}
