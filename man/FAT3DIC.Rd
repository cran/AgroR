% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FAT3DIC_function.R
\name{FAT3DIC}
\alias{FAT3DIC}
\title{Analysis: DIC experiments in triple factorial}
\usage{
FAT3DIC(
  f1,
  f2,
  f3,
  response,
  norm = "sw",
  alpha.t = 0.05,
  alpha.f = 0.05,
  quali = c(TRUE, TRUE, TRUE),
  mcomp = "tukey",
  grau = c(NA, NA, NA),
  grau12 = NA,
  grau13 = NA,
  grau23 = NA,
  grau21 = NA,
  grau31 = NA,
  grau32 = NA,
  grau123 = NA,
  grau213 = NA,
  grau312 = NA,
  transf = 1,
  constant = 0,
  names.fat = c("F1", "F2", "F3"),
  ylab = "Response",
  xlab = "",
  xlab.factor = c("F1", "F2", "F3"),
  sup = NA,
  fill = "lightblue",
  theme = theme_classic(),
  angulo = 0,
  family = "sans",
  addmean = TRUE,
  errorbar = TRUE,
  dec = 3,
  geom = "bar",
  textsize = 12,
  labelsize = 4,
  point = "mean_sd",
  angle.label = 0
)
}
\arguments{
\item{f1}{Numeric or complex vector with factor 1 levels}

\item{f2}{Numeric or complex vector with factor 2 levels}

\item{f3}{Numeric or complex vector with factor 3 levels}

\item{response}{Numerical vector containing the response of the experiment.}

\item{norm}{Error normality test (\emph{default} is Shapiro-Wilk)}

\item{alpha.t}{Significance level of the multiple comparison test (\emph{default} is 0.05)}

\item{alpha.f}{Level of significance of the F test (\emph{default} is 0.05)}

\item{quali}{Defines whether the factor is quantitative or qualitative (\emph{qualitative})}

\item{mcomp}{Multiple comparison test (Tukey (\emph{default}), LSD, Scott-Knott and Duncan)}

\item{grau}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with three elements.}

\item{grau12}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f1 x f2 and qualitative factor 2 and quantitative factor 1.}

\item{grau13}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f1 x f3 and qualitative factor 3 and quantitative factor 1.}

\item{grau23}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f2 x f3 and qualitative factor 3 and quantitative factor 2.}

\item{grau21}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f2 and qualitative factor 1 and quantitative factor 2.}

\item{grau31}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f3 and qualitative factor 1 and quantitative factor 3.}

\item{grau32}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f2 x f3 and qualitative factor 2 and quantitative factor 3.}

\item{grau123}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 1, in the case of interaction f1 x f2 x f3 and quantitative factor 1.}

\item{grau213}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 2, in the case of interaction f1 x f2 x f3 and quantitative factor 2.}

\item{grau312}{Polynomial degree in case of quantitative factor (\emph{default} is 1). Provide a vector with n levels of factor 3, in the case of interaction f1 x f2 x f3 and quantitative factor 3.}

\item{transf}{Applies data transformation (\emph{default} is 1; for log consider 0; `angular` for angular transformation)}

\item{constant}{Add a constant for transformation (enter value)}

\item{names.fat}{Allows labeling the factors 1, 2 and 3.}

\item{ylab}{Variable response name (Accepts the \emph{expression}() function)}

\item{xlab}{treatments name (Accepts the \emph{expression}() function)}

\item{xlab.factor}{Provide a vector with two observations referring to the x-axis name of factors 1, 2 and 3, respectively, when there is an isolated effect of the factors. This argument uses `parse`.}

\item{sup}{Number of units above the standard deviation or average bar on the graph}

\item{fill}{Defines chart color (to generate different colors for different treatments, define fill = "trat")}

\item{theme}{ggplot2 theme (\emph{default} is theme_classic())}

\item{angulo}{x-axis scale text rotation}

\item{family}{Font family}

\item{addmean}{Plot the average value on the graph (\emph{default} is TRUE)}

\item{errorbar}{Plot the standard deviation bar on the graph (In the case of a segment and column graph) - \emph{default} is TRUE}

\item{dec}{Number of cells}

\item{geom}{Graph type (columns or segments)}

\item{textsize}{Font size}

\item{labelsize}{Label Size}

\item{point}{This function defines whether the point must have all points ("all"), mean ("mean"), standard deviation (\emph{default} - "mean_sd") or mean with standard error ("mean_se") if quali= FALSE. For quali=TRUE, `mean_sd` and `mean_se` change which information will be displayed in the error bar.}

\item{angle.label}{label angle}
}
\value{
The analysis of variance table, the Shapiro-Wilk error normality test, the Bartlett homogeneity test of variances, the Durbin-Watson error independence test, multiple comparison test (Tukey, LSD, Scott-Knott or Duncan) or adjustment of regression models up to grade 3 polynomial, in the case of quantitative treatments. The column chart for qualitative treatments is also returned.For significant triple interaction only, no graph is returned.
}
\description{
Analysis of an experiment conducted in a completely randomized design in a triple factorial scheme using analysis of variance of fixed effects.
}
\note{
The order of the chart follows the alphabetical pattern. Please use `scale_x_discrete` from package ggplot2, `limits` argument to reorder x-axis. The bars of the column and segment graphs are standard deviation.

The function does not perform multiple regression in the case of two or more quantitative factors. The bars of the column and segment graphs are standard deviation.

In the final output when transformation (transf argument) is different from 1, the columns resp and respo in the mean test are returned, indicating transformed and non-transformed mean, respectively.
}
\examples{
library(AgroR)
data(enxofre)
with(enxofre, FAT3DIC(f1, f2, f3, resp))
}
\references{
Principles and procedures of statistics a biometrical approach Steel, Torry and Dickey. Third Edition 1997

Multiple comparisons theory and methods. Departament of statistics the Ohio State University. USA, 1996. Jason C. Hsu. Chapman Hall/CRC.

Practical Nonparametrics Statistics. W.J. Conover, 1999

Ramalho M.A.P., Ferreira D.F., Oliveira A.C. 2000. Experimentacao em Genetica e Melhoramento de Plantas. Editora UFLA.

Scott R.J., Knott M. 1974. A cluster analysis method for grouping mans in the analysis of variance. Biometrics, 30, 507-512.

Ferreira, E. B., Cavalcanti, P. P., and Nogueira, D. A. (2014). ExpDes: an R package for ANOVA and experimental designs. Applied Mathematics, 5(19), 2952.

Mendiburu, F., and de Mendiburu, M. F. (2019). Package ‘agricolae’. R Package, Version, 1-2.
}
\author{
Gabriel Danilo Shimizu, \email{shimizu@uel.br}

Leandro Simoes Azeredo Goncalves

Rodrigo Yudi Palhaci Marubayashi
}
\keyword{DIC}
\keyword{Factorial}
