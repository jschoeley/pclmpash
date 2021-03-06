% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pclm_pash.r
\name{pclm.general}
\alias{pclm.general}
\title{PCLM De-aggregation of Life-table}
\usage{
pclm.general(x, y, count.type = c("DX", "LX"), out.step = "auto",
  exposures = NULL, control = list())
}
\arguments{
\item{x}{Vector with start of the interval for age/time classes.}

\item{y}{Vector with counts, e.g. \code{ndx}. It must have the same length as
\code{x}.}

\item{count.type}{Type of the data, deaths(\code{"DX"})(default) or
exposures(\code{"LX"}.)}

\item{out.step}{Age interval length in output aggregated life-table. If set
to \code{"auto"} then the parameter is automatically set to the length of
the shortest age/time interval of \code{x}.}

\item{exposures}{Optional exposures to calculate smooth mortality rates. A
vector of the same length as \code{x} and \code{y}. See reference [1] for
further details.}

\item{control}{List with additional parameters. See
\code{\link{pclm.control}}.}
}
\value{
The output is of \code{"pclm"} class with the components:

\item{\code{grouped}}{Life-table based on aggregated PCLM fit and
  defined by \code{out.step}.}

\item{\code{raw}}{Life-table based on original (raw) PCLM fit.}

\item{\code{fit}}{PCLM fit used to construct life-tables.}

\item{\code{m}}{Interval multiple, see
  \code{\link{pclm.interval.multiple}}, \code{\link{pclm.compmat}}.}

\item{\code{x.div}}{Value of \code{x.div}, see
  \code{\link{pclm.control}}.}

\item{\code{out.step}}{Interval length of aggregated life-table, see
  \code{\link{pclm.control}}.}

\item{\code{control}}{Used control parameters, see
  \code{\link{pclm.control}}.}

\item{\code{warn.list}}{List with warnings.}
}
\description{
De-aggregates a life-table using the PCLM method.
}
\details{
The function has four major steps: \enumerate{ \item{Calculate
interval multiple (\code{\link{pclm.interval.multiple}} to remove fractional
parts from \code{x} vector. The removal of fractional parts is necessary to
build composition matrix.} \item{Calculate composition matrix using
\code{\link{pclm.compmat}}.} \item{Fit PCLM model using
\code{\link{pclm.opt}}.} \item{Calculate aggregated (grouped) life-table
using \code{\link{pclm.aggregate}}.} } More details for PCLM algorithm can be
found in reference [1], but see also \code{\link{pclm.general}} and
\code{\link{pclm.compmat}}.
}
\examples{
library(pash)

AU10 <- Inputlx(x = australia_10y$x, lx = australia_10y$lx,
                nax = australia_10y$nax, nx = australia_10y$nx,
                last_open = TRUE)

AU10pclm <- pclm.general(x = AU10$lt$x, y = AU10$lt$ndx*10000)
AU10pclm
}
\references{
\enumerate{ \item{Rizzi S, Gampe J, Eilers PHC. Efficient
estimation of smooth distributions from coarsely grouped data. Am J
Epidemiol. 2015;182:138?47.} \item{Rizzi S, Thinggaard M, Engholm G, et al.
Comparison of non-parametric methods for ungrouping coarsely aggregated data.
BMC Medical Research Methodology. 2016;16:59. doi:10.1186/s12874-016-0157-8.}
}
}
\seealso{
\code{\link{pclm.compmat}},
  \code{\link{pclm.interval.multiple}}, and \code{\link{pclm.nclasses}}.
}
\author{
Maciej J. Danko <\email{danko@demogr.mpg.de}>
  <\email{maciej.danko@gmail.com}>
}
