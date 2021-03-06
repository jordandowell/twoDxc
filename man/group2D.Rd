% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoDxc.R
\name{group2D}
\alias{group2D}
\title{Group 2D Peaks}
\usage{
group2D(
  object,
  mod.time,
  dead.time = 0,
  rt.tol = 180,
  rt2.tol = 5,
  ppm.tol = 20,
  parallelized = F
)
}
\arguments{
\item{object}{An xsAnnotate object returned by CAMERA that needs 2D peak
grouping}

\item{mod.time}{The modulation time for the comprehensive 2D run}

\item{dead.time}{The time before the injection reaches the detector, defaults
to 0}

\item{rt.tol}{The grouping tolerance in the RT dimension}

\item{ppm.tol}{The grouping tolerance in the m/z dimension}

\item{parallelized}{If true, runs grouping in parallel. Requires the future
package}
}
\value{
An xsAnnotate2D object with the added pspec2D slot
}
\description{
This function takes an xsAnnotate object produced by CAMERA and groups the
pseudospectrum groups according to their m/z's, 1D retention time, and 2D
retention time.
The default chromatographic/MS tolerances are 5s 2D RT tolerance and 20 ppm
m/z tolerance. Your system's mileage may vary.
}
