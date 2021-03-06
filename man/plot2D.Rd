% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twoDxc.R
\name{plot2D}
\alias{plot2D}
\title{Plot a 2D Chromatogram}
\usage{
plot2D(
  object,
  file = 1,
  mod.time,
  dead.time = 0,
  ion,
  mz.tol = c("ppm", "abs"),
  ppm.tol = 20,
  abs.tol = 0.5,
  log.scale = F,
  save.output = F,
  filename = "2dplot.png",
  filepath = ".",
  print.output = T,
  mz.digits = 4
)
}
\arguments{
\item{object}{An MSnExp object}

\item{file}{Which file (numeric, not filename) to pull data from. Default is
1}

\item{mod.time}{The modulation time of the 2D run}

\item{dead.time}{The time it takes for the injection to reach the detector.
Default is 0.}

\item{ion}{Select an ion to generate an extracted ion chromatogram (EIC). If
missing, the total ion chromatogram (TIC) is returned.}

\item{mz.tol}{Whether EIC m/z range is calculated using ppm (relative) or
absolute tolerance. Default is ppm}

\item{ppm.tol}{PPM tolerance, used to calculate m/z range for EIC. Default is
20}

\item{abs.tol}{Absolute tolerance, used to calculate m/z range for EIC.
Default is 0.5 m/z}

\item{log.scale}{Whether intensity coloring should be log scaled. Default is
FALSE}

\item{save.output}{Whether to save the output as a file. Default is FALSE}

\item{filename}{The name for saving the file if desired. Default is
'2dplot.png' plus the file selected plus the ion selected.}

\item{filepath}{Where to save the output file. Default is '.'}

\item{print.output}{Whether to show the plot output. Default is TRUE}

\item{mz.digits}{The number of digits after the decimal to show for m/z's
in the plot title. Default is 4}
}
\value{
A 2-dimensional plot
}
\description{
This function takes an MSnExp object and plots a 2D chromatogram using the
specified modulation time and optional ion extraction. The function relies on
ggplot2 for plotting.
}
