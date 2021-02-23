#packages necessary
BiocManager::install("xcms")
BiocManager::install("CAMERA")


library(dplyr)
library(future.apply)
library(ggplot2)
source('R/twoDxc.R')

mz.range.195 <- calc.mz.Window(195.0882, 20)
# m/z tolerance window for caffeine's accurate mass is:
print(mz.range.195)
#> [1] 195.0843 195.0921
eic.195 <- chromatogram(tea.data, mz = mz.range.195, rt = c(2500, 3000))
plot(eic.195)

#Now peak detection will be performed on the data as if it were 1D with xcms.
cwp <- CentWaveParam(snthresh = 100, peakwidth = c(6, 15),
                     prefilter = c(5, 1000))
tea.peaks <- findChromPeaks(tea.data %>%
                              filterRt(rt = c(4.8 * 60, 60 * 60)),
                            param = cwp)

pdp <- PeakDensityParam(sampleGroups = 1, minFraction = 1, bw = 3)
tea.peaks <- groupChromPeaks(tea.peaks, pdp)

feature.chroms <- featureChromatograms(tea.peaks, features = 1:4)
par(mar = rep(2, 4))
plot(feature.chroms)

#We will now group ions into pseudospectra using the CAMERA package.

tea.xsa <- xsAnnotate(as(tea.peaks, 'xcmsSet'))
# Group by fwhm
tea.xsaF <- groupFWHM(tea.xsa, perfwhm = 0.6)
#> Start grouping after retention time.
#> Created 562 pseudospectra.
# Label isotopes
tea.xsaI <- findIsotopes(tea.xsaF, mzabs = 0.01)
#> Generating peak matrix!
#> Run isotope peak annotation
#>  % finished: 10  20  30  40  50  60  70  80  90  100
#> Found isotopes: 247
# Group by correlation
tea.xsaC <- groupCorr(tea.xsaF)
#> Start grouping after correlation.
#>
#> Calculating peak correlations in 562 Groups...
#>  % finished: 10  20  30  40  50  60  70  80  90  100
#>
#> Calculating graph cross linking in 562 Groups...
#>  % finished: 10  20  30  40  50  60  70  80  90  100
#> New number of ps-groups:  906
#> xsAnnotate has now 906 groups, instead of 562
# Find adducts
tea.xsaFA <- findAdducts(tea.xsaC, polarity = 'positive')
#>
#> Calculating possible adducts in 906 Groups...
#>  % finished: 10  20  30  40  50  60  70  80  90  100
# Plot example spectrum
plotPsSpectrum(tea.xsaFA, pspec = c(1:1), maxlabel = 5)

filter.Ion <- function(pspectra, ion, ppm = 20){
  ion.range <- calc.mz.Window(ion, ppm)
  filtered.pspectra <- as.data.frame(pspectra) %>%
    filter(mz > ion.range[1] & mz < ion.range[2])
  return(filtered.pspectra)
}

print(getpspectra(tea.xsaFA, grp = c(1:10)) %>%
        filter.Ion(195.0882))




tea.xsa2D <- group2D(tea.xsaFA, 60, 30)

print(tea.xsa2D@pspec2D %>%
        filter.Ion(120))


plan(multiprocess)
timer <- proc.time()
xsa2D.multi.p <- group2D(tea.xsaFA, 60, 30, parallelized = T)

print(xsa2D.multi.p@pspec2D %>%
        filter.Ion(120))

system("say I am FUCKING done!")
View(xsa2D.multi.p@)

plot2D(tea.data, file = 1, mod.time = 60, dead.time = 30, ion = 137.5342)


#> Added  6778  spectra
#> Condensing...
#> Grouped 345 2D pseudospectra
