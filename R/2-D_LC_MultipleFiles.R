
#establish connection with the necessay directory. These are huge files!


teas.files <- dir(path='/Volumes/JordanDowell(MasonLab)/SAlix_2019/Salix Bark/BARK_mzML',
                  full.names = T)
#change sample group information if there are replicate samples. in this case there are none
teas.group.info <- data.frame(sample_name = sub(basename(teas.files),
                                                pattern = '.mzML',
                                                replacement = '',
                                                fixed = T),
                              sample_group = c(rep('lvl_1', 59)),
                              stringsAsFactors = F)
#read in MS data
teas.data <- readMSData(teas.files, pdata = new('NAnnotatedDataFrame',teas.group.info),mode = 'onDisk')


#establish chromatogram parameters
cwp <- CentWaveParam(snthresh = 100, peakwidth = c(6, 15),
                     prefilter = c(5, 1000))

#find peaks for each chromatogram


xchr.multi <- findChromPeaks(teas.data %>%
                               filterRt(rt = c(4.8 * 60, 60 * 60)),
                             param = cwp)

teas.data@phenoData@data

#process based on sample names this may need to be changed
pdp.multi <- PeakDensityParam(sampleGroups = teas.data$sample_group, minFraction = 0, bw = 3)



xchr.multi <- groupChromPeaks(xchr.multi, param = pdp.multi)



#fill approprate peaks

xchr.multi <- fillChromPeaks(xchr.multi)

system("say I am FUCKING done!")


xsa.multi <- xsAnnotate(as(xchr.multi, 'xcmsSet'))
xsaF.multi <- groupFWHM(xsa.multi, perfwhm = 0.6)
#start annotations
xsaI.multi <- findIsotopes(xsaF.multi)

xsaC.multi <- groupCorr(xsaI.multi)



xsaFA.multi <- findAdducts(xsaC.multi, polarity = 'positive')

system("say I am FUCKING done!")


plotPsSpectrum(xsaFA.multi, pspec = 1:4, maxlabel = 5)



View(xsaFA.multi@groupInfo)
dim(xsaFA.multi@groupInfo)
