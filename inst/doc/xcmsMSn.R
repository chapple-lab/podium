### R code from vignette source 'xcmsMSn.Rnw'

###################################################
### code chunk number 1: LibraryPreload
###################################################
library(xcms)
library(msdata)


###################################################
### code chunk number 2: RawFiles
###################################################
mzdatapath <- system.file("iontrap", package = "msdata")
list.files(mzdatapath, recursive = TRUE)


###################################################
### code chunk number 3: FileIO
###################################################
library(xcms)
mzdatafiles <- list.files(mzdatapath, pattern="extracted.mzData", recursive = TRUE, full.names = TRUE)
xraw <- xcmsRaw(mzdatafiles[1], includeMSn=TRUE)
xraw


###################################################
### code chunk number 4: SpectraSelection
###################################################
peaks <- findPeaks(xraw, method="MS1")


###################################################
### code chunk number 5: SpectraSelection
###################################################
xs <- xcmsSet(mzdatafiles, method="MS1")
xfrag <- xcmsFragments(xs)
xfrag


###################################################
### code chunk number 6: TreePlot
###################################################
plotTree(xfrag,xcmsFragmentPeakID=6)


###################################################
### code chunk number 7: TextTree
###################################################
plotTree(xfrag,xcmsFragmentPeakID=6, textOnly=TRUE)


