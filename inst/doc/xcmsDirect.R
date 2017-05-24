### R code from vignette source 'xcmsDirect.Rnw'

###################################################
### code chunk number 1: LoadLib
###################################################
library(xcms)
library(MassSpecWavelet)


###################################################
### code chunk number 2: LoadData
###################################################
library(msdata)
mzdatapath <- system.file("fticr", package = "msdata")
mzdatafiles <- list.files(mzdatapath, recursive = TRUE, full.names = TRUE)
cat("Starting xcmsDirect.Rnw")


###################################################
### code chunk number 3: ProcessData
###################################################
data.mean <- "data.mean"
xs <- xcmsSet(
        method="MSW",
        files=mzdatafiles,
        scales=c(1,4,9),
        nearbyPeak=T,
        verbose.columns = FALSE,
        winSize.noise=500,
        SNR.method="data.mean",
        snthr=10
)


###################################################
### code chunk number 4: CreateExample
###################################################

xs4 <- xcmsSet(
		method = "MSW",
		files = mzdatafiles[1],
		scales = c(1,4, 9),
		nearbyPeak = T,
		verbose.columns = FALSE,
		winSize.noise = 500,
		SNR.method = "data.mean",
		snthr = 10)

masslist <- xs4@peaks[c(1,4,7),"mz"]
xs4@peaks[,"mz"] <- xs4@peaks[,"mz"] + 0.00001*runif(1,0,0.4)*xs4@peaks[,"mz"] + 0.0001


###################################################
### code chunk number 5: xcmsDirect.Rnw:95-103
###################################################
xs4c <- calibrate(xs4,
		calibrants=masslist,
		method="edgeshift",
		mzabs=0.0001,
		mzppm=5,
		neighbours=3,
		plotres=TRUE
		)


###################################################
### code chunk number 6: MzClust
###################################################
xsg <- group(xs, method="mzClust")
xsg


###################################################
### code chunk number 7: ShowGroups
###################################################
groups(xsg)[1:10,]
peaks(xsg)[groupidx(xsg)[[1]]]


###################################################
### code chunk number 8: FillPeaks
###################################################
groupval(xsg)[1,]
xsgf <- fillPeaks(xsg, method="MSW")
groupval(xsgf, "medret", "into")[1:10,]


###################################################
### code chunk number 9: AnalysisVisualize
###################################################
reporttab <- diffreport(xsgf, "ham4", "ham5", "example", eicmax=4,
                        h=480, w=640)
reporttab[1:4,]


