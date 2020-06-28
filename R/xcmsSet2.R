#had to modify MPI
xcmsSet2 <-
function(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
                    profmethod = "bin", profparam = list(),
                    polarity = NULL, lockMassFreq=FALSE,
                    mslevel=NULL, nSlaves=0, progressCallback=NULL,
                    scanrange=NULL, ...) {

    object <- new("xcmsSet2")

    ## initialise progress information
    xcms.options <- getOption("BioC")$xcms
    xcms.methods <- c(paste("group", xcms.options$group.methods,sep="."), paste("findPeaks", xcms.options$findPeaks.methods,sep="."),
                      paste("retcor", xcms.options$retcor.methods,sep="."), paste("fillPeaks", xcms.options$fillPeaks.methods,sep="."))
    eval(parse(text=paste("object@progressInfo <- list(",paste(xcms.methods,"=0",sep="",collapse=","),")") ))

    if (is.function(progressCallback))
        object@progressCallback <- progressCallback

    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")

    if (is.null(files))
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)

    ## try making paths absolute
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]

    if(lockMassFreq==TRUE){
        ## remove the 02 files if there here
        lockMass.files<-grep("02.CDF", files)
        if(length(lockMass.files) > 0){
            files<-files[-lockMass.files]
        }
    }

    filepaths(object) <- files

    if (length(files) == 0)
        stop("No NetCDF/mzXML/mzData/mzML files were found.\n")

    ## determine experimental design
    fromPaths <- phenoDataFromPaths(files)
    if (is.null(snames)) {
        snames <- rownames(fromPaths)
    } else {
        rownames(fromPaths) <- snames
    }

    pdata <- phenoData
	if (is.null(pdata)) {
        pdata <- sclass
        if (is.null(pdata))
            pdata <- fromPaths
    }
    
    phenoData(object) <- pdata
    if (is.null(phenoData))
        rownames(phenoData(object)) <- snames

    rtlist <- list(raw = vector("list", length(snames)),
                   corrected = vector("list", length(snames)))

    if ("step" %in% names(profparam)) {
        if ("step" %in% names(list(...)) && profparam$step != list(...)$step) {
            stop("different step values defined in profparam and step arguments")
        }
        profstep <- profparam$step
        profparam <- profparam[names(profparam) != "step"]
    } else if ("step" %in% names(list(...))) {
        profstep <- list(...)$step
    } else {
        profstep <- 0.1
    }

    if ("method" %in% names(profparam)) {
        if (profparam$method != profmethod) {
            stop("different method values defined in profparam and profmethod arguments")
        }
        profmethod <- profparam$method
        profparam <- profparam[names(profparam) != "method"]
    }

    profinfo(object) <- c(list(method = profmethod, step = profstep), profparam)

    object@polarity <- as.character(polarity)
    includeMSn=FALSE

    ## implicitely TRUE if selecting MSn
    includeMSn <- !is.null(mslevel) &&  mslevel>1

    ## implicitely TRUE if MS1 parent peak picking
    xcmsSetArgs <- as.list(match.call())
    if (!is.null(xcmsSetArgs$method)) {
        if (xcmsSetArgs$method=="MS1") {
            includeMSn=TRUE
        }
    }

    parmode <- xcmsParallelSetup(nSlaves=nSlaves)
    runParallel <- parmode$runParallel
    parMode <- parmode$parMode    
    snowclust <- parmode$snowclust
    
        params <- list(...);
        params$profmethod <- profmethod;
        params$profparam <- profparam;
        params$includeMSn <- includeMSn;
        params$scanrange <- scanrange;

        params$mslevel <- mslevel; ## Actually, this is 
        params$lockMassFreq <- lockMassFreq;

        ft <- cbind(file=files,id=1:length(files))
        argList <- apply(ft,1,function(x) list(file=x["file"],id=as.numeric(x["id"]),params=params))

        if (parMode == "MPI") {
            res <- xcmsPapply(argList, findPeaksPar2)
            mpi.close.Rslaves()
        } else if (parMode == "SOCK") {
                res <- xcmsClusterApply(cl=snowclust, x=argList, fun=findPeaksPar2, msgfun=msgfun.featureDetection)
                stopCluster(snowclust)
            } else {
              ## serial mode
              res <- lapply(argList, findPeaksPar2)
            }

        peaklist <- lapply(res, function(x) x$peaks)
        rtlist$raw <-  rtlist$corrected <-  lapply(res, function(x) x$scantime)
        if(lockMassFreq){
            object@dataCorrection[1:length(files)]<-1
        }

    lapply(1:length(peaklist), function(i) {
        if (is.null(peaklist[[i]]))
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else  if (nrow(peaklist[[i]]) == 0)
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) == 1)
            warning("Only 1 peak found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) < 10)
            warning("Only ", nrow(peaklist[[i]]), " peaks found in sample",
                    snames[i], call. = FALSE)
    })

    peaks(object) <- do.call(rbind, peaklist)
    object@rt <- rtlist
	object@xcmsRaw <- lapply(res, function(x) x$xRaw) #store xcmsRaw objects

    object
}

setMethod("show", "xcmsSet2", function(object) {

    cat("An \"xcmsSet2\" object with", nrow(object@phenoData), "samples\n\n")

    cat("Time range: ", paste(round(range(object@peaks[,"rt"]), 1), collapse = "-"),
        " seconds (", paste(round(range(object@peaks[,"rt"])/60, 1), collapse = "-"),
        " minutes)\n", sep = "")
    cat("Mass range:", paste(round(range(object@peaks[,"mz"], na.rm = TRUE), 4), collapse = "-"),
        "m/z\n")
    cat("Peaks:", nrow(object@peaks), "(about",
        round(nrow(object@peaks)/nrow(object@phenoData)), "per sample)\n")
    cat("Peak Groups:", nrow(object@groups), "\n")
    cat("Sample classes:", paste(levels(sampclass(object)), collapse = ", "), "\n")
	cat("xcmsRaw Objects:",length(object@xcmsRaw),"\n\n")

    if (length(object@profinfo)) {
        cat("Profile settings: ")
        for (i in seq(along = object@profinfo)) {
            if (i != 1) cat("                  ")
            cat(names(object@profinfo)[i], " = ", object@profinfo[[i]], "\n", sep = "")
        }
        cat("\n")
    }

    memsize <- object.size(object)
    cat("Memory usage:", signif(memsize/2^20, 3), "MB\n")
})

c.xcmsSet2 <- function(...) {

    lcsets <- list(...)
    object <- new("xcmsSet2")

    peaklist <- vector("list", length(lcsets))
    namelist <- vector("list", length(lcsets))
    if (any(duplicated(unlist(namelist)))) {
        stop("Duplicated sample names\n")
    }

    classlist <- vector("list", length(lcsets))
    cdflist <- vector("list", length(lcsets))
	rawList = vector("list", 0)
    rtraw <- vector("list", 0)
    rtcor <- vector("list", 0)
    nsamp <- 0
    for (i in seq(along = lcsets)) {
        peaklist[[i]] <- peaks(lcsets[[i]])
        namelist[[i]] <- sampnames(lcsets[[i]])
        classlist[[i]] <- sampclass(lcsets[[i]])
        classlist[[i]] <- levels(classlist[[i]])[classlist[[i]]]
        cdflist[[i]] <- filepaths(lcsets[[i]])
        rtraw <- c(rtraw, lcsets[[i]]@rt$raw)
        rtcor <- c(rtcor, lcsets[[i]]@rt$corrected)
		rawList = c(rawList,lcsets[[i]]@xcmsRaw)

        sampidx <- seq(along = namelist[[i]]) + nsamp
        peaklist[[i]][,"sample"] <- sampidx[peaklist[[i]][,"sample"]]
        nsamp <- nsamp + length(namelist[[i]])
    }

    peaks(object) <- do.call(rbind, peaklist)
    sampnames(object) <- unlist(namelist)
    classlist <- unlist(classlist)
    sampclass(object) <- factor(classlist, unique(classlist))
    filepaths(object) <- unlist(cdflist)
    profinfo(object) <- profinfo(lcsets[[1]])
	object@xcmsRaw = rawList
    object@rt <- list(raw = rtraw, corrected = rtcor)

    invisible(object)
}

split.xcmsSet2 <- function(x, f, drop = TRUE, ...) {

    if (!is.factor(f))
        f <- factor(f)
    sampidx <- unclass(f)
    peakmat <- peaks(x)
    samples <- sampnames(x)
    classlabel <- sampclass(x)
    cdffiles <- filepaths(x)
    prof <- profinfo(x)
    rtraw <- x@rt$raw
    rtcor <- x@rt$corrected
	xraw = x@xcmsRaw
	

    lcsets <- vector("list", length(levels(f)))
    names(lcsets) <- levels(f)

    for (i in unique(sampidx)) {
        lcsets[[i]] <- new("xcmsSet2")

        samptrans <- numeric(length(f))
        samptrans[sampidx == i] <- rank(which(sampidx == i))
        samp <- samptrans[peakmat[,"sample"]]
        sidx <- which(samp != 0)
        cpeaks <- peakmat[sidx,, drop=FALSE]
        cpeaks[,"sample"] <- samp[sidx]
        peaks(lcsets[[i]]) <- cpeaks

        sampnames(lcsets[[i]]) <- samples[sampidx == i]
        sampclass(lcsets[[i]]) <- classlabel[sampidx == i, drop = TRUE]
        filepaths(lcsets[[i]]) <- cdffiles[sampidx == i]
        profinfo(lcsets[[i]]) <- prof
        lcsets[[i]]@rt$raw <- rtraw[sampidx == i]
        lcsets[[i]]@rt$corrected <- rtcor[sampidx == i]
		lcsets[[i]]@xcmsRaw = xraw[sampidx==i]
		
    }

    if (drop)
        lcsets <- lcsets[seq(along = lcsets) %in% sampidx]

    lcsets
}

#converts an xcmsSet2 object (xs2) to an xcmsSet object
#useful for compatibility with xcmsSet only methods
setGeneric("convert", function(xs2) standardGeneric("convert"))

setMethod("convert", "xcmsSet2", function(xs2)
{
	object <- new("xcmsSet")

	peaks(object) = peaks(xs2)
    object@groups = xs2@groups
	object@groupidx = xs2@groupidx
	object@filled = xs2@filled
	sampnames(object) = sampnames(xs2)
    sampclass(object) = factor(sampclass(xs2), unique(sampclass(xs2)))
    filepaths(object) = filepaths(xs2)
	object@phenoData= phenoDataFromPaths(filepaths(xs2))
    profinfo(object) = profinfo(xs2)
	object@rt$raw = xs2@rt$raw
	object@rt$corrected = xs2@rt$corrected
	object@dataCorrection=xs2@dataCorrection
	object@polarity = xs2@polarity
	object@progressInfo = xs2@progressInfo
	object@progressCallback = xs2@progressCallback
	
    invisible(object)
})




setMethod("group.density", "xcmsSet2", function(object, bw = 30, minfrac = 0.5, minsamp = 1,
                                               mzwid = 0.25, max = 50, sleep = 0) {

    samples <- sampnames(object)
    classlabel <- sampclass(object)
    classnames <- as.character(unique(sampclass(object)))
    classlabel <- as.vector(unclass(classlabel))
    classnum <- table(classlabel)

    peakmat <- peaks(object)
    porder <- order(peakmat[,"mz"])
    peakmat <- peakmat[porder,, drop=FALSE]
    rownames(peakmat) <- NULL
    retrange <- range(peakmat[,"rt"])

    mass <- seq(peakmat[1,"mz"], peakmat[nrow(peakmat),"mz"] + mzwid, by = mzwid/2)
    masspos <- findEqualGreaterM(peakmat[,"mz"], mass)

    groupmat <- matrix(nrow = 512, ncol = 7 + length(classnum))
    groupindex <- vector("list", 512)

    endidx <- 0
    num <- 0
    gcount <- integer(length(classnum))
    for (i in seq(length = length(mass)-2)) {
        if (i %% 500 == 0) {
            cat(round(mass[i]), "")
            flush.console()
        }
        startidx <- masspos[i]
        endidx <- masspos[i+2]-1
        if (endidx - startidx < 0)
            next
        speakmat <- peakmat[startidx:endidx,,drop=FALSE]
        den <- density(speakmat[,"rt"], bw, from = retrange[1]-3*bw, to = retrange[2]+3*bw)
        maxden <- max(den$y)
        deny <- den$y
        gmat <- matrix(nrow = 5, ncol = 2+length(classnum))
        snum <- 0
        while (deny[maxy <- which.max(deny)] > maxden/20 && snum < max) {
            grange <- descendMin(deny, maxy)
            deny[grange[1]:grange[2]] <- 0
            gidx <- which(speakmat[,"rt"] >= den$x[grange[1]] & speakmat[,"rt"] <= den$x[grange[2]])
            gnum <- classlabel[unique(speakmat[gidx,"sample"])]
            for (j in seq(along = gcount))
                gcount[j] <- sum(gnum == j)
            if (! any(gcount >= classnum*minfrac & gcount >= minsamp))
                next
            snum <- snum + 1
            num <- num + 1
### Double the size of the output containers if they're full
            if (num > nrow(groupmat)) {
                groupmat <- rbind(groupmat, matrix(nrow = nrow(groupmat), ncol = ncol(groupmat)))
                groupindex <- c(groupindex, vector("list", length(groupindex)))
            }
            groupmat[num, 1] <- median(speakmat[gidx, "mz"])
            groupmat[num, 2:3] <- range(speakmat[gidx, "mz"])
            groupmat[num, 4] <- median(speakmat[gidx, "rt"])
            groupmat[num, 5:6] <- range(speakmat[gidx, "rt"])
            groupmat[num, 7] <- length(gidx)
            groupmat[num, 7+seq(along = gcount)] <- gcount
            groupindex[[num]] <- sort(porder[(startidx:endidx)[gidx]])
        }
        if (sleep > 0) {
            plot(den, main = paste(round(min(speakmat[,"mz"]), 2), "-", round(max(speakmat[,"mz"]), 2)))
            for (i in seq(along = classnum)) {
                idx <- classlabel[speakmat[,"sample"]] == i
                points(speakmat[idx,"rt"], speakmat[idx,"into"]/max(speakmat[,"into"])*maxden, col = i, pch=20)
            }
            for (i in seq(length = snum))
                abline(v = groupmat[num-snum+i, 5:6], lty = "dashed", col = i)
            Sys.sleep(sleep)
        }
    }
    cat("\n")

    colnames(groupmat) <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                            "npeaks", classnames)

    groupmat <- groupmat[seq(length = num),,drop=FALSE]
    groupindex <- groupindex[seq(length = num)]

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])

    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)

    groups(object) <- groupmat[uindex,,drop=FALSE]
    groupidx(object) <- groupindex[uindex]

    object
})



setMethod("group", "xcmsSet2", function(object, method=getOption("BioC")$xcms$group.method,
                                       ...) {

    method <- match.arg(method, getOption("BioC")$xcms$group.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("group", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})

#had to modify MPI
setMethod("fillPeaks.chrom", "xcmsSet2", function(object, nSlaves=NULL,expand.mz=1,expand.rt=1,min.width.mz=0,min.width.rt=0) { #changed github
 
    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    files <- filepaths(object)
    samp <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    prof <- profinfo(object)
    rtcor <- object@rt$corrected

    ## Remove groups that overlap with more "well-behaved" groups
    numsamp <- rowSums(groupmat[,(match("npeaks", colnames(groupmat))+1):ncol(groupmat),drop=FALSE])
    uorder <- order(-numsamp, groupmat[,"npeaks"])
    uindex <- rectUnique(groupmat[,c("mzmin","mzmax","rtmin","rtmax"),drop=FALSE],
                         uorder)
    groupmat <- groupmat[uindex,]
    groupindex <- groupidx(object)[uindex]
    gvals <- groupval(object)[uindex,]

    peakrange <- matrix(nrow = nrow(gvals), ncol = 4)
    colnames(peakrange) <- c("mzmin","mzmax","rtmin","rtmax")

    mzmin <- peakmat[gvals,"mzmin"]
    dim(mzmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmin"] <- apply(mzmin, 1, median, na.rm = TRUE)
    mzmax <- peakmat[gvals,"mzmax"]
    dim(mzmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"mzmax"] <- apply(mzmax, 1, median, na.rm = TRUE)
    retmin <- peakmat[gvals,"rtmin"]
    dim(retmin) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmin"] <- apply(retmin, 1, median, na.rm = TRUE)
    retmax <- peakmat[gvals,"rtmax"]
    dim(retmax) <- c(nrow(gvals), ncol(gvals))
    peakrange[,"rtmax"] <- apply(retmax, 1, median, na.rm = TRUE)

    lastpeak <- nrow(peakmat)
    lastpeakOrig <- lastpeak

##    peakmat <- rbind(peakmat, matrix(nrow = sum(is.na(gvals)), ncol = ncol(peakmat)))

    cnames <- colnames(object@peaks)

# Making gvals environment so that when it is repeated for each file it only uses the memory one time
gvals_env <- new.env(parent=baseenv())
assign("gvals", gvals, envir = gvals_env)


    ft <- cbind(file=files,id=1:length(files))
    argList <- apply(ft,1,function(x) {
      ## Add only those samples which actually have NA in them
      if (!any(is.na(gvals[,as.numeric(x["id"])]))) {
        ## nothing to do.
        list()
      } else {
        list(file=x["file"],id=as.numeric(x["id"]), object = object,
             params=list(method="chrom",
               gvals=gvals_env, 
               prof=prof,
               dataCorrection=object@dataCorrection,
               polarity=object@polarity,
               rtcor=object@rt$corrected[[as.numeric(x["id"])]],
               peakrange=peakrange,
			   expand.mz=expand.mz, #added Github
                expand.rt=expand.rt, #added Github
				min.width.mz=min.width.mz, #added Github
                min.width.rt=min.width.rt))	#added Github
      }
    })

  nonemptyIdx <- (sapply(argList, length) > 0)

  if (!any(nonemptyIdx)) {
    ## Nothing to do
    return(invisible(object))
  }
    
  argList <- argList[nonemptyIdx]
    
  parmode <- xcmsParallelSetup(nSlaves=nSlaves)
  runParallel <- parmode$runParallel
  parMode <- parmode$parMode    
  snowclust <- parmode$snowclust

  if (parMode == "MPI") {
    newpeakslist <- xcmsPapply(argList, fillPeaksChromPar2)
    mpi.close.Rslaves()
  } else if (parMode == "SOCK") {
    newpeakslist <- xcmsClusterApply(cl=snowclust, x=argList,
                                     fun=fillPeaksChromPar2,
                                     msgfun=msgfunGeneric)    
    stopCluster(snowclust)
  } else {
    ## serial mode
    newpeakslist <- lapply(argList, fillPeaksChromPar2)
  }



  o <- order(sapply(newpeakslist, function(x) x$myID))
  newpeaks <- do.call(rbind, lapply(newpeakslist[o], function(x) x$newpeaks))

  ## Make sure colnames are compatible 
  newpeaks <- newpeaks[, match(cnames, colnames(newpeaks)), drop = FALSE]
  colnames(newpeaks) <- cnames

  peakmat <- rbind(peakmat, newpeaks)

  for (i in seq(along = files)) {
    naidx <- which(is.na(gvals[,i]))
    
    for (j in seq(along = naidx))
      groupindex[[naidx[j]]] <- c(groupindex[[naidx[j]]], lastpeak+j)

    lastpeak <- lastpeak + length(naidx)
  }

    peaks(object) <- peakmat
    object@filled <- seq((lastpeakOrig+1),nrow(peakmat))
    groups(object) <- groupmat
    groupidx(object) <- groupindex

    invisible(object)
})


setMethod("fillPeaks", "xcmsSet2", function(object, method=getOption("BioC")$xcms$fillPeaks.method,...) {
    method <- match.arg(method, getOption("BioC")$xcms$fillPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("fillPeaks", method, sep=".")
    invisible(do.call(method, alist(object, ...)))
})


setMethod("getEIC", "xcmsSet2", function(object, mzrange, rtrange = 200,
                                        groupidx, sampleidx = sampnames(object),
                                        rt = c("corrected", "raw"), mzExpMeth=c("minMax","ppm"), ppm=50 ){

    files <- filepaths(object)
    grp <- groups(object)
    samp <- sampnames(object)
    prof <- profinfo(object)

    rt <- match.arg(rt)

    if (is.numeric(sampleidx))
        sampleidx <- sampnames(object)[sampleidx]
    numsampidx <- match(sampleidx, sampnames(object))
	

    if (!missing(groupidx)) {
        if (is.numeric(groupidx))
            groupidx <- groupnames(object)[unique(as.integer(groupidx))]
        grpidx <- match(groupidx, groupnames(object, template = groupidx))
    }

    if (missing(mzrange)) {
        if (missing(groupidx))
            stop("No m/z range or groups specified")
        if (any(is.na(groupval(object, value = "mz"))))
            stop('Please use fillPeaks() to fill up NA values !')
        if(mzExpMeth == "minMax") #picks mz limits based on smallest mzmin and largest mzmax across all samples 
		{
			mzmin <- -rowMax(-groupval(object, value = "mzmin"))
			mzmax <- rowMax(groupval(object, value = "mzmax"))
			mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
		}
		else if(mzExpMeth == "ppm") #takes average of the mz values from each sample in the group and then does a symmetric ppm expansion around that value
		{
			cat("EIC_ppm= ",ppm,"\n")
			mz = rowMeans(groupval(object, value="mz"))
			mzmin = -(mz*ppm*1e-6)+mz
			mzmax = mz*ppm*1e-6+mz
			mzrange <- matrix(c(mzmin[grpidx], mzmax[grpidx]), ncol = 2)
		}
		
    } else if (all(c("mzmin","mzmax") %in% colnames(mzrange)))
        mzrange <- mzrange[,c("mzmin", "mzmax"),drop=FALSE]
    else if (is.null(dim(mzrange)))
        stop("mzrange must be a matrix")
    colnames(mzrange) <- c("mzmin", "mzmax")

    if (length(rtrange) == 1) {
        if (missing(groupidx))
            rtrange <- matrix(rep(range(object@rt[[rt]][numsampidx]), nrow(mzrange)),
                              ncol = 2, byrow = TRUE)
        else {
            rtrange <- retexp2(grp[grpidx,c("rtmin","rtmax"),drop=FALSE], rtrange)
        }
    } else if (is.null(dim(rtrange)))
        stop("rtrange must be a matrix or single number")
    colnames(rtrange) <- c("rtmin", "rtmax")

	###prevent Rt range from exceeding sample Rt bounds
	##Calculate bounds
	minRt = max(sapply(object@xcmsRaw,function(x) x@scantime[1])) #find largest minimum Rt across all samples
	maxRt = min(sapply(object@xcmsRaw,function(x) max(x@scantime))) #find smallest maximum Rt across all samples
	##Adjust rtRange matrix
	rtrange[ rtrange[ ,"rtmin"]<minRt, "rtmin"] = minRt
	rtrange[ rtrange[ ,"rtmax"]>maxRt, "rtmax"] = maxRt
	
    if (missing(groupidx))
        gnames <- character(0)
    else
        gnames <- groupidx

    eic <- vector("list", length(sampleidx))
    names(eic) <- sampleidx

    for (i in seq(along = numsampidx)) {

        cat(sampleidx[i], "\n")
        flush.console()
        lcraw <- object@xcmsRaw[[numsampidx[i]]]  #Need to make sure order is preserved!!! (ie. this is referencing the sample it should)
		#cat("\nProcessing Sample:",lcraw@filepath[1],"\n")
        if(length(object@dataCorrection) > 1){
            if(object@dataCorrection[numsampidx[i]] == 1)  ## Error in regular xcmsSet,  dataCorrection[i] is referenced instead of dataCorrection[sampidx[i]] 
                lcraw<-stitch(lcraw, AutoLockMass(lcraw))
        }
        if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[numsampidx[i]]]
        if (length(prof) > 2)
            lcraw@profparam <- prof[seq(3, length(prof))]
        currenteic <- getEIC(lcraw, mzrange, rtrange, step = prof$step)
        eic[[i]] <- currenteic@eic[[1]]
        rm(lcraw)
        gc()
    }
    cat("\n")

		
    invisible(new("xcmsEIC", eic = eic, mzrange = mzrange, rtrange = rtrange,
                  rt = rt, groupnames = gnames))
})


retexp2 <- function(peakrange, symExp = 20) 
{
    peakrange[,"rtmin"] <- peakrange[,"rtmin"]-symExp/2
    peakrange[,"rtmax"] <- peakrange[,"rtmax"]+symExp/2

    peakrange
}


setMethod("getSpec", "xcmsSet2", function(object, rawidx, rt=c("raw","corrected"), type=c("multiple","single","average"), ...) {
	## FIXME: unnecessary dependency on profile matrix?
	##generate spectra for each xcmsRaw object specified in rawidx
	#options(warn=-1)
	specHolder = vector("list",length=length(rawidx))
	count = 1
	for(xraw in rawidx)
	{
	
		lcraw = object@xcmsRaw[[xraw]]
		if (rt == "corrected")
            lcraw@scantime <- object@rt$corrected[[xraw]] 
		
		
		sel <- profRange(lcraw, ...)

		scans <- list(length(sel$scanidx))
		uniquemz <- numeric()
		for (i in seq(along = sel$scanidx)) {
			scans[[i]] <- getScan(lcraw, sel$scanidx[i], sel$mzrange)
			uniquemz <- unique(c(uniquemz, scans[[i]][,"mz"]))#!!! round(scans[[i]][,"mz"],digits=2)
		}
		uniquemz <- sort(uniquemz)

		intmat <- matrix(nrow = length(uniquemz), ncol = length(sel$scanidx))
		for (i in seq(along = sel$scanidx)) {
			scan <- getScan(lcraw, sel$scanidx[i], sel$mzrange)
			intmat[,i] <- approx(scan, xout = uniquemz)$y
		}

		pts <- cbind(mz = uniquemz, intensity = rowMeans(intmat))
		#make sure all spectra are the same length
		specHolder[[count]] = pts
		count = count+1
		rm(lcraw)
	}
	
	if(type=="multiple")
	{
		return(specHolder)
	}
	else if (type=="single")
	{
		return(specHolder[[1]])
	}
	else if(type=="average")
	{
		##average spectra to create composite spectra
		mzs = mat.or.vec(0,1)
		intens = mat.or.vec(0,1)
		#calculate values for padding specs
		initPts = sapply(specHolder,function(x) x[1 ,"mz"])
		rng = range(initPts)
		#cat("\n\nInit Pts: ",initPts,"\tRange: ",rng,"\n")
		minidx = which(initPts == rng[1])
		#minidx = minidx[1] #!!!!! for rounding 
		padNum = findRange(specHolder[[minidx]][ ,"mz"],c(rng[2],rng[2]),NAOK=T) #determine maximum amount of padding needed (first num returned)
		#specHolder_out <<-specHolder
		
		padHolder=list()
		#pad front of specs (to help  mz vals line up as best as possible for subsequent averaging)
		for(x in 1:length(specHolder))
		{
			temp1 = specHolder[[x]][,"mz"]
			temp2 = specHolder[[x]][,"intensity"]
			pad = findRange(specHolder[[x]][ ,"mz"],c(rng[2],rng[2]),NAOK=T)
			times = abs(padNum[1] - pad[1]) ###WARNING CLUGE!!!!!!!
			#cat("\n\nX: ",x,"\ttimes: ",times,"\n") 
			temp1 = c(rep(NA,times), temp1)
			temp2 = c(rep(NA,times), temp2)
			padHolder[[x]] = cbind(mz=temp1,intensity=temp2)
		}
		rm(specHolder)
		
		expand = max(sapply(padHolder,function(x) length(x)))/2 #determine max length, divide by 2 to get max col length
		
		#pad end of specs (to prevent recycling)
		for(x in 1:length(padHolder))
		{
			if(length(padHolder[[x]]) < (expand*2))
			{
				temp1 = padHolder[[x]][,"mz"]
				temp2 = padHolder[[x]][,"intensity"]
				length(temp1) = expand
				mzs = cbind(mzs,temp1)
				length(temp2) = expand
				intens = cbind(intens,temp2)
			}
			else
			{
				mzs = cbind(mzs, padHolder[[x]][ ,"mz"])
				intens = cbind(intens, padHolder[[x]][ ,"intensity"])
			}
			
		}
		mzs = rowMeans(mzs)#,na.rm=T)
		intens = rowMeans(intens)#,na.rm=T)
		matx = cbind(mz=mzs,intensity=intens)
		return(matx)
	}
	else
		stop("Invalid MS type requested.  Valid types are: 'multiple', 'average'")
	
})



setMethod("retcor.obiwarp", "xcmsSet2", function(object, plottype = c("none", "deviation"),
                                                profStep=1, center=NULL,
                                                col = NULL, ty = NULL,
                                                response=1, distFunc="cor_opt",
                                                gapInit=NULL, gapExtend=NULL,
                                                factorDiag=2, factorGap=1,
                                                localAlignment=0, initPenalty=0) {

    if (is.null(gapInit)) {
        if (distFunc=="cor") {gapInit=0.3}
        if (distFunc=="cor_opt") {gapInit=0.3}
        if (distFunc=="cov") {gapInit=0.0}
        if (distFunc=="euc") {gapInit=0.9}
        if (distFunc=="prd") {gapInit=0.0}
    }

    if (is.null(gapExtend)) {
        if (distFunc=="cor") {gapExtend=2.4}
        if (distFunc=="cor_opt") {gapExtend=2.4}
        if (distFunc=="cov") {gapExtend= 11.7}
        if (distFunc=="euc") {gapExtend= 1.8}
        if (distFunc=="prd") {gapExtend= 7.8}
    }

    peakmat <- peaks(object)
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    N <- length(samples)
    corpeaks <- peakmat
    plottype <- match.arg(plottype)

    if (length(object@rt) == 2) {
        rtcor <- object@rt$corrected
    } else {
        fnames <- filepaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            xraw <- object@xcmsRaw[[i]]
            rtcor[[i]] <- xraw@scantime
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }

    rtimecor <- vector("list", N)
    rtdevsmo <- vector("list", N)
    plength <- rep(0, N)

    if (missing(center)) {
        for(i in 1:N){
            plength[i] <- length(which(peakmat[,"sample"]==i))
        }
        center <- which.max(plength)
    }

    cat("center sample: ", samples[center], "\nProcessing: ")
    idx <- which(seq(1,N) != center)
    obj1 <- object@xcmsRaw[[center]] #xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0)
	
	## added t automatically find the correct scan range from the xcmsSet object
	if(length(obj1@scantime) != length(object@rt$raw[[center]])){
		##figure out the scan time range
		scantime.start	<-object@rt$raw[[center]][1]
		scantime.end	<-object@rt$raw[[center]][length(object@rt$raw[[center]])]
		
		scanrange.start	<-which.min(abs(obj1@scantime - scantime.start)) 
		scanrange.end	<-which.min(abs(obj1@scantime - scantime.end))
		scanrange<-c(scanrange.start, scanrange.end)
		obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0, scanrange=scanrange) ###May need to change this 
	} else{
		scanrange<-NULL	
	}

    for (si in 1:length(idx)) {
        s <- idx[si]
        cat(samples[s], " ")

        profStepPad(obj1) <- profStep ## (re-)generate profile matrix, since it might have been modified during previous iteration
		if(is.null(scanrange)){
			obj2 <- object@xcmsRaw[[s]] #xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0)	
		} else{
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0, scanrange=scanrange)
		}
        profStepPad(obj2) <- profStep ## generate profile matrix

        mzmin <-  min(obj1@mzrange[1], obj2@mzrange[1])
        mzmax <-  max(obj1@mzrange[2], obj2@mzrange[2])

        mz <- seq(mzmin,mzmax, by=profStep)
        mz <- as.double(mz)
        mzval <- length(mz)

        scantime1 <- obj1@scantime
        scantime2 <- obj2@scantime

        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        rtup1 <- c(1:length(scantime1))
        rtup2 <- c(1:length(scantime2))

        mst1 <- which(diff(scantime1)>5*mstdiff)[1]
        if(!is.na(mst1)) {
            rtup1 <- which(rtup1<=mst1)
            cat("Found gaps: cut scantime-vector at ", scantime1[mst1],"seconds", "\n")
        }

        mst2 <- which(diff(scantime2)>5*mstdiff)[1]
        if(!is.na(mst2)) {
            rtup2 <- which(rtup2<=mst2)
            cat("Found gaps: cut scantime-vector at ", scantime2[mst2],"seconds", "\n")
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]

        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        if(rtmaxdiff>(5*mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            rtup1 <- which(scantime1<=rtmax)
            rtup2 <- which(scantime2<=rtmax)
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)

        if(length(obj1@scantime)>valscantime1) {
            obj1@env$profile <- obj1@env$profile[,-c((valscantime1+1):length(obj1@scantime))]
        }
        if(length(obj2@scantime)>valscantime2) {
            obj2@env$profile <- obj2@env$profile[,-c((valscantime2+1):length(obj2@scantime))]
        }

        if(mzmin < obj1@mzrange[1]) {
            seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep))-1
            x <- matrix(0, seqlen,dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(x, obj1@env$profile)
        }
        if(mzmax > obj1@mzrange[2]){
            seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(obj1@env$profile, x)
        }
        if(mzmin < obj2@mzrange[1]){
            seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(x, obj2@env$profile)
        }
        if(mzmax > obj2@mzrange[2]){
            seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(obj2@env$profile, x)
        }

        intensity1 <- obj1@env$profile
        intensity2 <- obj2@env$profile

        if ((mzval * valscantime1 != length(intensity1)) ||  (mzval * valscantime2 != length(intensity2)))
            stop("Dimensions of profile matrices do not match !\n")

        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              response, distFunc,
                              gapInit, gapExtend,
                              factorDiag, factorGap,
                              localAlignment, initPenalty)

        if(length(obj2@scantime) > valscantime2) {
            object@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            object@rt$corrected[[s]] <- rtimecor[[s]]
        }

        rtdevsmo[[s]] <- round(rtcor[[s]]-object@rt$corrected[[s]],2)

        rm(obj2)
        gc()

        ## updateProgressInfo
        object@progressInfo$retcor.obiwarp <-  si / length(idx)
        xcms:::progressInfoUpdate(object)

    }

    cat("\n")
    rtdevsmo[[center]] <- round(rtcor[[center]] - object@rt$corrected[[center]], 2)

    if (plottype == "deviation") {

        ## Set up the colors and line type
        if (missing(col)) {
            col <- integer(N)
            for (i in 1:max(classlabel))
                col[classlabel == i] <- 1:sum(classlabel == i)
        }
        if (missing(ty)) {
            ty <- integer(N)
            for (i in 1:max(col))
                ty[col == i] <- 1:sum(col == i)
        }
        if (length(palette()) < max(col))
            mypal <- rainbow(max(col), end = 0.85)
        else
            mypal <- palette()[1:max(col)]

        rtrange <- range(do.call("c", rtcor))
        devrange <- range(do.call("c", rtdevsmo))

        layout(matrix(c(1, 2),ncol=2,  byrow=F),widths=c(1,0.3))
        par(mar=c(4,4,2,0))

        plot(0, 0, type="n", xlim = rtrange, ylim = devrange,
             main = "Retention Time Deviation vs. Retention Time",
             xlab = "Retention Time", ylab = "Retention Time Deviation")

        for (i in 1:N) {
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }

        plot.new() ;  par(mar= c(2, 0, 2, 0))
        plot.window(c(0,1), c(0,1))
        legend(0,1.04, basename(samples), col = mypal[col], lty = ty)
    }

    for (i in 1:N) {
        cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
        rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)

})


setMethod("retcor", "xcmsSet2", function(object, method=getOption("BioC")$xcms$retcor.method,
                                        ...) {

    ## Backward compatibility for old "methods"
    if (method == "linear" || method == "loess") {
        return(invisible(do.call(retcor.peakgroups, alist(object, smooth=method, ...))))
    }

    method <- match.arg(method, getOption("BioC")$xcms$retcor.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("retcor", method, sep=".")

    invisible(do.call(method, alist(object, ...)))
})


setMethod("peakTable", "xcmsSet2", function(object, filebase = character(), ...) {

    if (length(sampnames(object)) == 1) {
        return(object@peaks)
    }

    if (nrow(object@groups) < 1) {
        stop ('First argument must be an xcmsSet with group information or contain only one sample.')
    }

    groupmat <- groups(object)


    if (! "value" %in% names(list(...))) {
        ts <- data.frame(cbind(groupmat,groupval(object, value="into",  ...)), row.names = NULL)
    } else {
        ts <- data.frame(cbind(groupmat,groupval(object, ...)), row.names = NULL)
    }

    cnames <- colnames(ts)

    if (cnames[1] == 'mzmed') {
        cnames[1] <- 'mz'
    } else {
        stop ('mzmed column missing')
    }
    if (cnames[4] == 'rtmed') {
        cnames[4] <- 'rt'
    } else {
        stop ('mzmed column missing')
    }

    colnames(ts) <- cnames

    if (length(filebase))
        write.table(ts, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    ts
})


#this modified version of diffreport takes the group names as an input so that they can be correlated back to the pre-diffreport results
#this is needed because some of the groups names may change after removing all non-labeled groups
#also suppresses output of .tsv file
setMethod("diffreport", "xcmsSet2", function(object, groupnames=NULL, class1 = levels(sampclass(object))[1],
                                            class2 = levels(sampclass(object))[2],
                                            filebase = character(), eicmax = 0, eicwidth = 200,
                                            sortpval = TRUE, classeic = c(class1,class2),
                                            value = c("into","maxo","intb"), metlin = FALSE,
                                            h = 480, w = 640, mzdec=2, ppm=50, ...) {

    if ( nrow(object@groups)<1 || length(object@groupidx) <1) {
        stop("No group information. Use group().")
    }

    require(multtest) || stop("Couldn't load multtest")

    value <- match.arg(value)
    groupmat <- groups(object)
    if (length(groupmat) == 0)
        stop("No group information found")
    samples <- sampnames(object)
    n <- length(samples)
    classlabel <- sampclass(object)
    classlabel <- levels(classlabel)[as.vector(unclass(classlabel))]

    values <- groupval(object, "medret", value=value)
    indecies <- groupval(object, "medret", value = "index")

    if (!all(c(class1,class2) %in% classlabel))
        stop("Incorrect Class Labels")

    ## c1 and c2 are column indices of class1 and class2 resp.
    c1 <- which(classlabel %in% class1)
    c2 <- which(classlabel %in% class2)
    ceic <- which(classlabel %in% classeic)
    if (length(intersect(c1, c2)) > 0)
        stop("Intersecting Classes")

    ## Check against missing Values
    if (any(is.na(values[,c(c1,c2)]))) {
        stop("NA values in xcmsSet. Use fillPeaks()")
    }

    mean1 <- rowMeans(values[,c1,drop=FALSE], na.rm = TRUE)
    mean2 <- rowMeans(values[,c2,drop=FALSE], na.rm = TRUE)

    ## Calculate fold change.
    ## For foldchange <1 set fold to 1/fold
    ## See tstat to check which was higher
    fold <- mean2 / mean1
    fold[!is.na(fold) & fold < 1] <- 1/fold[!is.na(fold) & fold < 1]

    testval <- values[,c(c1,c2)]
    testclab <- c(rep(0,length(c1)),rep(1,length(c2)))

    if (min(length(c1), length(c2)) >= 2) {
        tstat <- mt.teststat(testval, testclab, ...)
        pvalue <- pval(testval, testclab, tstat)
    } else {
        message("Too few samples per class, skipping t-test.")
        tstat <- pvalue <- rep(NA,nrow(testval))
    }
    stat <- data.frame(fold = fold, tstat = tstat, pvalue = pvalue)
    if (length(levels(sampclass(object))) >2) {
        pvalAnova<-c()
        for(i in 1:nrow(values)){
            var<-as.numeric(values[i,])
            ano<-summary(aov(var ~ sampclass(object)) )
            pvalAnova<-append(pvalAnova, unlist(ano)["Pr(>F)1"])
        }
        stat<-cbind(stat, anova= pvalAnova)
    }
    if (metlin) {
        neutralmass <- groupmat[,"mzmed"] + ifelse(metlin < 0, 1, -1)
        metlin <- abs(metlin)
        digits <- ceiling(-log10(metlin))+1
        metlinurl <- paste("http://metlin.scripps.edu/metabo_list.php?mass_min=",
                           round(neutralmass - metlin, digits), "&mass_max=",
                           round(neutralmass + metlin, digits), sep="")
        values <- cbind(metlin = metlinurl, values)
    }
    if(length(groupnames)>0)
		twosamp <- cbind(name = groupnames, stat, groupmat, values)
	else
		twosamp <- cbind(name = groupnames(object), stat, groupmat, values)
		
    if (sortpval) {
        tsidx <- order(twosamp[,"pvalue"])
        twosamp <- twosamp[tsidx,]
        rownames(twosamp) <- 1:nrow(twosamp)
        values<-values[tsidx,]
    } else
        tsidx <- 1:nrow(values)

    #if (length(filebase))
        #write.table(twosamp, paste(filebase, ".tsv", sep = ""), quote = FALSE, sep = "\t", col.names = NA)

    if (eicmax > 0) {
        if (length(unique(peaks(object)[,"rt"])) > 1) {
            ## This looks like "normal" LC data

            eicmax <- min(eicmax, length(tsidx))
            eics <- getEIC(object, rtrange = eicwidth*1.1, sampleidx = ceic,
                           groupidx = tsidx[seq(length = eicmax)],mzExpMeth="ppm",ppm=ppm)

            if (length(filebase)) {
                eicdir <- paste(filebase, "_eic", sep="")
                boxdir <- paste(filebase, "_box", sep="")
                dir.create(eicdir)
                dir.create(boxdir)
                if (capabilities("png")){
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="png",  width=w, height=h)
                    png(file.path(eicdir, "%003d.png"), width = w, height = h)
                } else {
                    xcmsBoxPlot(values[seq(length = eicmax),],
                                sampclass(object), dirpath=boxdir, pic="pdf", width=w, height=h)
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }
            plot(eics, object, rtrange = eicwidth, mzdec=mzdec)

            if (length(filebase))
                dev.off()
        } else {
            ## This looks like a direct-infusion single spectrum
            if (length(filebase)) {
                specdir <- paste(filebase, "_spec", sep="")
                dir.create(specdir)
                if (capabilities("png")){
                    png(file.path(specdir, "%003d.png"), width = w, height = h)
                }else{
                    pdf(file.path(eicdir, "%003d.pdf"), width = w/72,
                        height = h/72, onefile = FALSE)
                }
            }

            plotSpecWindow(object, gidxs = tsidx[seq(length = eicmax)], borderwidth=1)

            if (length(filebase))
                dev.off()
        }
    }

    invisible(twosamp)
})

setGeneric("diffreport2", function(object, ...) standardGeneric("diffreport2"))

#does all xcmsSet pre-processing and then leverages diffreport to produce a rough statistical report regarding all labeled compounds
#returns resultant table after writing it to a .csv file in the results directory 
#removed modPhenoData parameter
setMethod("diffreport2","xcmsSet2", function(object,resultsPath=NULL,type,tTest,pheno,
											 class1=NULL, class2=NULL, filebase = character(),
											 eicmax = 0, eicwidth = 200, sortpval = TRUE, classeic =NULL,
                                             value = c("into","maxo","intb"), dataVal = c("into","maxo","intb"),
											 metlin = FALSE, h = 480, w = 640, mzdec=2,ppm=50, ...)

{
	
	##Prep for calling diffreport
	classes = sampclass(object)
	classes = levels(classes)#[as.vector(unclass(classes))]
	
	#if not specified, set parameters to default values
	if(length(class1)==0)
		class1=classes[1]
		
	if(length(class2)==0)
		class2=classes[2]
		
	if(length(classeic)==0)  
		classeic = c(class1,class2)
	
	#read in valid Labelled Group Data
	if(type=="suspect")
	{
		stop("Type must be \"valid\", make sure validation of type \"valid\" has already been run")
	}else if(type=="valid")
	{
		if(tTest==T)
		{
			labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_tTestFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
			fname = paste(pheno,value,class1,"vs",class2,"tTestFilter_FinalReport.csv",sep="_")
			
		} else {
			labelData = read.csv(file=file.path(resultsPath,paste(pheno,type,"Clusters_HeuristicFilter",paste(dataVal,".csv",sep=""),sep="_")),check.names=F,stringsAsFactors=F)
			fname = paste(pheno,value,class1,"vs",class2,"HeuristicFilter_FinalReport.csv",sep="_")
		}
	}else 
	{
		warning("invalid type argument")
	}
	#Remove all groups and related information from xcmsSet2 that are not labelled Groups
	labeledGpIdx = match(labelData[ ,"groupname"], groupnames(object))
	object@groups = object@groups[labeledGpIdx, ]
	object@groupidx = object@groupidx[labeledGpIdx]
	#if(modPhenoData&&(ncol(object@phenoData)>1))
	#{
		#Remove all sample class info except highest level phenotype
	#	object@phenoData[ ,-1] = NULL
	#}
	
	
	cat("\n\nGenerating Report for:",class1,"and",class2,"\n")
	report = diffreport(object,groupnames = labelData[,"groupname"],class1=class1,class2=class2,filebase = filebase, eicmax = eicmax, eicwidth = eicwidth,
                                             sortpval = sortpval, classeic = classeic,value = value, metlin = FALSE,
                                             h = 480, w = 640, mzdec=2,ppm=50,...)

	##Add labeling data to diffreport 
	reportIdx = match(report[ ,"name"],labelData[ ,"groupname"])
	finalResults = cbind(labelData[reportIdx,1:2],report)
	
	write.csv(finalResults,file=file.path(resultsPath,fname),row.names=F)
	return (finalResults)

})



#setMethod("sampnames","xcmsSet2", function(object) gsub("\\.[A-Za-z]*$","",rownames(object@phenoData))) #updated sampnames to remove file extension from sample names
