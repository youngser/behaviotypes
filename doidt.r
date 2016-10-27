## Author: youngser@jhu.edu
## Version: $Id: doidt.r,v 0.1 2012/02/30 20:55:01 parky Exp$
## Time-stamp: <Thu Oct 27, 2016 08:13:51 YP>
#####################################################################################
## Copyright 2014 Youngser Park
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##    http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing,
## software distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and limitations
## under the License.
#####################################################################################
##
##------------------------------------------------------------------------------
##' Run Iterative Denoising Tree
##'
##' @aliases doidt
##' @references
##' C.E. Priebe, D.J. Marchette, and D.M. Healy,
##' "Integrated Sensing and Processing Decision Trees,"
##' IEEE Transactions on Pattern Analysis and Machine Intelligence,
##' Vol. 26, No. 6, pp. 699-708, 2004.
##'
##' @note This code uses "mclust", a model-based clustering, but of course
##'     any clustering function can be used instead.
##'
##' @examples
##'    out <- doIDT(feat,
##'                 FUN="mclust",
##'                 Dmax=ncol(feat),
##'                 maxsamp=nrow(feat),
##'                 Kmax,
##'                 samp=1,
##'                 maxdepth,
##'                 minnum,
##'                 pure=FALSE,
##'                 verbose=FALSE)
##'
##' @param   feat:     n x d input feature matrix.
##' @param   FUN:      clustering function, has to be either "mclust" or "pamk" for now.
##' @param   Dmax:     An integer specifying the maximum number of dimensions
##'                  for which the clustering is performed.
##' @param   Kmax:     An integer specifying the maximum numbers of mixture components (clusters)
##'                   for which the BIC is to be calculated for Mclust.
##' @param   maxsamp:  An integer specifying the maximum number of samples for clustering.
##' @param   samp:     A logical or numeric vector specifying a subset of the data to be used
##'                   in the initial hierarchical clustering phase for Mclust.
##'                   samp=1 => no sampling, samp=2 => use n/2 for initialization, etc.
##' @param   maxdepth: (termination conditions) maximum depth of the tree
##' @param   minnum:   (termination conditions) minimum number of sample in branch
##' @param   pure:     (termination conditions) if all samples are the same classes or not
##' @param   verbose:  intermediate progress message
##'
##' @return  Returns a list of leaf (cluster) info with repeated fields of followings:
##'
##' @return    depth:    depth of the leaf, maximum number of these will be the depth of the tree
##' @return    bname:    branch name
##' @return    ids:      sample ids for the leaf
##' @return    isLeaf:   dummy element for counting the number of leaves
##'
##' @author Youngser Park \email{youngser@jhu.edu}


## main function
doIDT <- function(feat,FUN="mclust",Dmax=3, Kmax=2, maxsamp=200, samp=1,
                  modelnames=c("EEE","EEV","VEV","VVV"), knee=2, scree=FALSE,
                  maxdepth=5, minnum=50, pure=FALSE, scale=FALSE, verbose=FALSE)
{
    ## initial depth == "root" == depth=0
    depth <- 0

    ## number of samples & dimensions
    if (!is.matrix(feat)) stop("Input data has to be a matrix!")
    n <- nrow(feat)
    d <- ncol(feat)
    dmax <- min(d,Dmax)
#    feat <- feat[,1:dmax]

    ## sample ids & labels
    ids <- rownames(feat)
    if (is.null(ids)) ids <- 1:n
    blab <- rep(1,n)

    ## do idt
    idt <- dobranch(feat=feat, depth=depth, FUN=FUN, prev="", maxsamp=maxsamp,
                    bnum=1, bfeat=feat, bids=ids, blab=blab,
                    modelnames=modelnames,
                    dmax=dmax, Kmax=Kmax, samp=samp, scale=scale,
                    knee=knee, scree=scree,
                    maxdepth=maxdepth, minnum=minnum, pure=pure, verbose=verbose)

    idtlab <- rep(0,n)
    leaf.nodes <- pkg.env$idtall[which(sapply(pkg.env$idtall, function(x) x$isLeaf))]
    (nleaves <- length(leaf.nodes))
    cat("number of leaves (clusters) = ",nleaves,"\n")
    nids <- sapply(leaf.nodes, "[","ids")
    for (i in 1:nleaves) idtlab[nids[[i]]] <- i

    return(list(classification=idtlab,idtall=pkg.env$idtall))
}

#idtidx <- 1
#idtall <- NULL
pkg.env <- new.env()
pkg.env$idtidx <- 1
pkg.env$idtall <- NULL


## function for installing needed packages
installpkg <- function(x) {
    if (x %in% rownames(installed.packages())==FALSE) {
        if (x %in% rownames(available.packages())==FALSE) {
            paste(x,"is not a valid package - please check again...")
        } else {
            install.packages(x,dep=TRUE)
            require(x,character.only = TRUE)
        }

    } else {
        paste(x,"package already installed...")
        require(x,character.only = TRUE)
    }
}

## function for handling branches
dobranch <- function(feat, depth, FUN="mclust",prev="", maxsamp=maxsamp,
                     bnum, bfeat, bids, blab, modelnames="VVV",
                     dmax=40, Kmax=2, samp=1, scale=FALSE, knee=2, scree=FALSE,
                     maxdepth=5, minnum=1, pure=F, verbose=TRUE)
{
    depth <- depth+1
    bname <- paste(prev,bnum,sep="")

    if (verbose) {
        cat("===============================================\n")
        cat("Working on branch ",bname,", depth = ", depth, "(",maxdepth,")\n")
    }

    labz <- blab
    blab <- labz[labz==bnum]
    bids <- bids[labz==bnum]
    zL <- matrix(bfeat[labz==bnum,],ncol=ncol(bfeat))
    feat <- feat[labz==bnum,]
    vdmax <- min(dmax,ncol(zL))
    vKmax <- min(Kmax,nrow(zL))
    if (verbose) {
        cat("n = ", nrow(zL), ", dim = ", ncol(zL), ", dmax = ", vdmax, ", Kmax = ", vKmax, "\n")
    }

    ## stopping conditions
    leaf <- FALSE
    if(nrow(zL) <= minnum | depth > maxdepth | pure) {
        leaf <- TRUE
        cond1 <- pure
        cond2 <- nrow(zL) <= minnum
        cond3 <- depth > maxdepth
        if (verbose) {
            cat("***** LEAF: ", bname, ": (pure,small,deep)=(",cond1,",",cond2,",",cond3,")\n")
        }

        pkg.env$idtall[[pkg.env$idtidx]] <- list(depth=depth,
                                                 bname=bname,
                                                 pname=prev,
                                                 ids=bids,
                                                 isLeaf=leaf)
        pkg.env$idtidx <- pkg.env$idtidx + 1

        return()
    }

    ## reconfigure the data via pca
    set.seed(12345+bnum)
#    zL <- prcomp(zL,scale=scale)$x
    zL <- prcomp(feat,scale=scale)$x
#    elbow <- max(dim_select(apply(zL,2,sd)),vdmax)
#    elbow <- dim_select(apply(zL,2,sd))
    elbow <- getElbows(zL,2,plot=scree,main=bname)
    if (knee > 2) stop("knee cannot be > 2!")
    if (length(elbow) > 1) elbow <- elbow[knee]
    if (verbose) cat("Clustering in dim = ",elbow)
    Svout <- zL[,1:elbow,drop=FALSE]

#    modelname <- c("EEE","EEV","VEV","VVV") ## avoid "I" models!
    if (maxsamp < nrow(Svout)) {
        sampsize <- min(maxsamp,floor(nrow(Svout)/samp))
        sub <- sample(nrow(Svout),size=sampsize)
    } else {
        sub <- 1:nrow(Svout)
    }

    switch(FUN,
           kmeans2={
               kout <- kmeans(Svout,2)
               khat <- 2
               bclass <- kout$cluster
               if (verbose) cat(", Khat: ", khat, "\n")
           },
           mclust={
               moutb <- Mclust(Svout,G=1:vKmax,initialization=list(subset=sub),
                               modelNames=modelnames)
               khat <- moutb$G
               bclass <- moutb$class
               if (verbose) cat(", Khat: ", khat, ", ", moutb$modelName, "\n")
           },
           pamk={
               pamout <- pamk(Svout,krange=1:vKmax,usepam=FALSE,critout=FALSE)
               khat <- pamout$nc
               bclass <- pamout$pamobject$clustering
               if (verbose) cat(", Khat: ", khat, "\n")
           },
           stop("Enter either mclust or pamk!")
           )

    pure <- ifelse(khat==1,T,F)

    ## recursion
    bidt <- NULL
    for (b in 1:khat) {
        bout <- dobranch(feat,depth,FUN,bname,maxsamp,
                         b,Svout,bids,bclass,modelnames=modelnames,dmax,
                         Kmax,samp,scale,knee,scree,maxdepth,minnum,pure,verbose)
        bidt <- c(bidt,bout)
#        bidt <- append(bidt,bout)
        if (verbose) cat("\n")
    }

    pkg.env$idtall[[pkg.env$idtidx]] <- list(depth=depth,
                                             bname=bname,
                                             pname=prev,
                                             ids=bids,
                                             isLeaf=leaf)
    pkg.env$idtidx <- pkg.env$idtidx + 1
    return()
}

## test function
testIDT <- function()
{
    FUN <- "mclust"
#    FUN <- "kmeans2"
#    FUN <- "pamk"

    ## install necessary packages
    required_packages  <- c("igraph","mclust","fpc")
    lapply(required_packages,installpkg)

    dat <- as.matrix(iris[,-5])
    lab <- iris[,5]

    out <- doIDT(as.matrix(dat),
                 FUN=FUN,
                 Dmax=3, ## max dim for clustering
                 Kmax=3,  ## max K for clustering
                 maxsamp=nrow(dat), ## max n for clustering
#                 modelnames=c("EEE","EEV","VEV","VVV"),
                 samp=1,
                 maxdepth=5,
                 minnum=50,
                 pure=FALSE,
                 verbose=TRUE)

    cat("Done!\n")
    cat("************************************************\n")
    ## idtlab <- rep(0,nrow(dat))
    ## nleaves <- sum(names(out)=="isLeaf")
    ## cat("number of leaves (clusters) = ",nleaves,"\n")
    ## nids <- which(names(out)=="ids")
    ## for (i in 1:nleaves) {
    ##     idtlab[out[[nids[[i]]]]] <- i
    ## }
    idtlab <- out$class
    cat("maximum depth (root=1) = ",max(sapply(out[[2]], "[[", 1)),"\n")
    cat("number of leaves (clusters) = ",max(idtlab),"\n")

    print(table(idtlab,lab))
#      lab
#idtlab setosa versicolor virginica
#     1      0          3        50
#     2      0         47         0
#     3     50          0         0
## idtlab setosa versicolor virginica # pamk
##      1     32          0         0
##      2     18          1         0
##      3      0         23        13
##      4      0         26         1
##      5      0          0        36

    cat("\n ARI = ", adjustedRandIndex(idtlab,lab),"\n")
# ARI =  0.941045
# ARI =  0.5324059 pamk

    ## NB: without IDT
    cat("\n** With no IDT:\n")
    mout <- Mclust(dat)
    cat("number of clusters = ", mout$G,"\n")
    print(table(mout$class, lab))
#   lab
#    setosa versicolor virginica
#  1     50          0         0
#  2      0         50        50
    cat("\n ARI = ", adjustedRandIndex(mout$class,lab),"\n")
#[1] 0.5681159

}

getElbows <- function(dat, n = 3, threshold = FALSE, plot = TRUE, main="") {
  ## Given a decreasingly sorted vector, return the given number of elbows
  ##
  ## Args:
  ##   dat: a input vector (e.g. a vector of standard deviations), or a input feature matrix.
  ##   n: the number of returned elbows.
  ##   threshold: either FALSE or a number. If threshold is a number, then all
  ##   the elements in d that are not larger than the threshold will be ignored.
  ##   plot: logical. When T, it depicts a scree plot with highlighted elbows.
  ##
  ## Return:
  ##   q: a vector of length n.
  ##
  ## Reference:
  ##   Zhu, Mu and Ghodsi, Ali (2006), "Automatic dimensionality selection from
  ##   the scree plot via the use of profile likelihood", Computational
  ##   Statistics & Data Analysis, Volume 51 Issue 2, pp 918-930, November, 2006. 

#  if (is.unsorted(-d))


  if (is.matrix(dat)) {
      d <- sort(apply(dat,2,sd), decreasing=TRUE)
  } else {
      d <- sort(dat,decreasing=TRUE)
  }

  if (!is.logical(threshold))
      d <- d[d > threshold]
  
  p <- length(d)
  if (p == 0)
      stop(paste("d must have elements that are larger than the threshold ",
                 threshold), "!", sep="")

  lq <- rep(0.0, p)                     # log likelihood, function of q
  for (q in 1:p) {
      mu1 <- mean(d[1:q])
      mu2 <- mean(d[-(1:q)])              # = NaN when q = p
      sigma2 <- (sum((d[1:q] - mu1)^2) + sum((d[-(1:q)] - mu2)^2)) /
          (p - 1 - (q < p))
      lq[q] <- sum( dnorm(  d[1:q ], mu1, sqrt(sigma2), log=TRUE) ) +
          sum( dnorm(d[-(1:q)], mu2, sqrt(sigma2), log=TRUE) )
  }
  
  q <- which.max(lq)
  if (n > 1 && q < (p-1)) {
      q <- c(q, q + getElbows(d[(q+1):p], n-1, plot=FALSE))
  }
  
  if (plot==TRUE) {
      if (is.matrix(dat)) {
          sdv <- d # apply(dat,2,sd)
          plot(sdv,type="b",xlab="dim",ylab="stdev",main=main)
          points(q,sdv[q],col=2,pch=19)
      } else {
          plot(dat, type="b", main=main)
          points(q,dat[q],col=2,pch=19)
      }
  }
  
  return(q)
}
