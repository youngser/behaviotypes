## Author: youngser@jhu.edu
## Version: $Id: doidt.r,v 0.1 2012/02/30 20:55:01 parky Exp$
## Time-stamp: <Wed Oct 26, 2016 08:35:07 YP>
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
## DESCRIPTION: Code for Iterative Denoising Tree
##              NB: This code uses "mclust", a model-based clustering, but of course
##                  any clustering function can be used instead.
## USAGE:
##    out <- doIDT(feat,
##                 Dmax=ncol(feat), 
##                 maxsamp=nrow(feat),
##                 Kmax,  
##                 samp=1,
##                 maxdepth,
##                 minnum,
##                 pure=FALSE)    
##
##   feat:     n x d input feature matrix.
##   Dmax:     An integer specifying the maximum number of dimensions
##             for which the clustering is performed.
##   Kmax:     An integer specifying the maximum numbers of mixture components (clusters)
##             for which the BIC is to be calculated for Mclust.
##   maxsamp:  An integer specifying the maximum number of samples for clustering.
##   samp:     A logical or numeric vector specifying a subset of the data to be used
##             in the initial hierarchical clustering phase for Mclust.
##             samp=1 => no sampling, samp=2 => use n/2 for initialization, etc.
##   maxdepth: (termination conditions) maximum depth of the tree
##   minnum:   (termination conditions) minimum number of sample in branch
##   pure:     (termination conditions) if all samples are the same classes or not
##
## Output:  a list of leaf (cluster) info with repeated fields of followings:
##
##   depth:    depth of the leaf, maximum number of these will be the depth of the tree
##   bname:    branch name
##   ids:      sample ids for the leaf
##   isLeaf:   dummy element for counting the number of leaves

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

## main function
doIDT <- function(feat, Dmax=3, Kmax=2, maxsamp=200, samp=1,
                  maxdepth=5, minnum=50, pure=FALSE)
{
  ## initial depth == "root" == depth=0
  depth <- 0
  
  ## number of samples & dimensions
  if (!is.matrix(feat)) stop("Input data has to be a matrix!")
  n <- nrow(feat)
  d <- ncol(feat)
  dmax <- min(d,Dmax)

  ## sample ids & labels
  ids <- rownames(feat)
  if (is.null(ids)) ids <- 1:n
  blab <- rep(1,n)

  ## do idt
  idt <- dobranch(depth=depth, prev="", maxsamp=maxsamp,
                  bnum=1, bfeat=feat, bids=ids, blab=blab, 
                  dmax=dmax, Kmax=Kmax, samp=samp,
                  maxdepth=maxdepth, minnum=minnum, pure=pure)
  
  return(idt)
}

## function for handling branches
dobranch <- function(depth, prev="", maxsamp=maxsamp,
                     bnum, bfeat, bids, blab,
                     dmax=40, Kmax=2, samp=1,
                     maxdepth=5, minnum=1, pure=F)
{
  depth <- depth+1
  bname <- paste(prev,bnum,sep="")

  cat("===============================================\n")
  cat("Working on branch ",bname,", depth = ", depth, "(",maxdepth,")\n")

  labz <- blab
  blab <- labz[labz==bnum]
  bids <- bids[labz==bnum]
  zL <- matrix(bfeat[labz==bnum,],ncol=ncol(bfeat))
  vdmax <- min(dmax,ncol(zL))
  vKmax <- min(Kmax,nrow(zL))
  cat("n = ", nrow(zL), ", dim = ", ncol(zL), ", dmax = ", vdmax, ", Kmax = ", vKmax, "\n")

  ## stopping conditions
  leaf <- FALSE
  if(nrow(zL) <= minnum | depth > maxdepth | pure) {
    leaf <- TRUE
    cond1 <- pure
    cond2 <- nrow(zL) <= minnum
    cond3 <- depth > maxdepth
    cat("***** LEAF: ", bname, ": (pure,small,deep)=(",cond1,",",cond2,",",cond3,")\n")
    return(list(depth=depth,
                bname=bname,
                ids=bids,
                isLeaf=leaf))
  }

  ## reconfigure the data via pca
  set.seed(12345+bnum)
  zL <- prcomp(zL)$x
  elbow <- max(dim_select(apply(zL,2,sd)),vdmax)
  Svout <- zL[,1:elbow]

  cat("Clustering with dim = ",elbow)
  modelname <- c("EEE","EEV","VEV","VVV") ## avoid "I" models!
  sampsize <- min(maxsamp,floor(nrow(Svout)/samp))
  sub <- sample(nrow(Svout),size=sampsize)
  moutb <- Mclust(Svout,,G=1:vKmax,initialization=list(subset=sub),modelNames=modelname)

  khat <- moutb$G
  pure <- ifelse(khat==1,T,F)
  cat(", Khat: ", khat, ", ", moutb$modelName, "\n")

  ## recursion
  bidt <- NULL
  for (b in 1:moutb$G) {
    bout <- dobranch(depth,bname,maxsamp,
                     b,Svout,bids,moutb$class,
                     dmax,Kmax,samp,maxdepth,minnum,pure)
    bidt <- c(bidt,bout)
    cat("\n")
  }

  return(bidt)
}

## test function
testIDT <- function()
{
    ## install necessary packages
    required_packages  <- c("igraph","mclust")
    lapply(required_packages,installpkg)

    dat <- as.matrix(iris[,-5])
    lab <- iris[,5]
    
    out <- doIDT(as.matrix(dat),
                 Dmax=3, ## max dim for clustering
                 Kmax=3,  ## max K for clustering
                 maxsamp=nrow(dat), ## max n for clustering
                 samp=1,
                 maxdepth=5,
                 minnum=50,
                 pure=FALSE)    

    cat("Done!\n")
    cat("************************************************\n")
    idtlab <- rep(0,nrow(dat))
    nleaves <- sum(names(out)=="isLeaf")
    cat("number of leaves (clusters) = ",nleaves,"\n")
    nids <- which(names(out)=="ids")
    for (i in 1:nleaves) {
        idtlab[out[[nids[[i]]]]] <- i
    }

    print(table(idtlab,lab))
#      lab
#idtlab setosa versicolor virginica
#     1      0          3        50
#     2      0         47         0
#     3     50          0         0
    cat("\n ARI = ", adjustedRandIndex(idtlab,lab),"\n")
# ARI =  0.941045 

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
