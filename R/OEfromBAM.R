#' @importFrom grDevices bitmap dev.off png

#' @importFrom graphics axis barplot box par plot text

#' @importFrom stats lm predict

# Roxygen:
#' @title Compute OE scores from (a set of) bamfiles.
#' 
#' @description
#' \code{bamToOE} computes the observed-over-expected (OE) scores for a set of samples. The input is a bamfile for each sample.
#' 
#' @details
#' extended description of function
#' 
#' @param exp.name Character, used to create output file names.
#' @param bam.dir Character, directory which contains the bamfiles.
#' @param bam.fname Character vector, bamfile names for all samples.
#' @param sample.name Character vector, (short) names for all samples.
#' @param cell.count Integer vector, cell count for all samples.
#' @param sig Character, string with initials and date used for output file names (eg LP170619).
#' @param CHR Character vector, chromosome names to be used in the analysis.
#' @param bin.size integer, 0 (each bin is entire chr
#' length) or >0 specifying the width in bp of the bins used for computing the
#' OE scores.
#' @param fExp Numeric vector, length of twice the number of GATC fragments in
#' \code{GATC.mpbl.reads} (2 reads per fragment) specifying the (normalized)
#' expected readcount per fragment.
#' @param GATC.mpbl.reads GRanges object, granges specify genomic GATC fragments, mcols specifies whether forw/rev/both read count (0/1/2).
#' @param OE.norm.factor Numeric vector, normalization factor per window/chromosome to normalize raw OE scores.
#' @param outdir Character, specifies directory for output files.
#' @param VERBOSE Logical, print diagnostic msgs [default = FALSE].
#' @param do.norm Logical, perform normalization [default = FALSE].
#' @param plot Logical, generate plots [default = TRUE].
#' @param addToName Character, a substrings which will be added in fromt of the
#' chromsome names read from the bamfiles. This allows dealing with bamfiles
#' which lack the 'chr' in the name.
#' @param min.depth integer, minimal read depth per sample (samples with less reads are discarded) [default = 10000].
#' @param sample.depth integer, target read depth per sample (samples will be downsampled to target; must be >= \code{min.depth}; sampling is \emph{without} replacement; set \code{sample.depth} to 0 (default) to prevent sampling) [default = 0].
#' 
#' @return List with 3 elements;
#' \describe{
#'   \item{OE.norm.factor}{Numeric vector; normalization factor per window/chromosome to normalize raw OE scores}
#'   \item{OE.norm}{GRanges object; granges specify genomic GATC fragments, mcols specify normalized OE scores for each sample (ie each bamfile)}
#'   \item{OE.raw}{GRanges object; granges specify genomic GATC fragments, mcols specify raw OE scores for each sample (ie each bamfile)}
#' }
#' 
#' @examples
#' \dontrun{
#' res <- bamToOE(exp.name="myExperiment", bam.dir="./bamfiles/",
#'   bam.fname=c("1.bam","2.bam"), sample.name=c("smpl1","smpl2"),
#'   cell.count=c(20,1),CHR=paste0('chr', c(1:19, 'X')), fExp=fExp.mm10,
#'   GATC.mpbl.reads=GATC.mpbl.reads.mm10, outdir="myResults")
#' }
#' 
#' @export
bamToOE <- function(exp.name, bam.dir=getwd(), bam.fname, sample.name, cell.count, 
		    sig=sprintf("LP%s", format(Sys.time(), "%y%m%d")), CHR, 
		    bin.size=0, fExp, GATC.mpbl.reads, OE.norm.factor=NULL, outdir, 
		    VERBOSE=FALSE, do.norm=TRUE, plot=TRUE, plotCHR=FALSE, addToName="",
		    min.depth=1e3, sample.depth=0L) {

  # first, check input
  if ( ( length(bam.fname) != length(sample.name) ) || ( length(bam.fname) != length(cell.count) ) )
    stop("length of vectors 'bam files', 'experiment names', and 'cell counts' is not equal, aborting\n\n")

  bam.dir <- normalizePath(bam.dir)
  if (all (file.exists (bam.fname))) {
    bam.fname <- bam.fname 
  } else {
    if ( all (file.exists (file.path (bam.dir, bam.fname)))) {
      bam.fname <- file.path(bam.dir, bam.fname)
    } else {
      if (all (file.exists (file.path (bam.dir, basename(bam.fname))))) {
	bam.fname <- file.path (bam.dir, basename(bam.fname))
      } else {
	stop("(some) bamfiles do not exist, aborting\n\n")
      }
    }
  }
  # check existence of output directory, or create if not existing yet
  if (!dir.exists(outdir)) dir.create(outdir)

  ord <- order(sample.name)
  md <- S4Vectors::DataFrame(bam.fname=bam.fname[ord], sample.name=sample.name[ord], cell.count=cell.count[ord])

  if (VERBOSE) {
	print("processing bam files with following metadata:")
	print(md)
  }

  # run bam2Frag
  if (VERBOSE)
    print("starting bam2fragcounts")
  GATC.fragments.readcounts <- Bam2FragCounts(meta.data=md, bam.dir="", GenomicRanges::granges(GATC.mpbl.reads), maxgap=0, shift=0, VERBOSE=VERBOSE, addToName)
  # calculate total read count per sample (discard duplicate reads)
  depth <- sapply(rep(seq(1,ncol(S4Vectors::mcols(GATC.fragments.readcounts)),by=3),each=2)+0:1, # sum detected fragments in every 1st and 2nd column, ie discard every 3rd column. 3rd column contains the sum of 1st and 2nd but impossible to discard duplicates per forw/rev read.
		  function(i) sum(S4Vectors::mcols(GATC.fragments.readcounts)[,i]>0))
  depth <- tapply(depth, rep(1:(ncol(S4Vectors::mcols(GATC.fragments.readcounts))/3), each=2), sum) # sum counts of forw and rev reads per GATC
  names(depth) <- S4Vectors::mcols(S4Vectors::mcols(GATC.fragments.readcounts))$sample.name[c(T,F,F)]

  if(plot) {
    # generate image with barplot of readcounts
    ofname <- file.path(outdir, sprintf("%s_readcounts_%s.png", exp.name, sig))
    # bitmap(file=ofname, res=144, taa=4)
    png(filename=ofname, res=144, width=4*480, height=3*480)
    barplot(-depth, col=1:5, main='-log10(Readdepth) as proportion of mappable reads', ylab='proportion', xlab='sample index')
    dev.off()
  }

  ## optionally; discard samples with relative readdepth < 10^-4
  if (any (depth<min.depth)) {
    ok.idx <- which(! depth<min.depth)
    if (VERBOSE)
      print(sprintf("Some samples have too low read depth: \n%s", paste(depth[-ok.idx, drop=FALSE], collapse="\n", sep="")))
    ok.idx <- rep(ok.idx-1, each=3)*3+1:3
    S4Vectors::mcols(GATC.fragments.readcounts) <- S4Vectors::mcols(GATC.fragments.readcounts)[,ok.idx]
  }

  # downsample reads if reqeusted
  if (sample.depth) {
    if (VERBOSE)
      print(sprintf("downsampling reads to %.0e reads", sample.depth))
    # checks
    if (sample.depth < min.depth)
      stop("sample.depth < min.depth\nAborting")
    sample.depth <- as.integer(sample.depth)
    if (sample.depth<0)
      stop(sprintf("sample.depth should be an integer >= 0 (currently: %s)\nAborting",as.character(sample.depth)))
    if(any(depth<sample.depth))
      stop(sprintf("some samples have: sum(reads ) < sample.depth\ncurrent readdepth:\n%s\nAborting", paste(depth, sep="", collapse="\n")))

    # sampling; GATC.fragments.readcounts is a GRanges object, with readcounts
    # per GATC fragment in forw/rev/both columns per sample. Sampling needs to
    # be done at level of forw/rev column. I will make vectors with (repeated)
    # indices, corresponding to the GATC fragments with positive readcounts.
    # From this vector I will sample. The columns will then be re-populated via
    # tabulating the sampled index vector. Forw and rev reads will be
    # discriminated by taking the negative for the rev reads.
    n.samples <- length(S4Vectors::mcols(GATC.fragments.readcounts))/3
    # extract count matrix
    # readcounts.mat <- as.matrix(GenomicRanges::as.data.frame(S4Vectors::mcols(GATC.fragments.readcounts)))
    # iterate over samples
    for (s in seq.int(n.samples)) { 
      # column indices for forw/rev/both columns 
      forw.ind <- (s-1)*3+1
      rev.ind <- (s-1)*3+2
      both.ind <- s*3
      # collect row-numbers of readcount matrix and repeat them according to readcount
      # mark readcounts of rev column by taking negative
      inds <- c(rep(seq.int(length(GATC.mpbl.reads)), S4Vectors::mcols(GATC.fragments.readcounts)[,forw.ind]), 
		rep(-(seq.int(length(GATC.mpbl.reads))), S4Vectors::mcols(GATC.fragments.readcounts)[,rev.ind]))
      # inds <- c(rep(seq.int(length(forw.ind)), readcounts.mat[,forw.ind]), rep(-(seq.int(length(rev.ind))), readcounts.mat[,rev.ind]))
      # sample from row-numbers
      inds.smpl <- sample(inds, size=sample.depth, replace=FALSE)
      # tabulate sampled row-numbers
      cnts <- table(inds.smpl)
      # split forw and rev readcounts
      cnts.forw <- cnts[!grepl("^-",names(cnts))]
      cnts.rev <- cnts[grepl("^-",names(cnts))]
      # store sampled readcounts in matrix
      # first, initialize count matrix to 0
      S4Vectors::mcols(GATC.fragments.readcounts)[,c(forw.ind,rev.ind, both.ind)] <- 0
      S4Vectors::mcols(GATC.fragments.readcounts)[as.integer(names(cnts.forw)),forw.ind] <- cnts.forw
      S4Vectors::mcols(GATC.fragments.readcounts)[-as.integer(names(cnts.rev)),rev.ind] <- cnts.rev
      S4Vectors::mcols(GATC.fragments.readcounts)[,both.ind] <- 
	S4Vectors::mcols(GATC.fragments.readcounts)[,forw.ind] + 
	S4Vectors::mcols(GATC.fragments.readcounts)[,rev.ind]
#      readcounts.mat[T] <- 0
#      readcounts.mat[as.integer(names(cnts.forw)),forw.ind] <- cnts.forw
#      readcounts.mat[-as.integer(names(cnts.rev)),rev.ind] <- cnts.rev
#      readcounts.mat[,both.ind] <- readcounts.mat[,forw.ind] + readcounts.mat[,rev.ind]
    }
    # replace readcounts by downsampled counts
    # mcols(GATC) <- S4Vectors::DataFrame(readcounts.sampled.mat)
  }
  ofname <- file.path(outdir,sprintf("%s_GATC-readcounts_%s.RData", exp.name, sig))
  if(VERBOSE) print(sprintf("Saving GATC.fragments.readcounts to %s", ofname))
  save(file=ofname, x=GATC.fragments.readcounts)
  if (VERBOSE) print("done bam2fragcounts")
  # scale read depth by the total of possible reads, and logscale, historic reason for subsequent computations below
  depth <- log10(depth/sum(S4Vectors::mcols(GATC.mpbl.reads)$nreads))


  ################## ###
  # compute OE scores #
  #####################
  # set bin size of calculating OE scores
  if (bin.size == 0) {
    bin.size <- GenomeInfoDb::seqlengths(GATC.fragments.readcounts)
  } else {
    if ( any( bin.size>GenomeInfoDb::seqlengths(GATC.fragments.readcounts)[CHR] ) )
      stop("bin.size larger than a chromosomelength")
  }

  # coerce readcounts vectors to Rle (to speedup computations)
  for(i in seq.int(ncol(S4Vectors::mcols(GATC.fragments.readcounts)))) {
    S4Vectors::mcols(GATC.fragments.readcounts)[,i] <- S4Vectors::Rle(as.integer(S4Vectors::mcols(GATC.fragments.readcounts)[,i]))
  }
  # compute OE scores
  if (VERBOSE) print(paste("compute OE"))
  OE <- SummarizeReadcounts(fExp, GATC.fragments.readcounts, bin.size)
  # discard sequences not in CHR
  OE <- OE[as.character(GenomeInfoDb::seqnames(OE)) %in% CHR]
  GenomeInfoDb::seqlevels(OE) <- GenomeInfoDb::seqlevelsInUse(OE)
  # save OE scores
  ofname <- file.path(outdir, sprintf("%s_OE_raw_%s.RData", exp.name, sig))
  save(file=ofname, OE)

  ##########################
  ## plot raw OE profiles ##
  ##########################
  # extract data matrix from OE
  OE.mat <- as.matrix(GenomicRanges::as.data.frame(S4Vectors::mcols(OE)))
  # the following expression should not be necessary
  rownames(OE.mat) <- paste(as.character(GenomeInfoDb::seqnames(OE)), as.integer(GenomicRanges::start(OE)), sep="_")
  colnames(OE.mat) <- as.character(names(GenomicRanges::mcols(OE)))

  # generate image with profiles of all OE scores
  if(plot) {
    ofname <- file.path(outdir, sprintf("%s_OE_raw_profiles_%s.png", exp.name, sig))
    # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
    png(filename=ofname, res=144, width=4*480, height=3*480)
    opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat))),2), mar=c(0.5,1,0.5,1))
    for (i in seq.int(ncol(OE.mat))) {
      plot(OE.mat[,i], main='', ylim=c(0,2), axes=FALSE, pch=19)
      axis(2)
      box()
    }
    par(opar)
    dev.off()
  }

  # generate image of sample ID and cell count, as reference
  if (plot) {
    ofname <- file.path(outdir, sprintf("%s_OE_profiles_mdata_%s.png", exp.name, sig))
    # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
    png(filename=ofname, res=144, width=4*480, height=3*480)
    opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat))),2), mar=c(0.5,1,0.5,1))
    for (i in seq.int(ncol(OE.mat))) {
      plot(NA, main='', xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, pty='n')
      text(0,0,labels=sprintf("%s\n%d", colnames(OE.mat)[i], md$cell.count[i]), cex=1)
      # axis(2)
      box()
    }
    par(opar)
    dev.off()
  }
  # 
  if (plot && plotCHR) {
    CHRS <- paste0('chr',c(1:19,'X'))
    EXPS <- as.character(names(GenomicRanges::mcols(OE)))
    names <- sub("(.*RMA.*_)0(.*)_.*","\\1\\2", EXPS)
    names(names) <- EXPS

    for (exp in EXPS) {
      ofname <- file.path(outdir, sprintf("%s_OE_raw_profiles_%s_CHRS_%s.png", exp.name, exp, sig))
      # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
      png(filename=ofname, res=144, width=4*480, height=3*480)
      opar  <- par(mfrow=c(5,5), mar=c(0.5,1,0.5,1))
      for (chr in CHRS) {
	prof <- GenomicRanges::mcols(OE)[[exp]][as.character(GenomeInfoDb::seqnames(OE)) == chr]
	plot(prof, main=sprintf("%s; chr%s", names[exp], chr), ylim=c(0,2), axes=FALSE, pch=19)
	axis(2)
	box()
      }
      par(opar)
      dev.off()
    }
  }


  if (!do.norm) {
    OE.norm.factor <- NULL
    OE.norm <- NULL
  } else {
    #################
    # NORMALIZATION #
    #################
    # if null model of OE scores is not given, compute here
    if (is.null(OE.norm.factor)) {
      # compute avg profile on basis of 'trimmed' subset of OE scores
      # trimming is based on following thresholds: median(data) +/- 3*mad(data)
      OE.mat.cleaned <- OE.mat
      for (r in 1:nrow(OE.mat)) {
        median <- median(OE.mat[r,])
        mad <- mad(OE.mat[r,])
        idx <- which(OE.mat[r,] < median-3*mad | OE.mat[r,] > median+3*mad)
        OE.mat.cleaned[r,idx] <- NA
      }
      OE.null <- rowMeans(OE.mat.cleaned, na.rm=TRUE)
      rm(OE.mat.cleaned)
  
      # compute AT per GATC fragment
      seq.compos <- Biostrings::letterFrequency(
             x=IRanges::Views(BSgenome.Mmusculus.UCSC.mm10::Mmusculus, 
        		            GenomicRanges::granges(GATC.mpbl.reads)), 
  	     letters=c('AT','CG'))/IRanges::width(GATC.mpbl.reads)
      # compute mean seq.compos per unit (ie chromosome of window)
      # 1st compute overlap between OE bins and GATC fragments 
      #   (here seq.compos per OE bin is the mean of seq.compos of all
      #   GATC.fragments overlapping with the OE.bin, which means some GATC
      #   fragments may overlap with multiple OE.bins)
      ovl <- GenomicRanges::findOverlaps(GenomicRanges::granges(GATC.mpbl.reads), GenomicRanges::granges(OE))
      seq.compos.binned <- aggregate(seq.compos[S4Vectors::queryHits(ovl),1], list(S4Vectors::subjectHits(ovl)), mean)$x
      OE.norm.factor <- predict(lm(y~x, data.frame(y=OE.null,x=seq.compos.binned)),
                                newdata=data.frame(x=seq.compos.binned))
    }
    OE.mat.norm <- OE.mat/OE.norm.factor
  
    # save normalized OE data
    ofname <- file.path(outdir, sprintf("%s_OE_norm_%s.RData", exp.name, sig))
    save(file=ofname, OE.mat.norm)

    # generate profiles of normalized OE scores, for all samples
    if (plot) {
      ofname <- file.path(outdir, sprintf("%s_OE_norm_profiles_%s.png", exp.name, sig))
      # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
      png(filename=ofname, res=144, width=4*480, height=3*480)
      opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat.norm))),2), mar=c(0.5,1,0.5,1))
      for (i in seq.int(ncol(OE.mat.norm))) {
	plot(OE.mat.norm[,i], main="", ylim=c(0,2), axes=FALSE, pch=19)
	axis(2)
	box()
      }
      par(opar)
      dev.off()
    }
  
    # coerce OE.mat.norm to GRanges object
    OE.norm <- OE
    S4Vectors::mcols(OE.norm) <- S4Vectors::DataFrame(OE.mat.norm)
  }

  return(list(OE.norm.factor=OE.norm.factor, OE.norm=OE.norm, OE.raw=OE))
}
