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
		    VERBOSE=FALSE, do.norm=TRUE) {

  # first, check input
  if ( ( length(bam.fname) != length(sample.name) ) || ( length(bam.fname) != length(cell.count) ) )
    stop("length of vectors 'bam files', 'experiment names', and 'cell counts' is not equal, aborting\n\n")

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
  GATC.fragments.readcounts <- Bam2FragCounts(meta.data=md, bam.dir=".", GenomicRanges::granges(GATC.mpbl.reads), maxgap=0, shift=0, VERBOSE=VERBOSE)
  ofname <- file.path(outdir,sprintf("%s_GATC-readcounts_%s.RData", exp.name, sig))
  if(VERBOSE) print(sprintf("Saving GATC.fragments.readcounts to %s", ofname))
  save(file=ofname, x=GATC.fragments.readcounts)
  if (VERBOSE) print("done bam2fragcounts")

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

  # calculate total read count per sample (discard duplicate reads)
  depth <- sapply(rep(seq(1,ncol(S4Vectors::mcols(GATC.fragments.readcounts)),by=3),each=2)+0:1, # sum detected fragments in every 1st and 2nd column, ie discard every 3rd column. 3rd column contains the sum of 1st and 2nd but impossible to discard duplicates per forw/rev read.
		  function(i) sum(S4Vectors::mcols(GATC.fragments.readcounts)[,i]>0))
  depth <- tapply(depth, rep(1:(ncol(S4Vectors::mcols(GATC.fragments.readcounts))/3), each=2), sum) # sum counts of forw and rev reads per GATC
  # scale read depth by the total of possible reads, and logscale
  depth <- log10(depth/sum(S4Vectors::mcols(GATC.mpbl.reads)$nreads))

  # generate image with barplot of readcounts
  ofname <- file.path(outdir, sprintf("%s_readcounts_%s.png", exp.name, sig))
  # bitmap(file=ofname, res=144, taa=4)
  png(file=ofname, res=144, width=4*480, height=3*480)
  barplot(-depth, col=1:5, main='-log10(Readdepth) as proportion of mappable reads', ylab='proportion', xlab='sample index')
  dev.off()
  ## optionally; discard samples with relative readdepth < 10^-4
  if (any (! (depth>-4) ) ) {
    ok.idx <- which(depth>-4)
    ok.idx <- rep(ok.idx-1, each=3)*3+1:3
    S4Vectors::mcols(GATC.fragments.readcounts) <- S4Vectors::mcols(GATC.fragments.readcounts)[,ok.idx]
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
  ofname <- file.path(outdir, sprintf("%s_OE_raw_profiles_%s.png", exp.name, sig))
  # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
  png(file=ofname, res=144, width=4*480, height=3*480)
  opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat))),2), mar=c(0.5,1,0.5,1))
  for (i in seq.int(ncol(OE.mat))) {
    plot(OE.mat[,i], main='', ylim=c(0,2), axes=FALSE, pch=19)
    axis(2)
    box()
  }
  par(opar)
  dev.off()
  # generate image of sample ID and cell count, as reference
  ofname <- file.path(outdir, sprintf("%s_OE_profiles_mdata_%s.png", exp.name, sig))
  # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
  png(file=ofname, res=144, width=4*480, height=3*480)
  opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat))),2), mar=c(0.5,1,0.5,1))
  for (i in seq.int(ncol(OE.mat))) {
    plot(NA, main='', xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, pty='n')
    text(0,0,labels=sprintf("%s\n%d", colnames(OE.mat)[i], md$cell.count[i]), cex=1)
    # axis(2)
    box()
  }
  par(opar)
  dev.off()

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
    ofname <- file.path(outdir, sprintf("%s_OE_norm_profiles_%s.png", exp.name, sig))
    # bitmap(file=ofname, res=144, taa=4, width=28, height=21)
    png(file=ofname, res=144, width=4*480, height=3*480)
    opar  <- par(mfrow=rep(ceiling(sqrt(ncol(OE.mat.norm))),2), mar=c(0.5,1,0.5,1))
    for (i in seq.int(ncol(OE.mat.norm))) {
      plot(OE.mat.norm[,i], main="", ylim=c(0,2), axes=FALSE, pch=19)
      axis(2)
      box()
    }
    par(opar)
    dev.off()
  
    # coerce OE.mat.norm to GRanges object
    OE.norm <- OE
    S4Vectors::mcols(OE.norm) <- S4Vectors::DataFrame(OE.mat.norm)
  }

  return(list(OE.norm.factor=OE.norm.factor, OE.norm=OE.norm, OE.raw=OE))
}
