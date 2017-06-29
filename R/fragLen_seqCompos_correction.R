
################################################################################
# Ludo Pagie, 2014-06-04, fragLen_seqCompos_correction.R
#
# DESCRIPTION:
#   a script to apply correction for the length and sequence composition of
#   GATC fragments to read coun data of single-cell DamID-seq data
#   Largely based on a previous script
#   (SingleCellDamID_preprocessing_BvS140218b.Rhtml) by Bas v. Steensel
#
# INPUT:
#   GATC.readcounts: a GRanges object with in the metadata object (mcols(...))
#     3 columns per sample; 1st with readcounts corresponding to the 5'-end of
#     the fragment, 2nd with readcounts of 3' end of fragment, 3rd with sum of
#     columns 1 and 2
#     Data of all samples is used. If samples have very low readcounts these
#     should be discarded beforehand.
#   GATC.fragments.mpbl: GRanges object with ....
#   binsize: integer, specifying the width of the genomic windows in which the
#     data is summarized
#
# OUTPUT:
#   A list with 5 elements:
#   - OE: a GRanges object with ranges covering the entire genome, each of
#     length 'binsize', and 1 column per sample in the mcols metadata slot. The
#     values in the mcols slot are the sum of the readcounts in overlapping
#     GATC fragments, normalized by the sum of the expected counts in the GATC
#     fragments 
#   - fragLenObs; a table object with the distribution of occurences of fragment lengths in the data
#   - fragLenObs.smth; same as fragLenObs but smoothed (numeric vector)
#   - fragLenExp; a tableobject with the distribution of occurences of fragment lengths in the mappable genome
#   - fragLenExp.smth; smoothed version (numeric vector) of fragLenExp
#   Using the last 4 elements plots of the Obs/Exp ratio can be plotted
#
# VERSIONS:
#   140604:
#     - first version
#   150112:
#     - small adaptation of function SummarizeReadcounts.new so that it takes
#       multiple binsizes, one for each chromosome. This adaptation is
#       necessary for the chromocounting project because we are only interested
#       in values per chromosome and not per bin with fixed size.


# Roxygen:
#' Preprocess readcounts in GATC.readcounts
#' 
#' \code{PreprocReadcounts} preprocesses readcounts, ie flatten counts to 0/1, etc
#' 
#' Details: extended description of function
#' 
#' @param GATC.readcounts GRanges object with GATC regions and readcounts (forw/rev/both) per sample.
#' 
#' @return GRanges object; granges specify genomic GATC fragments, mcols specify readcounts (forw/rev/both) per sample
#' 
#' 
#' 
#' 
#' 
PreprocReadcounts <- function(GATC.readcounts) {
  ##############################################################################
  # processing readcounts; 
  #   - remove columns containing sum of fwd and rev readcounts
  #   - coerce columns into Rle
  #   - flatten redcounts to [0,1]
  ##############################################################################
  # work only with read counts per fragment end; remove every 3rd column which
  # contains the sum of the fragment ends
  NCOL <- ncol(S4Vectors::mcols(GATC.readcounts))
  NROW <- nrow(S4Vectors::mcols(GATC.readcounts))
  idx.a <- seq(3, NCOL, by=3) # idx of columns with total counts per fragment (sum of forw and rev)
  # remove the columns which have the sums of the fragment ends
  S4Vectors::mcols(GATC.readcounts) <- S4Vectors::mcols(GATC.readcounts)[,-idx.a]
  NCOL <- ncol(S4Vectors::mcols(GATC.readcounts))

  # compress the Data frames for memory efficiency
  if (class(S4Vectors::mcols(GATC.readcounts)[,1]) != 'Rle') {
    for(i in seq.int(NCOL)) {
      S4Vectors::mcols(GATC.readcounts)[,i] <- S4Vectors::Rle(as.integer(S4Vectors::mcols(GATC.readcounts)[,i]))
    }
  }

  # flatten the readcounts per fragment end
  for(i in seq.int(NCOL)) {
    S4Vectors::runValue(S4Vectors::mcols(GATC.readcounts)[,i]) <- as.integer(S4Vectors::runValue(S4Vectors::mcols(GATC.readcounts)[,i]) > 0)
  }

  # now mcols(GATC) has values [0,1] to indicate detection of a
  # fragment. forw and reverse reads are in separate (subsequent)
  # columns. columns containing sum of forw and reverse are removed
  return(GATC.readcounts)
}

# Roxygen:
#' Compute OE scores from read counts
#' 
#' \code{SummarizeReadcounts} computes OE scores from read counts
#' 
#' Details: extended description of function
#' 
#' @param fExp Numeric vector, length 2x number of GATC fragments, specifies
#'   for each possible forward and reverse read whether it is expected (ie,
#'   mappable).
#' @param GATC.readcounts GRanges object with GATC regions and readcounts (forw/rev/both) per sample.
#' @param binsize Integer, either a single value or a vector of same length as number of chromosomes.
#' 
#' @return GRanges object; granges specify genomic regions of size \code{bin.size}, mcols specify computed OE scores
#' 
SummarizeReadcounts <- function(fExp, GATC.readcounts, binsize) {
  ##############################################################################
  # get (reformatted) readcounts from GRanges object
  ##############################################################################
  GATC.readcounts <- PreprocReadcounts(GATC.readcounts)

  ##############################################################################
  # make sure al readcounts for fragment ends which have 0 fExp are set to 0
  # fExp is 0 if:
  # - read is unmappable
  # - read is from too long fragment
  # - no fragment of this length or no fragment with this AT content was observed
  #     in data (and smoothing of fExp didn't results in fExp > 0 either)
  # Only the latter case may result in 'true' read counts being discarded; highly
  #   unlikely (unless the fExp was based on severly undersampled data)
  ##############################################################################
  tolerance <- 1e-10
  # fExp is a vector with length = 2xnrow(GATC.readcounts); 2 elements
  # per fragment, 1 forward and 1 reverse.
  # GATC.readcounts has same setup for columns
  # forw reads:
  forw.logi <- c(TRUE, FALSE)
  S4Vectors::mcols(GATC.readcounts)[abs(fExp[forw.logi]) < tolerance, forw.logi] <- 0L
  # reverse reads:
  rev.logi <- c(FALSE, TRUE)
  S4Vectors::mcols(GATC.readcounts)[abs(fExp[rev.logi]) < tolerance, rev.logi] <- 0L
  # cleanup 
  rm(tolerance, forw.logi, rev.logi)

  #############################################
  # Calculate genome-wide read depth per sample
  #############################################
  # calculate genome-wide frequency of detected fragments per cell, (ie 'read
  # depth' after flattening the read counts):
  # The read depth is scaled by the number of fragments
  sample.depth <- sapply(seq(1,ncol(S4Vectors::mcols(GATC.readcounts)),by=2),
                         function(i) sum(S4Vectors::mcols(GATC.readcounts)[,i]>0) + 
                                     sum(S4Vectors::mcols(GATC.readcounts)[,i+1]>0)
                        )
  sample.depth <- sample.depth/nrow(S4Vectors::mcols(GATC.readcounts))

  ############################
  # check structure of binsize
  ############################
  # added LP150112
  # binsize should either be a single integer (original situation) or an
  # integer vector with a length equal to the number of chromosomes in fExp and
  # GATC.readcounts. The values in binsize should be less than or equal to the
  # length of each chromosome
  # If binsize is a single integer than it should be multiplied into a vector
  # with a length equal to the number of chromosomes in the other arguments.
  if (length(binsize) !=1 & length(binsize) != length(GenomeInfoDb::seqlengths(GATC.readcounts)))
    stop("binsize in 'SummarizeReadcounts.new(..)' should be of length 1 or of same length as the number of chromosomes")
  if (length(binsize) == 1)
    binsize <- rep(binsize, times=length(GenomeInfoDb::seqlengths(GATC.readcounts)))
  chr.ends <- GenomeInfoDb::seqlengths(GATC.readcounts)
  if (!is.null(names(binsize)) & !all(names(binsize) %in% names(chr.ends)))
    stop("in 'SummarizeReadcounts.new(..)': binsize should either have no names or the names should be the same as in GATC.readcounts")
  if (is.null(names(binsize)))
    names(binsize) <- names(chr.ends)
  # at this point binsize should be a named vector with names equal to the sequence names in GATC.readcounts
  # end LP150112

  ###############################################################################################
  # compute the final observed over expected readcounts for all genomic windows of size 'binsize'
  ###############################################################################################
  # A GATC fragment will be assigned to a bin if it's central base pair overlaps with the bin.
  # Hence, create new GRange object GATC.fragments.mid that contains only center positions of F:
  GATC.fragments.mid <- GenomicRanges::GRanges(GenomeInfoDb::seqnames(GATC.readcounts), 
                                IRanges::ranges(GATC.readcounts), 
                                seqinfo=GenomeInfoDb::seqinfo(GATC.readcounts))
  IRanges::start(GATC.fragments.mid) <- IRanges::end(GATC.fragments.mid) <- IRanges::mid(IRanges::ranges(GATC.fragments.mid))
  # make bins; suppressWarnings() because I know bins are defined outside
  # sequence bounds, which is why they are trimmed in the last statement
  # LP150112
  # chr.ends <- seqlengths(GATC.readcounts)
  chr.bins <- lapply(names(chr.ends), function(n) seq(from=1, to=chr.ends[n], by=binsize[n])); names(chr.bins) <- names(chr.ends)
  chr.bins <- lapply(names(chr.bins), function(n) IRanges::IRanges(start=chr.bins[[n]], width=binsize[n])); names(chr.bins) <- names(chr.ends)
  chr.bins <- as(do.call(IRanges::IRangesList, chr.bins),"GRanges")
  # chr.bins <- lapply(chr.ends, seq, from=1, by=binsize)
  # chr.bins <- lapply(chr.bins, IRanges, width=binsize)
  # chr.bins <- as(do.call(IRangesList, chr.bins),"GRanges")
  # LP150112 end
  suppressWarnings(GenomeInfoDb::seqinfo(chr.bins) <- GenomeInfoDb::seqinfo(GATC.readcounts))
  # temporarily set circularity of chrM to FALSE for function trim(..) to work as expected
  suppressWarnings(GenomeInfoDb::isCircular(chr.bins)[TRUE] <- FALSE) 
  chr.bins <- IRanges::trim(chr.bins)
  GenomeInfoDb::isCircular(chr.bins) <- GenomeInfoDb::isCircular(GATC.readcounts)
  # clean up
  rm(chr.ends)

  # get overlap of GATC fragments and genomic bins
  ovl <- GenomicRanges::findOverlaps(query=chr.bins, subject=GATC.fragments.mid)
  # sum expected (relative) readcounts per bin; fExp has 1 element per
  # fragment-end: expand ovl by factor 2
  binexp <- aggregate(fExp, by=list(rep(S4Vectors::queryHits(ovl), each=2)), FUN=sum)$x

  qhit.rle <- S4Vectors::Rle(S4Vectors::queryHits(ovl))
  binobs.agg <- sapply(seq.int(ncol(S4Vectors::mcols(GATC.readcounts))), 
                       function(idx) S4Vectors::aggregate.Rle(S4Vectors::mcols(GATC.readcounts)[,idx], by=qhit.rle, FUN=sum))
  # sum forw columns and reverse columns
  binobs.agg <- binobs.agg[,c(TRUE,FALSE)] + binobs.agg[,c(FALSE,TRUE)]

  # create a matrix to store the obs/Exp values
  OEmat <- matrix(NA,nrow=length(chr.bins), ncol=ncol(S4Vectors::mcols(GATC.readcounts))/2, 
                  dimnames=list(NULL,colnames(S4Vectors::mcols(GATC.readcounts))[c(TRUE,FALSE)]))
  OEmat[S4Vectors::runValue(qhit.rle),] <- binobs.agg/binexp
  OEmat <- t(t(OEmat)/sample.depth)

  # create GRanges object to hold the OE values by making a copy of
  # chr.bins
  OE <- chr.bins
  S4Vectors::mcols(OE) <- S4Vectors::DataFrame(OEmat, check.names=FALSE)
  # and copy the sample meta data (but only every 2nd row)
  S4Vectors::mcols(S4Vectors::mcols(OE)) <- S4Vectors::mcols(S4Vectors::mcols(GATC.readcounts))[c(TRUE,FALSE),,drop=FALSE]
  # correct column names of mcols(OE)
  names(S4Vectors::mcols(OE)) <- sub('_forw$', '', names(S4Vectors::mcols(OE)))

  # return OE
  return (OE=OE)
}
###############
