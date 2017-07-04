# Roxygen:
#' @title Import bam files and generate read counts per GATC fragment
#' 
#' @description
#' \code{BamToFragCounts} generates read counts per GATC fragment for a set of bam files.
#' 
#' @details
#' extended description of function
#' 
#' @param meta.data Data frame.
#' @param bam.dir Character, directory which contains the bamfiles.
#' @param GATC.fragments GRanges object with GATC regions.
#' @param maxgap Integer, maximum allowed gap between read alignment position
#' and GATC position.
#' @param shift Integer, exact number of bases by which reads are shifted
#' before comparing to GATC position.
#' @param VERBOSE Logical, print diagnostic msgs [default = FALSE].
#' 
#' @return GRanges object; granges specify genomic GATC fragments, mcols
#' specify readcounts (forw/rev/both) per sample
#' 
#' 
#' 
#' 
#' 
Bam2FragCounts <- 
  function(meta.data, bam.dir, GATC.fragments, maxgap=0,
			   shift=0, VERBOSE=FALSE) {
  # Ludo Pagie, 140521, Bam2FragCounts,
  #
  # a function to read a set of bam files containing reads from DamID-seq
  # experiment. The reads are aligned to GATC sites in the genome, and counted
  # per GATC fragment. The result is returned a sa GRanges object.
  #
  # NOTES:
  #   This function treats circular chromosomes as regular non-circular
  #   chromosomes. This is necessary to compute overlap between reads and
  #   fragments
  #
  # ARGUMENTS:
  #   sample.file: filename of sample meta file containing at least a column
  #     named 'file' which contains for each sample the name of the bamfile
  #   bam.dir: path of the directory containing the bam files
  #   GATC.fragments: GRanges object containing the all possible fragments from
  #     which the reads originate. The fragments are expected to include the
  #     restriction motif (eg, GATC) at both the start and the end of the
  #     fragments. Reads are expected to align to exactly the start (or end, in
  #     case of the reverse reads) of the fragment. If reads are expected to
  #     align some bases away from the start/end of the fragment use argument
  #     'shift' to adjust this.
  #   maxgap:integer, specifies the extend by which a read can be offset relative
  #     to the expected alignment position, due to, for instance, inaccurate
  #     ligation of adapter sequences.
  #   shift: integer, specifiying the fixed number of positions by which the
  #     reads need to be shifted in order to align precisely with the start/end
  #     of the fragments (see also 'GATC.fragments')
  #
  # RETURN VALUE:
  #   a GRanges object which is a copy of the user supplied GATC.fragments
  #   object. The returned object in addition contains a column (in the
  #   mcols(GATC.fragments) part) for each sample with the number of reads per
  #   fragment in that sample


  bam.files <- as.character(meta.data$bam.fname)
  sample.names <- as.character(meta.data$sample.name)
  # check if all bam files exist
  notfound <- which( !file.exists(file.path(bam.dir, bam.files)) )
  if(length(notfound) > 0)
    stop(sprintf("some bam files could not be found:\n%s\naborting\n", 
                 paste(bam.files[notfound], sep='\n')))
  # read bamfiles. readGAlignmentsFromBam returns objects of class GAlignments.
  # These objects can readily be used in findOverlaps function
  # mapreads <- lapply(file.path(bam.dir,bam.files), readGAlignmentsFromBam)
  mapreads <- lapply(file.path(bam.dir,bam.files), GenomicAlignments::readGAlignments)
  names(mapreads) <- sample.names

  ## if reads need to be shited relative to the start of the GATC site (shift !=
  ## 0) do it dependent on the strand of the read
  if (shift != 0)
    mapreads <- sapply(mapreads, function(elm) GenomicRanges::shift(elm, shift*ifelse(GenomicRanges::strand(elm) == '+', 1, -1)))

  # mapping the reads to whole GATC fragments, either to start or end depending on the strand
  GenomeInfoDb::isCircular(GATC.fragments)[TRUE] <- FALSE # take care of potential circular chromosomes
  # coerce GATC fragments into an interval tree to speed up findOverlaps 
  # GATC.tree <- GIntervalTree(GATC.fragments)
  GATC.gncl <- GenomicRanges::GNCList(GATC.fragments)
  # for all samples count reads flanking GATC in forw or reverse
  # direction
  for (nm in sample.names) {
    # print(paste("processing sample",nm))
    cat("\r",paste("processing sample",nm))
    # separate reads in strand direction
    fwd.str.idx <- GenomicRanges::strand(mapreads[[nm]]) == '+'
    S4Vectors::mcols(GATC.fragments)[,paste(nm, 'forw',sep='_')] <- 
      GenomicRanges::countOverlaps(GATC.gncl, mapreads[[nm]][fwd.str.idx], type='start', maxgap=maxgap)
    S4Vectors::mcols(GATC.fragments)[,paste(nm, 'rev',sep='_')] <- 
      GenomicRanges::countOverlaps(GATC.gncl, mapreads[[nm]][!fwd.str.idx], type='end', maxgap=maxgap)
    S4Vectors::mcols(GATC.fragments)[,nm] <- 
      S4Vectors::mcols(GATC.fragments)[,paste(nm, 'forw',sep='_')] +
      S4Vectors::mcols(GATC.fragments)[,paste(nm, 'rev',sep='_')]
  }
  cat("\r\n")

  # add sample meta data to 'GATC.fragments'
  # there are 3 columns per sample:
  meta.data <- meta.data[rep(ncol(meta.data), each=3,drop=FALSE),]
  # add column with forw/rev/both direction
  meta.data$direction = c('forw','rev','both')
  S4Vectors::mcols(S4Vectors::mcols(GATC.fragments)) <- meta.data

  # return
  return (GATC.fragments)
}
