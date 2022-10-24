## Identify potential mispriming sites by aligning sequences along ITR-seq peak sites to primer sequence
MatchGenomicLocusToPrimer <- function(gref, gchr, gstart, gend, pseq, base_ext=2*nchar(pseq[[1]]),
                                      misp.req=list(match.start=1, match.base=10, match.percent=0.3, match.perfect=FALSE)) {
  # gref      BSgenome object of reference genome
  # gchr      A vector of characters with chromosome names of the genomic at the junction sites
  # gstart    A vector of integer with the start positions of the genomic at the junction sites; same length as gchr
  # gend      A vector of integer with the end positions of the genomic at the junction sites; same length as gchr
  # pseq      Full primer sequence
  # base_ext  No. of bases to be extended on both sides
  # misp.req  Conditions to call a mispriming
  #           match.start: the start position of the match in the primer; position from the end if strand is '-'
  #           match.base: minimum number of matching bases
  #           match.percent: minimum percent of matching bases relative to full primer length
  #           match.perfect: whether mismatch and INDEL are allowed in the matching
  
  require(Biostrings);
  require(GenomicRanges);

  ##############################################################################################
  ## Run pairwiseAlignment to find the best local alignment b/w primer and genomic sequences
  bestLocalAlign <- function(gref, rng, s0, str) {
    # gref: BSgenome of reference genome
    # rng:  Genomic regions
    # s0:   Primer sequence
    # str:  Strand to align to; reverse-complement primer sequence if '-'
    
    #######################
    ## internal function ##
    summarizeAlignRange <- function(pa, rng) {
      rng0 <- pa@subject@range;
      rng1 <- pa@pattern@range;
      
      # Align genomic sequences to primer to report best local alignment
      if (str != '-') pos0 <- cbind(pstart=start(rng0), pend=end(rng0)) else
        pos0 <- cbind(pstart=nchar(s)-end(rng0)+1, pend=nchar(s)-start(rng0)+1);
      pos1 <- cbind(gstart=start(rng)+start(rng1)-1, gend=start(rng)+end(rng1)-1);
      
      # Summarize alignment
      tbl <- data.frame(score=score(pa), pos0, pos1, stringsAsFactors = FALSE);
      tbl$total <- nchar(pa);
      tbl$match <- nmatch(pa);
      tbl$mismatch <- nmismatch(pa);
      tbl$indel <- nindel(pa)@insertion[, 2] + nindel(pa)@deletion[, 2];
      tbl$primer <- as.character(subject(pa));
      tbl$genomic <- as.character(pattern(pa));
      if (str == '-') {
        tbl$primer <- as.character(reverseComplement(DNAStringSet(tbl$primer)));
        tbl$genomic <- as.character(reverseComplement(DNAStringSet(tbl$genomic)));
      }
      rownames(tbl) <- names(rng);
      
      tbl;
    }
    ###############
    
    s0 <- DNAStringSet(s0);
    s1 <- getSeq(gref, rng);
    
    if (str != '-') s <- s0 else s <- reverseComplement(s0); 
    pa  <- pairwiseAlignment(s1, s, type='local', gapOpening=2);
    tbl <- summarizeAlignRange(pa, rng);
    
    ###################################################
    ## Fine-tuning alignment to the 1st base of primer

    ## Allow 1 base gap in the primer
    sub0 <- subseq(rep(s, length(s1)), rep(1, length(s1)), end(pa@subject@range));
    stt1 <- pmax(1, start(pa@pattern@range) - start(pa@subject@range));
    end1 <- end(pa@pattern@range);
    rng1 <- rng;
    start(rng1) <- start(rng) + stt1 - 1;
    end(rng1) <- start(rng1) + (end1 - stt1);
    sub1 <- subseq(s1, stt1, end1);
    pwa1 <- pairwiseAlignment(sub1, sub0, gapOpening=0);
    tbl1 <- summarizeAlignRange(pwa1, rng1);

    ## Allow 1 base gap in the genomic sequence
    sub0 <- subseq(rep(s, length(s1)), rep(1, length(s1)), end(pa@subject@range))
    stt2 <- pmax(1, start(pa@pattern@range) - start(pa@subject@range) + 2);
    end2 <- end(pa@pattern@range);
    rng2 <- rng;
    start(rng2) <- start(rng) + stt2 - 1;
    end(rng2) <- start(rng2) + (end2 - stt2);
    sub1 <- subseq(s1, stt2, end2)
    pwa2 <- pairwiseAlignment(sub1, sub0, gapOpening=0);
    tbl2 <- summarizeAlignRange(pwa2, rng2);

    ## Better alignment with 1 base added or removed
    tbl0 <- tbl1;
    tbl0[tbl2[, 1]>=tbl1[, 1], ] <- tbl2[tbl2[, 1]>=tbl1[, 1], ];

    ## quite arbitrary selection criteria
    whh0 <- tbl0[, 1]>=10 & tbl0$pstart==1 & tbl0$match>tbl$match & (tbl0$match/tbl0$total)>=0.8;
    tbl[whh0, ] <- tbl0[whh0, , drop=FALSE];
    ###################################################
    
    tbl; 
  }
  ##############################################################################################
  
  ## Re-format primer sequence
  pseq <- DNAString(as.character(pseq)[1]);
  
  ## Reconcile difference in chromosome names
  chr0 <- unique(gchr);
  chr1 <- seqlevels(gref); 
  names(chr0) <- sub('MT', 'M', sub('^chr', '', chr0, ignore.case=TRUE)); ## Simpliest form of chromosome names
  names(chr1) <- sub('MT', 'M', sub('^chr', '', chr1, ignore.case=TRUE));
  cmap <- sapply(names(chr0), function(cnm) as.vector(chr1[names(chr1)==cnm])[1]);
  cmap <- cmap[!is.na(cmap)];
  
  ## Retrieve genomic sequences at given loci
  rng <- GRanges(cmap[gchr], IRanges(gstart-base_ext, gend+base_ext));
  names(rng) <- 1:length(rng);

  ## Align genomic sequences to both strands of primer sequence; report the best local alignment
  aln1 <- bestLocalAlign(gref, rng, pseq, '+');
  aln2 <- bestLocalAlign(gref, rng, pseq, '-');
  
  ## Summary alignment results and use the strand with higher alignment score
  bstr <- c('+', '-')[1+as.integer(aln2$score>aln1$score)]; # which strand has higher score
  aln <- aln1;
  aln[bstr=='-', ] <- aln2[bstr=='-', , drop=FALSE];
  aln <- data.frame(chr=gchr, start=gstart, end=gend, strand=bstr, aln, stringsAsFactors=FALSE); # done with alignment
  
  ###############################
  ## Determine mispriming site ##
  ###############################
  misp <- rep(0, nrow(aln));  # Whether a mispriming
  
  ## Requirement 1: Best alignment starts before given base position in primer (default is the 1st base)
  req1 <- rep(FALSE, nrow(aln));
  req1 <- aln$pstart<=misp.req$match.start;

  ## Requirement 2:  Alignment is long enough to include no fewer than given number of matching bases (default is 10)
  req2 <- aln$match>=misp.req$match.base;
  
  ## Requirement 3:  Alignment is long enough to include no fewer than given % of full primer length (default is 30%)
  req3 <- (aln$match/nchar(pseq)) >= misp.req$match.percent;
  
  ## Requirement 4 (optional): Alignment is perfect (no mismatch or INDEL)
  if (misp.req$match.perfect) req4 <- aln$total==aln$match else req4 <- rep(TRUE, nrow(aln));
  
  ## Sites satisfying all requrements
  reqs <- req1 & req2 & req3 & req4;
  
  ## Whether the alignment positions on genome is overlapping the input junction sites
  olap <- aln$start<=aln$gend & aln$end>=aln$gstart;

  # 1 - meet all the requirements, and the alignment at genomic region overlapping the input loci
  # 2 - meet all the requirements, but the alignment at genomic region not overlapping the input loci
  misp[reqs &  olap] <- 2;
  misp[reqs & !olap] <- 1;
  
  dist <- pmin(abs(aln$end-aln$gstart), abs(aln$start-aln$gend));
  dist[aln$end>=aln$gstart & aln$start<=aln$gend] <- 0;
  
  aln$distance   <- dist;
  aln$mispriming <- misp;
  
  invisible(aln);
};



###############
## DEMO code ##
###############

## Rhesus 
# require(BSgenome.Mmulatta.UCSC.rheMac10);
# pseq <- readRDS('/Users/zhangz/Documents/GTP/Project/2022-08_ITRseq2_Diag/R/B050/primer_seq.rds');
# cln1 <- readRDS('/Users/zhangz/Documents/GTP/Project/2022-08_ITRseq2_Diag/R/B050/filtered_clone.rds')[[1]];
# gref <- Mmulatta;
# gchr <- cln1[, 1];
# gstart <- cln1[, 2];
# gend <- cln1[, 3];
# base_ext <- 2*nchar(pseq);
# misp.req <- list(match.start=1, match.base=10, match.percent=0.3, match.perfect=FALSE)

## Canine
require(BSgenome.Cfamiliaris.Ensembl.107);
pseq <- readRDS('/Users/zhangz/Documents/GTP/Project/2022-08_ITRseq2_Diag/R/B050/primer_seq.rds'); # Primer sequence
# cln1 <- readRDS('/Users/zhangz/Documents/GTP/Project/2022-08_ITRseq2_Diag/R/B407/filtered_clone.rds')[[1]];
gref <- Cfamiliaris;
gchr <- cln1[, 1];
gstart <- cln1[, 2];
gend <- cln1[, 3];
base_ext <- 2*nchar(pseq);
misp.req <- list(match.start=1, match.base=10, match.percent=0.3, match.perfect=FALSE)

aln <- MatchGenomicLocusToPrimer(gref, gchr, gstart, gend, pseq, base_ext, misp.req);
