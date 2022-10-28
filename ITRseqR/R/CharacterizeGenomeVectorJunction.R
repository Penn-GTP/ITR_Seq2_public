#### This function characterizes sequences at the genome-vector junction sites where integration happens, by reads
#### It requires a bam file of paired reads aligned to the host genome and reference sequences
#### This 1st version assumes insertion starts within ITR and 5ITR and 3ITR are perfect RevComp to each other
CharacterizeGenomeVectorJunctionRead <- function(fbam, gref, vref, pitr) {
  # fbam:   bam file of paired reads aligned to the reference genome
  # gref:   an BSgenome object of reference genome sequence
  # vref:   an DNAString object of vector sequence
  # pitr:   3 integer positions in vector: start of 5ITR, end of 3ITR, and ITR length
  #         ex. p0146 c(1, 2801, 168)
  
  require(Biostrings);
  require(GenomicAlignments);
  
  ######################################################################################################
  ## Read in all read pairs with both ends aligned to the reference genome and all of their SAM fields
  rds <- readGAlignmentPairs(fbam, use.names = TRUE, param = ScanBamParam(what=scanBamWhat()));
  read.count <- list(original=length(rds)); # No. of reads at each step
  ######################################################################################################
  
  ## Match chromosome names between bam file and reference genome
  chr0 <- seqlevels(rds);             # original chromosome names in bam file
  chr1 <- seqlevels(gref);            # original chromosome names in reference genome
  names(chr0) <- sub('MT', 'M', sub('^chr', '', chr0, ignore.case=TRUE)); ## Simpliest form of chromosome names
  names(chr1) <- sub('MT', 'M', sub('^chr', '', chr1, ignore.case=TRUE));
  cmap <- sapply(names(chr0), function(cnm) as.vector(chr1[names(chr1)==cnm])[1]);
  cmap <- cmap[!is.na(cmap)];
  cat(length(cmap), 'chromosomes in the bam file were found in the reference genome.\n');
  rds <- rds[seqnames(rds) %in% names(cmap)];  # Remove reads aligned to chromosomes not found in the reference genome
  seqlevels(rds) <- as.vector(cmap);        # Rename chromosome names of aligned reads
  read.count$chromosome.found <- length(rds);
  
  ######################################################################################################
  ######################################################################################################
  #### Filtering read pairs
 
  ## 1. select properly paired reads (same chromosome)
  chr1 <- as.vector(seqnames(rds@first)); # same chromosome
  chr2 <- as.vector(seqnames(rds@last));
  str1 <- as.vector(strand(rds@first));   # opposite strands
  str2 <- as.vector(strand(rds@last));
  stt1 <- start(rds@first);               # correct order along chromosome
  stt2 <- start(rds@last);
  end1 <- end(rds@first);
  end2 <- end(rds@last);
  flg1 <- chr1==chr2 & str1=='+' & str2=='-' & end1<=end2;
  flg2 <- chr1==chr2 & str1=='-' & str2=='+' & stt1>=stt2;
  rds <- rds[flg1 | flg2];
  read.count$proper.pair <- length(rds);
  
  ## 2. select read pairs with high mapping quality
  min.mapq <- 32;
  mpq1 <- rds@first@elementMetadata$mapq;
  mpq2 <- rds@last@elementMetadata$mapq;
  rds <- rds[mpq1>=min.mapq & mpq2>=min.mapq];
  read.count$high.mapq <- length(rds);
  
  ## 3. select read pair with soft-clipping in Read 2 (at the beginning/end if Read 2 mapped to plus/minus strand
  ## [Based on the experiment protocol that captures ITR at the ends of vector and goes outwards into genome]
  str  <- as.vector(strand(rds@last));
  cig  <- cigar(rds@last);
  flg1 <- str=='+' & grepl('^[0-9]+S', cig);  # softclipping at the end
  flg2 <- str=='-' & grepl('[0-9]+S$', cig);  # softclipping at the end
  rds <- rds[flg1 | flg2];
  read.count$soft.clipped <- length(rds);
  
  ## 4. remove duplicates if read pairs aligned to multiple locations 
  if (length(rds) > length(unique(names(rds)))) {
    tmp <- rds;
    tmp <- tmp[order(tmp@last@elementMetadata$flag)];      # keep one with the lowest Read 2 SAM flag value
    tmp <- tmp[rev(order(tmp@last@elementMetadata$mapq))]; # keep one with the highest Read 2 mapping quality
    tmp <- tmp[!duplicated(names(tmp))]; # remove duplicated read pairs
    rds <- rds[names(rds) %in% names(tmp)];
  }
  read.count$duplicate.removed <- length(rds);
  
  #### Done with filtering
  ######################################################################################################
  ######################################################################################################
  
  ## Prepare summary table of alignment
  fst <- rds@first;
  lst <- rds@last;
  smm <- data.frame(Chromosome=as.vector(seqnames(rds)), 
                    StrandR1=as.vector(strand(fst)), StrandR2=as.vector(strand(lst)),
                    StartR1=start(fst), EndR1=end(fst), StartR2=start(lst), EndR2=end(lst),
                    ReadLenR1=nchar(fst@elementMetadata$seq), ReadLenR2=nchar(lst@elementMetadata$seq), 
                    CigarR1=cigar(fst), CigarR2=cigar(lst), stringsAsFactors = FALSE);
  rownames(smm) <- names(rds);
  
  
  ######################################################################################################  
  ######################################################################################################
  #### Map Read 2 and genome at the junction to to vector and find the best local match
  
  ##############################
  ## Map Read 2 to vector
  seq <- rds@last@elementMetadata$seq; # Full sequences of Read 2
  dup <- duplicated(seq);              # Collase reads with the same sequence to save time
  seq <- seq[!dup];

  ###############################################################################
  ## Find the best match in Read 2 to vector sequence b/w 2 ITR
  pwa <- pairwiseAlignment(seq, substr(vref, pitr[1], pitr[2]), type='local');
  ###############################################################################
  # [Not reverse reads mapped to "+" strand as the 2 ITR are RevCom to each other
  # And we are only interested in the ITR integration here]
  
  cnt0 <- nmatch(pwa);
  cnt1 <- nmismatch(pwa) + nindel(pwa)@insertion[, 2] + nindel(pwa)@deletion[, 2];
  
  rng0 <- pwa@subject@range;  # Alignment location in the vector 
  rng1 <- pwa@pattern@range;  # Alignment location in Read 2 
  
  ## Summary of local best matches per read
  r2v <- cbind(start0=start(rng0), end0=end(rng0), start1=start(rng1), end1=end(rng1), length=cnt0, mismatch=cnt1);
  rownames(r2v) <- as.character(seq);
  r2v <- r2v[as.character(rds@last@elementMetadata$seq), ];
  rownames(r2v) <- names(rds);
  
  ## Re-format and transform the read-to-vector table
  rtbl <- r2v;
  # Transform vector positions if R2 mapped to plus strand
  rtbl[smm$StrandR2=='+', 1] <- pitr[2] + 1 - rtbl[smm$StrandR2=='+', 1]; 
  rtbl[smm$StrandR2=='+', 2] <- pitr[2] + 1 - rtbl[smm$StrandR2=='+', 2];
  rtbl[smm$StrandR2=='+', 1:2] <- rtbl[smm$StrandR2=='+', 2:1, drop=FALSE]; 
  # Transform reads positions if R2 mapped to plus strand
  rtbl[smm$StrandR2=='+', 3] <- smm[smm$StrandR2=='+', 'ReadLenR2'] + 1 - rtbl[smm$StrandR2=='+', 3]; 
  rtbl[smm$StrandR2=='+', 4] <- smm[smm$StrandR2=='+', 'ReadLenR2'] + 1 - rtbl[smm$StrandR2=='+', 4]; 
  rtbl[smm$StrandR2=='+', 3:4] <- rtbl[smm$StrandR2=='+', 4:3, drop=FALSE]; 
  colnames(rtbl) <- paste0('R2V_', c('StartV', 'EndV', 'StartR2', 'EndR2', 'NMatch', 'NMismatch'));
  ##############################

  ##############################
  ## Map genome sequence at the junction to vector
  grs <- GRanges(seqnames(rds), IRanges(start(rds@last), end(rds@last)), strand=strand(rds));
  ref <- getSeq(gref, grs);            # Retrieve junction sequences from reference genome
  dup <- duplicated(ref);              # Collapse same sequences to save time
  seq <- ref[!dup];

  ##############################################################################
  ## Find the best match in reference sequence to vector sequence b/w 2 ITR
  pwa <- pairwiseAlignment(seq, substr(vref, pitr[1], pitr[2]), , type='local');
  ##############################################################################
  
  cnt0 <- nmatch(pwa);
  cnt1 <- nmismatch(pwa) + nindel(pwa)@insertion[, 2] + nindel(pwa)@deletion[, 2];
  
  rng0 <- pwa@subject@range;  # Alignment location in the vector 
  rng1 <- pwa@pattern@range;  # Alignment location in Read 2 
  
  ## Summary of local best matches per read
  g2v <- cbind(start0=start(rng0), end0=end(rng0), start1=start(rng1), end1=end(rng1), length=cnt0, mismatch=cnt1);
  rownames(g2v) <- as.character(seq);
  g2v <- g2v[as.character(ref), ];
  rownames(g2v) <- names(rds);
  
  ## Genomic coordinates of the best local match
  str <- as.vector(strand(rds@last));
  stt <- start(grs) + g2v[, 3] - 1;
  end <- start(grs) + g2v[, 4] - 1;
  stt[str=='+'] <- end(grs[str=='+']) - g2v[str=='+', 3] + 1;
  end[str=='+'] <- end(grs[str=='+']) - g2v[str=='+', 4] + 1;
  g2v[, 3] <- pmin(stt, end);
  g2v[, 4] <- pmax(stt, end);
  
  ## Re-format and transform the read-to-vector table
  gtbl <- g2v;
  colnames(gtbl) <- paste0('G2V_', c('StartV', 'EndV', 'StartG', 'EndG', 'NMatch', 'NMismatch'));
  ##############################
  
  ######################################################################################################
  #### Prepare result tables
  tbl <- cbind(smm, rtbl, gtbl);
  spl <- split(tbl, tbl$Chromosome);
  chr <- suppressWarnings(as.integer(sub('^chr', '', names(spl))));
  spl <- spl[order(chr)]; 
  tbl <- do.call('rbind', spl);
  rownames(tbl) <- unlist(lapply(spl, rownames));
  
  out <- list(read.count=read.count, read.summary=tbl);
  
  invisible(out);
}

#### This function characterizes sequences at the genome-vector junction sites where integration happens, by sites
#### It uses output of CharacterizeGenomeVectorJunctionRead() to summarize reads mapped to the same junction sites
CharacterizeGenomeVectorJunctionSite <- function(smmr, gref, vref, pseq) {
  # smmr:   Read-level summary table as part of outputs from CharacterizeGenomeVectorJunctionRead
  # gref:   an BSgenome object of reference genome sequence
  # vref:   an DNAString object of vector sequence
  # pseq:   PCR primer sequence
  
  require(Biostrings);
  require(GenomicRanges);
  
  map2 <- GRanges(smmr$Chromosome, IRanges(smmr$StartR2, smmr$EndR2), smmr$StrandR2); # Genomic position Read 2 mapped to
  bjnc <- resize(map2, 1);                                                            # Genomic position next to the junction site
  pjnc <- paste0(smmr$Chromosome, ':', start(bjnc));                                  # Junction site ID
  site <- split(smmr, pjnc)[unique(pjnc)];                                            # Split reads by junction sites
  
  gr <- GRanges()
  
  ## Determine orientation of the junction sites
  npls <- sapply(site, function(s) nrow(s[s$StrandR2=='+', , drop=FALSE])); # No. of Read 2 mapped to + strand
  nmns <- sapply(site, function(s) nrow(s[s$StrandR2=='-', , drop=FALSE])); # No. of Read 2 mapped to - strand
  jstr <- rep('*', length(site));
  jstr[npls>nmns] <- '+'; # Genomic 
  jstr[npls<nmns] <- '-';
  
  ## Summary by sites
  gr  <- GRanges(sub(':[0-9]+$', '', names(site)), IRanges(as.integer(sub('^.+:', '', names(site))), width=1), strand=jstr);
  tbl <- data.frame(Chromosome=as.vector(seqnames(gr)), Strand=jstr, PositionG=start(gr), stringsAsFactors = FALSE);
  names(gr) <- rownames(tbl) <- names(site);
  tbl$NRead <- npls + nmns;
  tbl$NReadPlus <- npls;
  tbl$NReadMinus <- nmns; 
  tbl$StrandExcl <- as.integer(pmin(npls, nmns)==0)
  
  ## Map genomic sequences next to junction sites to PCR primer
  gseq <- ge
}
  