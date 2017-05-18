## SOFTWARE COPYRIGHT NOTICE AGREEMENT
## This software and its documentation are copyright (2013) by 
## Cancer Research UK (CRUK). This work is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License (http://creativecommons.org/licenses/by-nc-sa/3.0/). 
##
## This software is supplied as is without any warranty or guaranteed support
## whatsoever. CRUK can not be responsible for its use, misuse, or functionality.
# http://www.ploscollections.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1002598


### varSLR ###
# v1.0
# 28/07/2014
# Max Salm


#################
### Libraries ###
#################
if( TRUE ){
	## Change boolean if sourcing this file
	library(Rsamtools) # Going to have use the Bioconductor version for now
	# library(GenomicAlignments) # For readGAlignmentsFromBam(), if not found in Rsamtools
	library(GenomicRanges)	
	library(rtracklayer)
	library(compiler) # Byte Code compiler, to speed up functions
	library(BMA) # for bic.glm(), http://www.medicine.mcgill.ca/epidemiology/joseph/courses/EPIB-621/logselect.pdf
	# library(logistf) # for Firth logistic regression,  http://www.r-bloggers.com/example-8-15-firth-logistic-regression/	 
	library(brglm) ## Faster library for logistf
	library(BSgenome.Hsapiens.UCSC.hg19) ## Load reference genome
	# library(BSgenome.Hsapiens.UCSC.hg19.masked)
	# library(BSgenome.Hsapiens.UCSC.mm9)
	# library(BSgenome.Mmusculus.UCSC.mm9.masked)	

	SnpHsapiens <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20120608") ## Inject SNPs into reference genome (used in getBAMdata() but slow)	
}

#################
### Functions ###
#################

#########################
### Utility functions ###
#########################
prepareBlacklist <- function() {
	## Read and merge UCSC table data: genomicSuperDups, wgEncodeDacMapability & RepeatMasker. Save to an R object

	lcrs <- read.delim("/farm/home/salm01/NGS_gen/Ref_Data/hg19_genomicSuperDups_UCSC_12032012.bed",
						header=FALSE, skip=1, 
						colClasses=c("character", "integer", "integer", "NULL")) ## hg19 mapping low-copy repeats, downloaded from UCSC on 12/03/2012
	bla <- read.delim("/farm/home/salm01/NGS_gen/Ref_Data/wgEncodeDacMapabilityConsensusExcludable.bed",
						header=FALSE, 
						colClasses=c("character", "integer", "integer", "NULL", "NULL", "NULL"))
	blb <- read.delim("/farm/home/salm01/NGS_gen/Ref_Data/wgEncodeDukeMapabilityRegionsExcludable.bed",
					header=FALSE,
					colClasses=c("character", "integer", "integer", "NULL", "NULL", "NULL")) 					
	blacklist <- rbind(lcrs, bla, blb) ## merge "blacklisted" regions (ie. those that are highly repetitive)
	## Include RepeatMasked regions
	strs <- read.delim("/farm/home/salm01/NGS_gen/Ref_Data/hg19_RepeatMasker_UCSC_12032012.sorted.bed",
					header=FALSE, skip=1,
					colClasses=c("character", "integer", "integer", "NULL"))	
	blacklist <- rbind(blacklist, strs) 
	
	## Motifs that induce sequencing errors, http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3622629/
	# Skip, as mapped relative to GRCh37
	cse_GAII <- read.delim("/farm/home/salm01/NGS_gen/Ref_Data/hg19_RepeatMasker_UCSC_12032012.sorted.bed",
					header=FALSE, skip=1,
					colClasses=c("character", "integer", "integer", "NULL"))	

	
	
	blacklist_gr <- GRanges(seqnames = blacklist[, 1], ranges = IRanges(start=blacklist[, 2], end=blacklist[, 3])) ## convert to GRange object
	blacklist_hg19 <- reduce(blacklist_gr) # Merge overlapping intervals
	save(blacklist_hg19, file="blacklist_hg19.rda") # Save object to dir()
	tmp <- as.data.frame(blacklist_hg19)
	write.table(tmp[, 1:3], file="hg19_blacklist.bed", col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE)

}
# prepareBlacklist()

filterBlacklist <-function(x=snv, genome="hg19"){
	## Remove loci that map to ENCODE "blacklisted" regions (e.g. LCRs), alignments are unreliable in these regions a priori.
	## See, http://www.ncbi.nlm.nih.gov/pubmed/22178994
	# Args:
	#   x: A data.frame produced by variantConstructor() 
	#	genome: Genome build (currently only hg19 supported)
	# Returns:
	#   A data.frame with the "is_blacklist" column variable set.

	if (genome == "hg19") {
		load("/farm/home/salm01/projects/main_lib/SNV_QC/VarSLR/blacklist_hg19.rda") # see prepareBlacklist()
	}else{
		stop("No other genome build supported yet\n")
	}
	
	x_gr <- GRanges(seqnames = x$chromosome, 
			ranges = IRanges(start=x$position, end=x$position) ) ## convert x in GRange object
	ols <- as.matrix(suppressWarnings(findOverlaps(x_gr, blacklist_hg19))) ## find Overlaps
	cat(nrow(ols), "/", nrow(x), "loci in blacklisted regions.\n")
	
	if (length(ols) > 0) {	
		x$is_blacklist[ ols[, 1] ] <- TRUE
	}
	return(x)
}


filterMappability <-function(x=snv, genome_build="hg19"){
	## Remove loci that map to the following:
	##		simple repeats; 
	##		microsatellites;
	##		RepeatMasked regions;
	##		Segmental Duplications (aka LCRs);
	## 		ENCODE blacklisted regions;
	## alignments are unreliable in these regions a priori.
	## NB This requires network access, as it connects to UCSC.
	## See, http://www.ncbi.nlm.nih.gov/pubmed/22178994
	# Args:
	#   x: A data.frame produced by variantConstructor() 
	#	genome_build: Genome build, used for call to UCSC
	# Returns:
	#   A data.frame with the "is_blacklist" column variable set.
	
	mySession = browserSession("UCSC")
	genome(mySession) <- genome_build
	
	table_name <- "wgEncodeDacMapabilityConsensusExcludable"    
	
	x_gr <- GRanges(seqnames = x$chromosome, 
			ranges = IRanges(start=x$position, end=x$position) ) ## convert x in GRange object


	table_simple <- getTable(
					ucscTableQuery(mySession, track="Simple Repeats", range=reduce(x_gr), table="simpleRepeat")
						)

	table_microsat <- getTable(
					ucscTableQuery(mySession, track="Microsatellite", range=reduce(x_gr), table="microsat")
						)
						
	table_rmsk <- getTable(
					ucscTableQuery(mySession, track="RepeatMasker", range=reduce(x_gr), table="rmsk")
						)

	table_LCR <- getTable(
					ucscTableQuery(mySession, track="Segmental Dups", range=reduce(x_gr), table="genomicSuperDups")
						)
						
	table_DAC_Blacklisted_Regions <- getTable(
					ucscTableQuery(mySession, track="Mappability", range=reduce(x_gr), table="wgEncodeDacMapabilityConsensusExcludable")
						)
	
	table_Duke_Excluded_Regions <- getTable(
					ucscTableQuery(mySession, track="Mappability", range=reduce(x_gr), table="wgEncodeDukeMapabilityRegionsExcludable")
						)		
		
	targ_col <- c("chrom", "chromStart", "chromEnd")
	collated_tables <- rbind(as.matrix(table_simple[, targ_col]), 
							as.matrix(table_microsat[, targ_col]),
							as.matrix(table_rmsk[, c("genoName", "genoStart", "genoEnd")]), 
							as.matrix(table_LCR[, targ_col]), 
							as.matrix(table_DAC_Blacklisted_Regions[, targ_col]), 
							as.matrix(table_Duke_Excluded_Regions[, targ_col]))
	
	blacklist_hg19 <- GRanges(seqnames = collated_tables[, 1], 
		ranges = IRanges(start=as.numeric(collated_tables[, 2]), end=as.numeric(collated_tables[, 3]) ) ) ## convert into GRange object
		
	ols <- as.matrix(suppressWarnings(findOverlaps(x_gr, reduce(blacklist_hg19)))) ## find Overlaps
	
	cat(nrow(ols), "/", nrow(x), "loci in repetitive regions.\n")
	
	if (length(ols) > 0) {	
		x$is_blacklist[ ols[, 1] ] <- TRUE
	}
	return(x)
}				  

# test_data <- data.frame(
				# chromosome=rep("chr8", 6),
				# position=c(8996010, 8996032, 9000465, 9001392, 8964467, 8090413),
				# ref_allele=c("A", "C", "G", "C", "G", "G"),
				# var_allele=c("C", "G", "A", "T", "A", "A")
				# )				  		  
# test_data <- variantConstructor(test_data)	  
# tmp <- filterMappability(test_data)			  


weightedHomopolymerRate <- function(xx = "GGGTGCCCCCAAAATATT") {
	## http://www.broadinstitute.org/crd/wiki/index.php/Homopolymer
	## "The weighted homopolymer rate (WHR) of a sequence is a measure of the frequency of homopolymers in the sequence... The lowest possible WHR of a sequence is 1; the highest possible is the square of the sequence length (if N = 1). A randomly-generated sequence has an expected WHR of 20/9 ˜ 2.222. Most genomes have WHRs higher than the random value, due to imbalances in GC-content and the presence of junk DNA."
	xx_rle <- rle(unlist(strsplit(xx, "")))
	yy <- sum(xx_rle$lengths ^ 2) / length(xx_rle$lengths)
	return(yy)
}

# weightedHomopolymerRate("TGATTCAAGCATTCGATC") == 1.6
# weightedHomopolymerRate("GGGTGCCCCCAAAATATT") == 7.25
# weightedHomopolymerRate("AAAAAAAAAA") == 100
# weightedHomopolymerRate("AAAAATTTTT") == 25
				  
filterHomopolymer <-function(x=variants, flank_window=5, whr_threshold=10){
	## Remove loci that map to homopolymer runs, as variants are unreliable in these regions a priori, particularly indels.
	# Args:
	#   x: A data.frame produced by variantConstructor() 
	#	flank_window: The number of bases flanking the variant to screen
	#	whr_threshold: The weighted homopolymer rate threshold above which variants are filtered out (see http://www.broadinstitute.org/crd/wiki/index.php/Homopolymer)
	# Returns:
	#   A data.frame with the "is_blacklist" column variable set.
		
	flank_seq <- getSeq(Hsapiens, 
						names=x$chromosome, 
						start=x$position - flank_window, 
						end=x$position + flank_window, 
						as.character=FALSE)	

	whr <- unlist(lapply(X=as.character(flank_seq), FUN=weightedHomopolymerRate))
	ols <- which(whr >= whr_threshold)
	cat(length(ols), "/", nrow(x), "loci in homopolymer regions.\n")
	
	if (length(ols) > 0) {	
		x$is_blacklist[ ols ] <- TRUE
	}
	return(x)
}				  

ascii2base <- function(w="!") { 
	## Convert ASCII-encoded base (+33) to the phred-scaled base error probability.
	## char -> int
	## See, http://samtools.sourceforge.net/SAM1.pdf
	# Args:
	#   w: A single ASCII character 
	# Returns:
	#   Phred-scaled base error probability (i.e. -log10Pr(base is wrong)).
	if (length(w) == 1 && is.na(w)) {
		return(w)
	}else if(length(w) != 1) {
		stop("Only one element at a time please!")
	}else if(nchar(w) != 1){
		stop("Only one character at a time please!")
	}else{
		return( strtoi(charToRaw(w), 16L) - 33 )
	}
}

getCigar <- function(x="40M1I1P9M11=") {
	## Returns data.frame with listed cigar items.
	## string -> data.frame
	x <- as.character(x)
	counter <- unlist(strsplit(x, '[A-Z,=]', perl=TRUE))
	desc <- unlist(strsplit(x, '[0-9]', perl=TRUE))
	desc <- desc[which(desc != "")]
	out <- data.frame(Length=as.numeric(counter), Type=desc)
	test <- all(out$Type %in% c("M", "I", "N", "D", "S", "H", "P", "=", "X"))
	if( !(test) ) {
		stop("Unknown CIGAR Operation detected\n")
	}
	return(out)
}

QDiff <- function(x=test_ID$qual) {
	## Return the average inter-base quality scores across a complete read
	# Args:
	#   x: A 'qual' entry of a SAM
	# Returns:
	#   data.frame with various features derived from each read covering the specified locus
	tmp <- unlist(strsplit(as.character(x), ""))
	tmp1 <- unlist(lapply(X=tmp, FUN=ascii2base))
	#return(var(diff(tmp1))) ## Older version used inter-base quality variance, not mean.
	return( mean( abs( tmp1 ) ) ) ## Average base quality
}

getPositionInExpandedString <- function (x, views) {
   ## A hack of padAndClip()
	x_width <- width(x)
    x_names <- names(x)
    x_mcols <- mcols(x)

	views_start <- start(views)
	views_width <- width(views)	
		
    Lmargin <- 1L - views_start
    Rmargin <- views_width - x_width - Lmargin
    
	Lclipping <- pmin(pmax(-Lmargin, 0L), x_width)
    Rclipping <- pmin(pmax(-Rmargin, 0L), x_width)
	
    # clipped_x <- narrow(x, start = 1L + Lclipping, end = x_width - Rclipping, use.names = FALSE)

    # ans <- clipped_x
    # names(ans) <- x_names
    # mcols(ans) <- x_mcols
    # return(ans)
	
	return(x_width - Rclipping)	
}

getBAMdata <- function (file=LiveBamConnection, param=params){
 	## Gather raw data at a single base position from a BAM for the classification step
	## connection, ScanBamParam -> data.frame
	# Args:
	#   file: A BAM connection 
	#	param: ScanBamParam object passed to readGAlignmentsFromBam()
	
	## Code based on stackStringsFromBam()
	
    region <- unlist(bamWhich(param), use.names = FALSE) ## IRanges object (genomic interval)
	if( width(region) != 1 ){
		stop("Designed for a single nucleotide position only.\n")
	}
	
	gal <- readGAlignmentsFromBam(file, use.names = FALSE, param = param) ##  class 'GAlignments' 
    gal_mcols <- mcols(gal) # Get the metadata columns
	
	if(length(gal) == 0){stop("getBAMdata(): No alignments retrieved by readGAlignmentsFromBam()")}
	
	
	## 'seq' @ position
    seq_col <- gal_mcols[[ match('seq', colnames(gal_mcols)) ]] ## Subset to data of interest (ie seq/qual strings)
	
	## 'qual' @ position
    qual_col <- gal_mcols[[ match('qual', colnames(gal_mcols)) ]] ## Subset to data of interest (ie seq/qual strings)	
	qual_col <- BStringSet(qual_col) ## Convert to BStringSet datatype if 'qual'

	## Lay the sequences alongside the reference space, using their CIGARs (essentially applying CIGAR transformation of sequence).
    layed_seq <- sequenceLayer(seq_col, cigar(gal), D.letter = "-", N.letter = "-") 
    layed_qual <- sequenceLayer(qual_col, cigar(gal), D.letter = "!", N.letter = "!") 	
	
	## Stack sequences on the specied region.
    ans_seq <- stackStrings(x=layed_seq, from=start(region), to=end(region), shift = start(gal) - 1L,  		Lpadding.letter="+", Rpadding.letter="+") 
    ans_qual <- stackStrings(x=layed_qual, from=start(region), to=end(region), shift = start(gal) - 1L,  Lpadding.letter="!", Rpadding.letter="!") 	

	## Get genomic position index in read (from 5')
	## Calculate DETPE: "the distance to the effective 3' end (DETPE), which equals the smallest distance between the variant base and either the 3' end of the read, the 3' clipped end or the 5' end of the Q2 run (ie a run of '#'s in qual at 3' end ) normalized by the unclipped read length." http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268238/	
	genome_pos <- IRanges(start(region) - (start(gal) - 1L), end(region) - (start(gal) - 1L))
	string_pos <- getPositionInExpandedString(x=layed_seq, views=genome_pos)
	test <- read_idx <- detpe <- rep(NA, length(layed_seq))
	for(i in 1:length(layed_seq)) {
		test[i] <- as.character(layed_seq[[i]][ string_pos[i]  ]) # sanity check

		## TODO: Unit test: Is this implementation confounded by hard- or soft-clipping? H/S clips are removed from layed_seq, so I guess so...
		read_idx[i] <- sum(alphabetFrequency( layed_seq[[i]][ 1:string_pos[i]  ])[c('A', 'T', 'G', 'C')])
			
		## For DETPE...
		qual_rle <- rle(as.vector(as.matrix(layed_qual[i])))
		if ( rev(qual_rle$values)[1] == "#" ) {
			## "Beginning with version 1.5, the Illumina pipeline uses the Read Segment Quality Control Indicator to identify 3' portions of reads that should be discarded. This is done by setting the base quality for the region to 2."
			clip_from_3 <- rev(qual_rle$lengths)[1]
			detpe[i] <- abs((width(gal)[i] - clip_from_3) - read_idx[i]) / width(gal)[i]
		}else{
			detpe[i] <- abs(width(gal)[i] - read_idx[i]) / width(gal)[i]	
		}
	}
	if( !(all(as.matrix(ans_seq) == test))){stop("getBAMdata(): Indexing problem\n")} ## Sanity check
	
	## Adjust 'position in read' for strandness
	negs <- gal_mcols$strand == "-"
	read_idx[ negs ] <- qwidth(gal)[ negs ] - read_idx[ negs ] + 1
	detpe[ negs ]	<- 1 - detpe[ negs ]

	## Calculate MMQS (from layed_seq/layed_qual objects (taking SNPs into account))
	## "To identify reads that might support a paralog, we quantified the number of mismatches in a given read by summing the base qualities of the mismatches up for a given read the mismatch quality sum (MMQS)." http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268238/
	# For each read in layed_seq, retrieve corresponding reference genome window and create an basic object representing alignment	
	ref_seqs <- getSeq(SnpHsapiens, names=seqnames(gal), start=start(gal), end=end(gal), as.character=FALSE)
	mmqs <- function(xx=1) {
		## At non-dbSNP differences (and excluding query variant), return summed base qualities
		## NB. Ns in the read are included
		test_mat <- rbind( as.matrix(ref_seqs[xx]),
							as.matrix(layed_seq[xx]), 
							as.matrix(layed_qual[xx])
						)
		test_mat <- test_mat[, -string_pos[xx]]	## Remove query locus
		test_mat <- test_mat[, test_mat[1, ] %in% c("A", "T", "G", "C")] ## Remove SNPs
		diff_idx <- which( test_mat[1, ] != test_mat[2, ] )
		if (length(diff_idx) > 0) {
			return( sum( unlist(lapply(test_mat[3, diff_idx], FUN=ascii2base)) ) )	
		}else{
			return(0)
		}
	}
	mmqs_ls <- lapply(X=c(1:length(layed_seq)), FUN=mmqs)	

	## Prepare output data.frame
	qual_at_pos <- unlist(lapply(as.matrix(ans_qual), FUN=ascii2base))	
	qual_at_pos[ as.matrix(ans_seq) == '-' ] <- 0 ## Zero-out deletions
   	dfr <- data.frame(
				qname=gal_mcols$qname, 
				position=gal_mcols$pos, 
				revstrand=gal_mcols$strand == "-",
				seq=as.character(gal_mcols$seq),
				qual=as.character(gal_mcols$qual),
				cigar=as.character(gal_mcols$cigar),
				mapq=gal_mcols$mapq,
				B=as.matrix(ans_seq),
				Q=qual_at_pos,
				I=read_idx, 
				MMQS=unlist(mmqs_ls), 
				DETPE=detpe,
				stringsAsFactors=FALSE
				)	
   
	return(dfr)
}

getSNVData <- function(bam_idx_2=NA, locus_idx_3=NA, connection_list_3=NA) {
	## Gather feature data at a single base position for the classification step
	## connection, variant -> data.frame
	# Args:
	#   bam_idx_2: Index to a BAM connection [int]
	#	locus_idx_3: A row index of the object created by variantConstructor()
	
	# Returns:
	#   data.frame with various features derived from each read covering the specified locus

	locus <- VARIANTS[locus_idx_3, ]

	empty_output <- data.frame(Allele=NA,
								SampleIdx=NA,
								CyclePosition=NA,
								Qual=NA,
								Mapq=NA,
								Strand=NA,
								NumberOfCigarOps=NA,
								QDiff=NA,
								MMQS=NA,
								DETPE=NA)
	# Gatekeeper, ensuring ref/alt allele correctly specified
	test <-  all(c( locus$ref_allele, locus$var_allele) %in% c('A', 'T', 'G', 'C') )
	if( !test ) {
		warning(paste('REF/ALT alleles not {A, T, G, C} for SNV at', locus$chromosome, locus$position, "\n"))
		return(empty_output)
	}
	
	### Start...
	LiveBamConnection <- connection_list_3[[bam_idx_2]] # Get a BAM connection 
	
	# Use scanBamWhat() to get availbale feature vector, for "what" below.
	WHAT <- c("qname","rname", "strand", "pos", "mapq", "cigar", "seq", "qual")
	POS <- GRanges(seqnames = as.character(locus$chromosome), IRanges(start=as.numeric(locus$position), end=as.numeric(locus$position)))
	# TAG <- c("NM") # can be added, but no point for now
	# Would love to add MAPQ of paired read, but massively slows things down as have to use: 
	# mates <- findMateAlignment(readBamGappedAlignments(LiveBamConnection, use.names=TRUE, param=ScanBamParam(what=c("flag", "mrnm", "mpos")), which=POS))
	# Could use FLAG info instead (i.e. is pair mapped or not)
	params <- ScanBamParam(what=WHAT, 
							flag=scanBamFlag(isUnmappedQuery = FALSE),
							simpleCigar = FALSE, 
							reverseComplement = FALSE,
							which=POS)

	### Data retrievals from the BAM 	
	dfr <- getBAMdata(file=LiveBamConnection, param=params)
	
	### Remove overlapping reads (to avoid unfair double counting): not as useful as expected.
	# overlapping_read <- which(duplicated(dfr$qname)) 
	# if( length(overlapping_read) > 0 ) {
		# dfr <- dfr[-overlapping_read, ]
	# }
	
	## Next stage...
	if (nrow(dfr) == 0) {
		## If no data returned from call to BAM, return NAs
		warning(paste('No data retrieved', bam_idx_2, locus_idx_3 ))
		return(empty_output)	
	
	}else{		
		dfr <- na.omit(dfr) # Remove reads with missing data 	
		nCigar <- unlist(lapply(X=dfr$cigar, FUN=function(u){return(nrow(getCigar(u)))}))
		q_diff <- unlist(lapply(X=dfr$qual, FUN=QDiff)) # The mean read QUAL

		## Output
		out <- data.frame(Allele=dfr$B,
						SampleIdx=rep(bam_idx_2, nrow(dfr)),
						CyclePosition=dfr$I,
						Qual=dfr$Q,
						Mapq=dfr$mapq,
						Strand=dfr$revstrand,
						NumberOfCigarOps=nCigar,
						QDiff=q_diff,
						MMQS=dfr$MMQS,
						DETPE=dfr$DETPE
						)
		return(out)
	}
}

cigarIndel <- function(y="+TGCT") {
	## Translate an indel variant into a CIGAR representation
	## string -> string
	# Args:
	#   y: An indel variant
	# Returns:
	#   A CIGAR string 
	
	type <- substr(y, 1, 1)
	len <- nchar(substr(y, 2, nchar(y)))
	
	test_alphabet <- toupper(unique(unlist(strsplit(substr(y, 2, nchar(y)), ""))))
	if( any(!( test_alphabet %in% c("A", "T", "G", "C") )) ){
		stop(cat("cigarIndel: Incorrect indel allele coding, non-{A/T/G/C} character found, ", y, "\n"))
	}
	
	if (type == "+") {
		return(paste(len, "I", sep=""))
	}else if(type == "-"){
		return(paste(len, "D", sep=""))	
	}else{
		stop(cat("cigarIndel: Incorrect indel allele coding, unexpected symbol for indel, ", y, "\n{+/-} only please.\n"))
	}
}

posIdxIndel <- function(r=dfr[1, ], x) {
	## Retrieve base specific features at Indels
	## (matrix.row, int) -> list	
	# Args:
	#   r: An alignment object (see getINDELdata)
	#	x: A single entry of a "variant"-class object
	# Returns:
	#   A list, B=BASE and I=INDEX values for the specified position
	
	### Indel description/expansion
	ear <- expandAlign(r)				# Convert read into matrix
	var_i <- cigarIndel(x$var_allele)	#  Translate an indel variant into a CIGAR representation
	is_variant <- length(grep(var_i, r$cigar)) > 0  # Is the CIGAR representation present in CIGAR string?
	is_insertion <- length(grep("I", var_i)) > 0 	# Is it an insertion?
	is_deletion <- length(grep("D", var_i)) > 0		# Is it a deletion?
	five_prime <- which(ear$idx_s == x$position) + 1 # Where are we looking genomically?
	b <- "0" # Lets start with the read representing reference variant...update expected
	qi <- NA # Initialise qi
	
	if (is_variant && is_insertion && length(five_prime) != 0 && any(ear$idx_q == five_prime,  na.rm = TRUE) ) {
		## Indel of correct size detected in CIGAR...
		## Does the insertion match nucleotides at correct position?
		nucs <- substr(x$var_allele, 2, nchar(x$var_allele))
		if (five_prime == max(ear$idx_q, na.rm=TRUE) ) {
			three_prime <- five_prime # If it's the last read position...
		}else{
			three_prime <- which(ear$idx_s == x$position+1) - 1
			if (length(three_prime) == 0) {three_prime <- five_prime}		
		}
		insert <- paste(ear$B[five_prime:three_prime], collapse="")
		if(nucs == insert) {
			b <- "1"
			# Mean Quality of bases in insertion
			qa <- as.character(ear$Q[five_prime:three_prime])
			qb <- unlist(lapply(X=qa, FUN=ascii2base))
			qi <- round(mean(qb, na.rm=TRUE))	
		}
		if (r$revstrand) {
			five_prime <- nchar(as.character(r$seq)) - five_prime + 1 # If on reverse strand, flip index...
		}		
		return(list(B=b, Q=qi, I=five_prime))
	}else if (is_variant && is_deletion && length(five_prime) != 0 && any(ear$idx_q == five_prime,  na.rm = TRUE) ) {
		## Does the deletion position match expected position?
		## Can't check nucleotides, as they are not expected to be in read.
		three_prime <- five_prime + nchar(x$var_allele) - 2
		del_i <- is.na(ear$B[five_prime:three_prime]) # Is the deletion interval missing bases?
		del_i <- c(del_i, !( is.na( ear$B[five_prime-1] ) ) ) # Are the flanking borders not missing bases?
		del_i <- c(del_i, !( is.na( ear$B[three_prime + 1]) ) )
		if( all(del_i) ) {
			b <- "1" 				# The read contains the deletion
			# Mean Quality @ deletion flanking bases
			qa <- as.character(ear$Q[five_prime-1])
			qb <- as.character(ear$Q[three_prime + 1])
			qi <- round(mean(c(ascii2base(qa), ascii2base(qb)), na.rm=TRUE))
		}
		if (r$revstrand) {
			five_prime <- nchar(as.character(r$seq)) - five_prime + 1 # If on reverse strand, flip index...
		}		
		return(list(B=b, Q=qi, I=five_prime))
	}else{
		## Indel of correct size not detected in CIGAR	
		qi = as.character(ear[five_prime, "Q"]) # Quality @ position i
		if (r$revstrand) {
			five_prime <- nchar(as.character(r$seq)) - five_prime + 1 # If on reverse strand, flip index...
		}		
		if (length(five_prime) == 0) {b<-NA; qi<-NA; five_prime <- NA} # Position not covered.
		return(list(B=b, Q=ascii2base(qi), I=five_prime)) 
	}
}

getBAMdataIndel <- function(file=LiveBamConnection, param=params, is_deletion=TRUE, locus=NA){
 	## Gather raw data at a deletion (or insertion) position from a BAM for the classification step
	## connection, ScanBamParam -> data.frame
	# Args:
	#   file: A BAM connection 
	#	param: ScanBamParam object passed to readGAlignmentsFromBam()
	
	## Code based on stackStringsFromBam()
	
	# ptm <- proc.time() ## Start a timer
	
    region <- unlist(bamWhich(param), use.names = FALSE) ## IRanges object (genomic interval)
	
	gal <- readGAlignmentsFromBam(file, use.names = FALSE, param = param) ##  class 'GAlignments' 
    gal_mcols <- mcols(gal) # Get the metadata columns
	
	if(length(gal) == 0){stop("getBAMdata(): No alignments retrieved by readGAlignmentsFromBam()")}
	
	
	## 'seq' @ position
    seq_col <- gal_mcols[[ match('seq', colnames(gal_mcols)) ]] ## Subset to data of interest (ie seq/qual strings)
	
	## 'qual' @ position
    qual_col <- gal_mcols[[ match('qual', colnames(gal_mcols)) ]] ## Subset to data of interest (ie seq/qual strings)	
	qual_col <- BStringSet(qual_col) ## Convert to BStringSet datatype if 'qual'

	## Lay the sequences alongside the reference space, using their CIGARs (essentially applying CIGAR transformation of sequence). This will remove the substrings associated with insertions to the reference (I operations) and soft clipping (S operations), and will inject new substrings (filled with "-") where deletions from the reference (D operations) and skipped regions from the reference (N operations) occurred during the alignment process:
    layed_seq <- sequenceLayer(seq_col, cigar(gal), from="query", to="reference", D.letter = "-", N.letter = "-") 
    layed_qual <- sequenceLayer(qual_col, cigar(gal), from="query", to="reference", D.letter = "!", N.letter = "!") 	
	
	# cat("[1]", proc.time() - ptm, "\n"); ptm <- proc.time();
	
	## Stack sequences on the specied region.
	if( is_deletion ) {
		ans_seq <- stackStrings(x=layed_seq, from=start(region), to=end(region) + 2, shift = start(gal) - 1L,  		Lpadding.letter="+", Rpadding.letter="+") ## Extra 2-bp padding, leads to one base flanking each side of deletion
		ans_qual <- stackStrings(x=layed_qual, from=start(region), to=end(region) + 2, shift = start(gal) - 1L,  Lpadding.letter="!", Rpadding.letter="!") 	
	}else{
		ans_seq <- stackStrings(x=layed_seq, from=start(region), to=end(region), shift = start(gal) - 1L,  		Lpadding.letter="+", Rpadding.letter="+")
		ans_qual <- stackStrings(x=layed_qual, from=start(region), to=end(region), shift = start(gal) - 1L,  Lpadding.letter="!", Rpadding.letter="!") 		
	}
	
	## Get Base/Quality Info
	if( is_deletion ) {
		## For DEL...
		base_indel <- as.matrix(ans_seq)
		base_indel <- base_indel[, -c(1, ncol(base_indel)), drop=FALSE] ## Trim off flanking bases
		base_indel <- apply(X=base_indel, MARGIN=1, FUN=paste, collapse="")	

		qualIndel <- function(xx=1) {
			tmp <- unlist(lapply(X=unlist(strsplit(as.character(ans_qual[xx]), "")), FUN=ascii2base))
			return(mean(tmp[tmp != 0])) ## Ignore zero-entries for DEL
		}
		qual_ls <- lapply(X=c(1:length(layed_seq)), FUN=qualIndel)
	
	}else{
		## For INS...
		base_indel <- as.matrix(ans_seq)
		base_indel <- apply(X=base_indel, MARGIN=1, FUN=paste, collapse="")	
		qual_tmp <- as.character(ans_qual)
		
		bait_cigar <- paste((nchar(locus$var_allele) - 1), "I", sep="") ## Which reads carry an insertion of correct size? 
		bait_hooked <- grep(bait_cigar, gal_mcols$cigar) ## Find candidate insertions in CIGAR
		
		for (n in bait_hooked) {
			cigar_ops <- cigarRangesAlongQuerySpace(gal_mcols$cigar[n], with.ops=TRUE)
			targ <- paste(width(cigar_ops[[1]]), names(cigar_ops[[1]]), sep="") == bait_cigar
			if( any(targ) ) {
				tmp <- cigar_ops[[1]][targ]
				potential_nucmatch <- rep("N", length(tmp))
				for(j in 1:length(tmp)){
					## To handle the case when there is more than one candidate insertion...
					potential_nucmatch[j] <- as.character( narrow(seq_col[n], start=start(tmp[j]), end=end(tmp[j]))) ## Retrieve/update nucleotide string in insertion
				}
				nuc_targ <- which(potential_nucmatch == substring(locus$var_allele, 2))
				if(length(nuc_targ) > 0) {
					## Performe the update if an insertion nucloetide match was found
					tmp <- tmp[nuc_targ[1]]
					base_indel[n] <- as.character( narrow(seq_col[n], start=start(tmp), end=end(tmp))) ## Retrieve/update nucleotide string in insertion (only works one range at a time)
					qual_tmp[n] <- as.character( narrow(qual_col[n], start=start(tmp), end=end(tmp)) ) ## Retrieve/update quality string in insertion (only works one range at a time)
				}
			}
		}
		
		qualIndel <- function(xx=1) {
			tmp <- unlist(lapply(X=unlist(strsplit(as.character(qual_tmp[xx]), "")), FUN=ascii2base))
			return(mean(tmp)) 
		}
		qual_ls <- lapply(X=c(1:length(qual_tmp)), FUN=qualIndel)	
	}
	
	# cat("[2]", proc.time() - ptm, "\n"); ptm <- proc.time();
	
	## Get genomic position index in read (from 5', and start of deletion)
	## Calculate DETPE: "the distance to the effective 3' end (DETPE), which equals the smallest distance between the variant base and either the 3' end of the read, the 3' clipped end or the 5' end of the Q2 run (ie a run of '#'s in qual at 3' end ) normalized by the unclipped read length." http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268238/	
	genome_pos <- IRanges(start(region) - (start(gal) - 1L), end(region) - (start(gal) - 1L))
	string_pos <- getPositionInExpandedString(x=layed_seq, views=genome_pos)
	read_idx <- detpe <- rep(NA, length(layed_seq))
	
	for(i in 1:length(layed_seq)) {
		## TODO: Unit test: Is this implementation confounded by hard- or soft-clipping? H/S clips are removed from layed_seq, so I guess so...could also use 'query-before-hard-clipping' in sequenceLayer()
		read_idx[i] <- sum(alphabetFrequency( layed_seq[[i]][ 1:string_pos[i]  ])[c('A', 'T', 'G', 'C')])
			
		## For DETPE...
		qual_rle <- rle(as.vector(as.matrix(layed_qual[i])))
		if ( rev(qual_rle$values)[1] == "#" ) {
			## "Beginning with version 1.5, the Illumina pipeline uses the Read Segment Quality Control Indicator to identify 3' portions of reads that should be discarded. This is done by setting the base quality for the region to 2."
			clip_from_3 <- rev(qual_rle$lengths)[1]
			detpe[i] <- abs((width(gal)[i] - clip_from_3) - read_idx[i]) / width(gal)[i]
		}else{
			detpe[i] <- abs(width(gal)[i] - read_idx[i]) / width(gal)[i]	
		}
	}
	
	## Adjust 'position in read' for strandness
	negs <- gal_mcols$strand == "-"
	read_idx[ negs ] <- qwidth(gal)[ negs ] - read_idx[ negs ] + 1
	detpe[ negs ]	<- 1 - detpe[ negs ]

	## Calculate MMQS (from layed_seq/layed_qual objects (taking SNPs into account))
	## "To identify reads that might support a paralog, we quantified the number of mismatches in a given read by summing the base qualities of the mismatches up for a given read the mismatch quality sum (MMQS)." http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3268238/
	# For each read in layed_seq, retrieve corresponding reference genome window (with dbSNP alleles having been inject into it) and create an basic object representing alignment	
	
	
	# cat("[3]", proc.time() - ptm, "\n"); ptm <- proc.time();
	ref_seqs <- getSeq(SnpHsapiens, names=seqnames(gal), start=start(gal), end=end(gal), as.character=FALSE) ## SLOW, takes ~8s for 150 reads - weirdly not much faster if retirive whole window once (to then use coordinates to subset). Perhaps I have to load into memory, rather than relying on disk access?

	# ref_seqs <- getSeq(SnpHsapiens, names=seqnames(gal)[1], start=min(start(gal)), end=max(end(gal)), as.character=FALSE) ## Get whole window once and use coordinates to subset
	
	# cat("[4]", proc.time() - ptm, "\n"); ptm <- proc.time();

	mmqs <- function(xx=1) {
		## At non-dbSNP differences (and excluding query variant), return summed base qualities
		## NB. Ns in the read are included
		
		test_mat <- rbind( 
						as.matrix(c(ref_seqs[xx], layed_seq[xx])), 
						as.matrix(layed_qual[xx]) 
						) ## For INS, target insertion implicitly excluded from object, so no adjustment required
						
		test_mat <- test_mat[, test_mat[1, ] %in% c("A", "T", "G", "C")] ## Remove SNPs
		diff_idx <- which( test_mat[1, ] != test_mat[2, ] )
		if (length(diff_idx) > 0) {
			return( sum( unlist(lapply(test_mat[3, diff_idx], FUN=ascii2base)) ) )	## DEL doesn't add to MMQS due to padding character chosen, so no adjustment required. 
		}else{
			return(0)
		}
	}
	mmqs_ls <- lapply(X=c(1:length(layed_seq)), FUN=mmqs) ## SLOW, takes ~10s for ~150 reads	

	# cat("[5]", proc.time() - ptm, "\n"); ptm <- proc.time();
	
	## Prepare output data.frame
	dfr <- data.frame(
				qname=gal_mcols$qname, 
				position=gal_mcols$pos, 
				revstrand=gal_mcols$strand == "-",
				seq=as.character(gal_mcols$seq),
				qual=as.character(gal_mcols$qual),
				cigar=as.character(gal_mcols$cigar),
				mapq=gal_mcols$mapq,
				B=base_indel,
				Q=unlist(qual_ls),
				I=read_idx, 
				MMQS=unlist(mmqs_ls), 
				DETPE=detpe,
				stringsAsFactors=FALSE
				)	
   ## Remove entries in the "padded zone"
   rem <- grep("\\+", dfr$B)
   if(length(rem) > 0) {
		dfr <- dfr[-rem, ]
   }
   
	return(dfr)
}

getINDELdata <- function(bam_idx_2=NA, locus_idx_3=NA, connection_list_3=NA) {
	## Gather feature data at a base position for the classification step	
	## connection, variant -> data.frame
	# Args:
	#   bam_idx_2: Index to a BAM connection [int]
	#	locus_idx_3: A row index of the object created by variantConstructor()
	#	connection_list_3: 
	# Returns:
	#   data.frame with various features derived from each read covering the specified locus
	### From VarScan2 site
	# A Note on Indel Positions
	#"VarScan calls consensus bases, SNPs, and indels at the position reported by SAMtools in the pileup file. For SNPs and consensus bases, this is the 1-based position of the site or variant. INDELS, HOWEVER, ARE REPORTED AT THE BASE IMMEDIATELY UPSTREAM OF WHERE THEY OCCUR. THUS, THE FIRST INSERTED BASE OCCURS AT (POSITION + 1) AND THE FIRST DELETED BASE OCCURS AT (POSITION + 1). THIS IS WHY, FOR DELETIONS, THE REFERENCE BASE (E.G. A) MAY NOT MATCH THE DELETED SEQUENCE (E.G. */C). "

	locus <- VARIANTS[locus_idx_3, ]

	empty_output <- data.frame(Allele=NA,
								SampleIdx=NA,
								CyclePosition=NA,
								Qual=NA,
								Mapq=NA,
								Strand=NA,
								NumberOfCigarOps=NA,
								QDiff=NA)

	
	### Start...
	LiveBamConnection <- connection_list_3[[bam_idx_2]] # Get a BAM connection 
	
	# Use scanBamWhat() to get availbale feature vector, for "what" below.
	WHAT <- c("qname","rname", "strand", "pos", "mapq", "cigar", "seq", "qual")
	POS <- GRanges(seqnames = as.character(locus$chromosome), IRanges(start=as.numeric(locus$position), end=as.numeric(locus$position)+1))	
	TAG <- c("NM", "MD") # can be added, but what do we do when tag is not present? Sorted.
	# TAG = character(0)
	params <- ScanBamParam(what=WHAT, 
							flag=scanBamFlag(isUnmappedQuery = FALSE),
							simpleCigar = FALSE, 
							reverseComplement = FALSE,
							tag = TAG,
							which=POS)

	## retrievals from the BAM 
	# RawData <- scanBam(LiveBamConnection, param=params) # The scanBam,character-method returns a list of lists. The outer list groups results from each Ranges list of bamWhich(param); the outer list is of length one when bamWhich(param) has length 0. Each inner list contains elements named after scanBamWhat(); elements omitted from bamWhat(param) are removed.
	# dfr <- data.frame(qname=RawData[[1]]$qname, 
				# position=RawData[[1]]$pos, 
				# revstrand=RawData[[1]]$strand == "-",
				# seq=as.character(RawData[[1]]$seq),
				# qual=as.character(RawData[[1]]$qual),
				# cigar=as.character(RawData[[1]]$cigar),
				# mapq=RawData[[1]]$mapq,
				# stringsAsFactors=FALSE)	
	# if( !is.null(RawData[[1]]$tag$NM) ){
		# dfr['editDist'] <- RawData[[1]]$tag$NM
	# }
	
	indel_type <- substr(locus$var_allele, 1, 1)
	if ( indel_type == "-" ) {
		dfr <- getBAMdataIndel(file=LiveBamConnection, param=params, is_deletion=TRUE, locus=locus)
	}else if ( indel_type == "+" ){
		dfr <- getBAMdataIndel(file=LiveBamConnection, param=params, is_deletion=FALSE, locus=locus)
	}else{
		stop(cat("getINDELdata(): Indel allele encoding not recognised for", paste(locus[1:4], collapse="|"), "\n"))
	}
	
	#########################################################################################################
	
	## Remove overlapping reads (to avoid unfair double counting)
	# overlapping_read <- which(duplicated(dfr$qname)) 
	# if( length(overlapping_read) > 0 ) {
		# dfr <- dfr[-overlapping_read, ]
	# }
						
	if (nrow(dfr) == 0) {
		## If no data returned from call to BAM, return NAs
		warning(paste('No data retrieved', bam_idx_2, locus_idx_3 ))
		return(empty_output)	
		
	}else{					
		dfr <- na.omit(dfr) # Remove reads with missing data	

		#dfr$position <- dfr$position + 1 # "The Fencepost Problem": express genome position in 1-based encoding. In this way, it matches BLAT alignments at UCSC (easier for debugging indels)
		# Add a few extra features to above data.frame (dfr)
		
		nCigar <- unlist(lapply(X=dfr$cigar, FUN=function(u){return(nrow(getCigar(u)))}))

		q_diff <- unlist(lapply(X=dfr$qual, FUN=QDiff)) # The mean read QUAL
		# q_diff <- rep(NA, nrow(dfr))
		
		## Output
		out <- data.frame(Allele=dfr$B,
						SampleIdx=rep(bam_idx_2, nrow(dfr)),
						CyclePosition=dfr$I,
						Qual=dfr$Q,
						Mapq=dfr$mapq,
						Strand=dfr$revstrand,
						NumberOfCigarOps=nCigar,
						QDiff=q_diff,
						MMQS=dfr$MMQS,
						DETPE=dfr$DETPE
						)		

		### For indels, mapq is a confounder rather than a useful feature, as the mismatch by the indel causes a bias proportional to indel length. Therefore, will use edit distance (adjusted for presence/absence of target indel) instead. This can only be done if edit distance is included the BAM.
		# if( !is.null(RawData[[1]]$tag$NM) ){
			# tmp_dist <- dfr$editDist - (as.numeric(as.character((out$Allele))) * (nchar(x$var_allele)-1))
			# out$Mapq <- tmp_dist
		# }
		return(out)
	}
}

gatherFeatures <- function(bam_idx_1=NA, locus_idx_2=NA, connection_list_2=NA) {
	## Simple wrapper for getSNVData & getINDELdata
	# Args:
	#   bam_idx_1: Index to a BAM connection [int]
	#	locus_idx_2: A row index of the object created by variantConstructor()
	#	bam_set: Output from makeBamSet()	
	# Returns:
	#   data.frame with various features derived from each read covering the specified locus
	
	if( VARIANTS[locus_idx_2, 'is_SNV'] ){
		## SNV section
		out <- getSNVData(bam_idx_2=bam_idx_1, locus_idx_3=locus_idx_2, connection_list_3=connection_list_2)
		return(out)
		
	}else if( !(VARIANTS[locus_idx_2, 'is_SNV']) ){
		## Indel section
		out <- getINDELdata(bam_idx_2=bam_idx_1, locus_idx_3=locus_idx_2, connection_list_3=connection_list_2)
		return(out)
		
	}else{
		stop("Some Error in gatherFeatures()")
		
	}
}


########################
### Public functions ###
########################

assign2run <- function(germline_bam, tumour_bams, variants) {
	## Convenience wrapper to assign two objects (BAMSET/VARIANTS) to global environment for varSLR run 
	# (string, vector, data.frame) -> NULL
	# Args:
	#	germline_bam: Full path to germline BAM [string]
	#   tumour_bams: Full path to tumour BAM file(s) [vector of strings]
	#	variants: A data.frame prepared by variantConstructor()
	# Returns:
	#   Nothing


	if(length(germline_bam) != 1){
		stop('assign2run(): Please specify a single reference BAM file\n')
	}
	if(length(germline_bam) < 1) {
		stop('assign2run(): Please specify at least one tumour BAM file\n')
	}
	
	if(any(ls(envir=.GlobalEnv, pattern="BAMSET") %in% "BAMSET")){
		warning("assign2run(): An object named BAMSET already existed in GlobalEnv and has been overwritten!\n")
		rm("BAMSET", envir=.GlobalEnv) ## Remove existing target objects
	}

	if(any(ls(envir=.GlobalEnv, pattern="VARIANTS") %in% "VARIANTS")){
		warning("assign2run(): An object named VARIANTS already existed in GlobalEnv and has been overwritten!\n")
		rm("VARIANTS", envir=.GlobalEnv) ## Remove existing target objects
	}
	
	SnpHsapiens <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20120608") ## Inject SNPs into reference genome (used in getBAMdata() but slow)	

	bams <- c(germline_bam, tumour_bams) ## This ensures that 'germline' is the baseline sample for logistic regression later (ie always sample == 1)
	## Test if BAM/BAI exist
	missing_bam <- !(file.exists(bams))
	if( any( missing_bam ) ){
		cat("Missing file:", bams[missing_bam], "\n")
		stop("assign2run(): BAM not found\n")
	}
	
	bais <- paste(bams, ".bai", sep="")
	missing_bai <- !(file.exists( bais ))
	if( any( missing_bai ) ){
		cat("Missing file:", bais[missing_bai], "\n")
		stop("assign2run(): BAI not found\n")
	}
	
	## Return list (and assign representation to global env as well)
	BAMSET <- list( bams=bams,
				 bai=bais)
	assign("BAMSET", BAMSET, envir=globalenv())
	
	## Assign representation to global env as well
	VARIANTS <- variants
	assign("VARIANTS", VARIANTS, envir=globalenv())

}

variantConstructor <- function(x=snv) {
	## A data.frame constructor to organise the variant data
	## x: a matrix 'x'
	## where the first four columns contain: chromosome, genomic position, reference allele, variant allele.
	
	if(class(x) == "data.frame") {
		## Convert data.frame to matrix
		x <- as.matrix(x)
	}
	
	if (class(x) == "matrix") {
		if(length(grep("chr", x[, 1])) > 0) {
			cat("Assuming UCSC chromosome naming convention\n")
		}else{
			cat("Assuming NCBI chromosome naming convention\n")
		}
		
		if(any(is.na(as.numeric(x[, 2])))) {
			stop("variantConstructor(): Non-numeric values in 'position' column")
		}
		x <- x[order(x[, 1], as.numeric(x[, 2])), ]## Order data

		z <- data.frame(chromosome=as.character(x[, 1]),
				position=as.numeric(x[, 2]),
				ref_allele=toupper(as.character(x[, 3])),
				var_allele=toupper(as.character(x[, 4])), 
				stringsAsFactors=FALSE)
		z['is_SNV'] <- z[, 'ref_allele'] %in% c("A", "T", "G", "C") & z[, 'var_allele'] %in% c("A", "T", "G", "C")
		z['is_blacklist'] <- rep(FALSE, nrow(z))
		
		
		## Test data validity [position data]
		if( class( z$position ) != "numeric" ){stop("variantConstructor(): Wrong data-type for positions (ie non-numeric).\n")}
		test <- range(z$position)
		if(test[1] < 0){stop("variantConstructor: Negative genome positions found.\n")}
		if(test[2] > 500000000){warning("variantConstructor: Unrealistically large genome position found (>500M). I've been trained on human chromosomes, where unless things have changed a lot, the largest chromosome was ~250M base-pairs long.\n")}
		
		## Test data validity [Reference alleles - SNV]			
		singleton <- unique(z$ref_allele[z$is_SNV])
		unexpected_nucleotide <- which(!(singleton %in% c("A", "T", "G", "C")))
		if(length(unexpected_nucleotide) > 0) {
			warning(cat("variantConstructor: SNV reference alleles must be A/T/G/C, but I found {", singleton[unexpected_nucleotide], "}. Please either include the correct reference allele [A/T/G/C] or delete query\n"))
			targ_black <- which(z$ref_allele %in% singleton[unexpected_nucleotide])    ## Set 'blacklist' flag
			z$is_blacklist[targ_black] <- TRUE
		}
		
		## Test data validity [alternative alleles - SNV]			
		singleton <- unique(z$var_allele[z$is_SNV])
		unexpected_nucleotide <- which(!(singleton %in% c("A", "T", "G", "C")))
		if(length(unexpected_nucleotide) > 0) {
			warning(cat("variantConstructor: SNV alternative alleles must be A/T/G/C, but I found {", singleton[unexpected_nucleotide], "}. Please either include the correct alternative allele [A/T/G/C] or delete query\n"))
			targ_black <- which(z$var_allele %in% singleton[unexpected_nucleotide])    ## Set 'blacklist' flag
			z$is_blacklist[targ_black] <- TRUE
		}		
		
		## Test data validity [Reference alleles - Indel]			
		nonsingleton <- unique(z$ref_allele[!(z$is_SNV)])
		unexpected_nucleotide <- which(!(nonsingleton %in% c("A", "T", "G", "C")))
		if(length(unexpected_nucleotide) > 0) {
			warning(cat("variantConstructor: INDEL reference alleles must be A/T/G/C. Please either include the correct reference allele or delete query.\n"))
			targ_black <- which(z$ref_allele %in% nonsingleton[unexpected_nucleotide])    ## Set 'blacklist' flag
			z$is_blacklist[targ_black] <- TRUE	
		}
		
		## Test data validity [Alternative alleles - Indel]			
		nonsingleton <- unique(z$var_allele[!(z$is_SNV)])		
		wrong_encoding <- which(nchar(nonsingleton) <= 1)
		if(length(wrong_encoding) > 0) {
			warning(cat("variantConstructor: That is a mighty short indel! Variant alleles are expected to be greater than 1-bp in length. \nProblematic alleles:", nonsingleton[wrong_encoding], "\n"))	
			targ_black <- which(z$var_allele %in% nonsingleton[wrong_encoding])    ## Set 'blacklist' flag
			z$is_blacklist[targ_black] <- TRUE
		}
		plus_minus <- unique(unlist(lapply(nonsingleton, FUN=substr, 1, 1)))
		wrong_encoding <- which(!(plus_minus %in% c("+", "-")))
		if(length(wrong_encoding) > 0) {
			warning(cat("variantConstructor: Uh-oh. Indel variant alleles follow the VarScan2 output format, and thus expect to be prepended with either '+' (insertion) or '-' (deletion). This was not uniformly the case. \nProblematic characters found:", plus_minus[wrong_encoding], "\n"))	
		}		
		tmp <- unique(unlist(strsplit(nonsingleton, "")))
		tmp <- tmp[!(tmp %in% c("+", "-"))]
		unexpected_nucleotide <- which(!(tmp %in% c("A", "T", "G", "C")))
		if(length(unexpected_nucleotide) > 0) {
			warning(cat("variantConstructor: All nucleotides in INDEL variant alleles must be A/T/G/C. \nProblematic characters found:", tmp[unexpected_nucleotide], "\n"))
		}
		
		## Test data validity [Does Reference allele != Alternative allele]	
		z$is_blacklist[z$ref_allele == z$var_allele] <- TRUE
		
		## Return object 
		VARIANTS <- z
		return(VARIANTS)
	
	}else{
		stop("variantConstructor(): expect x to be of class matrix")
	}
}

gatherFeaureSet <- function(locus_idx_1=x, d=10000, min_allele_count=4) {
	## Apply gatherFeatures to a set of BAM indexes, preparing output for logistic regression
	## data.frame -> data.frame 
	# Args:
	#   locus: A index (row) of the object created by variantConstructor()
	#	bam_set: Output from makeBamSet()
	#	d: Maximum average read depth
	# 	min_allele_count: Minimum number of variant alleles required for whole dataset
	# Returns:
	#   list of data.frame with various features derived from each read covering the specified locus, and a string	

	locus <- VARIANTS[locus_idx_1, ]
	
	tmp <- data.frame(Allele=as.character(NA),
				SampleIdx=as.integer(NA),
				CyclePosition=as.integer(NA),
				Qual=as.integer(NA),
				Mapq=as.integer(NA),
				Strand=as.logical(NA),     
				NumberOfCigarOps=as.integer(NA), 
				QDiff=as.numeric(NA))
	output_msg <- NA
	
	## Gatekeeper 1: Is the locus 'blacklisted'	
	if( locus$is_blacklist ){
		output_msg <- 'Skip:Blacklisted'
		return( list(output_df=tmp, output_msg=output_msg) ) ## Blacklisted locus, so returning NA
	}
	
	## Due to the way Rsamtools is set-up, we MUST initialise seperate connections to the BAM file set for each locus, because when run in parallel traffic conflicts on the same connection appear to occur.
	## NB When 'yieldSize' is set, scanBam() will iterate through the file in chunks, but is a pain to implement
	connection_list_1 <- BamFileList( BAMSET$bams, index=BAMSET$bai )
	open(connection_list_1)
	
	## Gatekeeper 2: Check for 'excessive' mean coverage
	POS <- GRanges(seqnames = as.character(locus$chromosome), IRanges(start=locus$position, end=locus$position))
	depth_at_position <- countBam(connection_list_1, 
								param=ScanBamParam(what='pos', simpleCigar = FALSE, 
												reverseComplement = FALSE, which=POS))
	if(mean(depth_at_position$records) > d) {
		output_msg <- paste('Skip:Excessive mean read depth=', round(mean(depth_at_position$records)), sep="")
		return( list(output_df=tmp, output_msg=output_msg) ) # Mean read number at position exceeds threshold 'd', return NAs	
	}

	## Start...
	gathered_data <- new_list <- vector("list", length=length(BAMSET$bams))
	for(bam_idx_0 in 1:length(BAMSET$bams)){
		depth_test <- depth_at_position[BAMSET$bams[ bam_idx_0 ], "records"] > 1 ## Has to be at least 2 reads covering position in BAM...
		if(depth_test){
			gathered_data[[bam_idx_0]] <- gatherFeatures(bam_idx_1=bam_idx_0, locus_idx_2=locus_idx_1, connection_list_2=connection_list_1)
		}
	}
	
	## Close connections to BAMs - VITAL FOR PARALLELISM TO WORK
	close( connection_list_1 )

	## If data gathering failed for single BAM, abort (except if there is insufficient coverage in BAM). Otherwise, convert to df.
	if( any( unlist(lapply(X=gathered_data, FUN=nrow)) == 0 ) ) {
		output_msg <- 'Skip:Data missing from some BAMs'
		return( list(output_df=tmp, output_msg=output_msg) ) # Mean read number at position exceeds threshold 'd', return NAs	
	}else{
	  gathered_data <- do.call("rbind", gathered_data)
	}
	
	## Gatekeeper 3: For SNVs, ensure < 50% of overlapping reads don't contain a deletion (ie skip compound variants) 	
	if( locus$is_SNV ){
		if( !all(is.na( gathered_data )) ) {
			del_af <- (length(which(is.na(gathered_data$Allele))) / nrow(gathered_data))
			if(del_af >= 0.5) {
				output_msg <- 'Skip:Overlapping deletion'
				return( list(output_df=tmp, output_msg=output_msg) ) 
			}
		}
	}
	
	## QC of gathered data, and return()
	if ( all(is.na( gathered_data )) ) {
		return( list(output_df=gathered_data, output_msg='Skip:No data')  )	## No data retrieved, so returning NA
	}else{
		## code ref/var allele as 0/1, and blank out non-expected calls
		if ( locus$is_SNV ) {
			coded <- as.character(gathered_data$Allele)
			coded[ which(coded != locus$ref_allele & coded != locus$var_allele) ] <- NA
			coded[which(coded == locus$var_allele)] <- 1
			coded[which(coded == locus$ref_allele)] <- 0
			gathered_data$Allele <- as.factor(coded)
		}else{
			coded <- as.character(gathered_data$Allele)
			indel_type <- substr(locus$var_allele, 1, 1)
			if (indel_type == "-") {
				coded <- substring(coded, 1, 1)
				coded[which( coded == "-" )] <- "1" ## Assumes that the first base captures the full del
				coded[which(coded != "1")] <- "0"	
				gathered_data$Allele <- as.factor(coded)
			
			}else if(indel_type == "+"){
				coded[which(coded == substring(locus$var_allele, 2))] <- "1" ## NB needs full nucleotide match to count
				coded[which(coded != "1")] <- "0"	
				gathered_data$Allele <- as.factor(coded)	
			
			}else{
				stop(cat("gatherFeaureSet(): Indel allele encoding not recognised for", paste(locus[1:4], collapse="|"), "\n"))
			}
			
		
		}
	
		## Remove MAPQ values that indicate mapping is unreliable ('A value 255 indicates that the mapping quality is not available')
		if( any( colnames(gathered_data) == 'Mapq' ) ){
			gathered_data$Mapq[gathered_data$Mapq > 250] <- NA
		}
		
		## Remove empty columns
		empty_col <- apply(X=gathered_data, MARGIN=2, FUN=function(v){return(all(is.na(v)))})
		if ( any(empty_col) ) {
			gathered_data <- gathered_data[, !empty_col]
		}
		
		## # Remove observations containing missing values
		gathered_data <- na.omit(gathered_data)
	
		## Gatekeeper 4: If at after filtering no variant alleles still exist...	
		if( !(any( gathered_data$Allele == 1 )) ){
			output_msg <- 'Skip:No variant alleles to test after filtering'
			return( list(output_df=tmp, output_msg=output_msg) )		
		}
	
		## Gatekeeper 5: If there are fewer than N variant alleles in the whole dataset...	
		if(  length(which(gathered_data$Allele == 1)) < min_allele_count ){
			output_msg <- paste("Skip:Less than", min_allele_count, "var alleles")
			return( list(output_df=tmp, output_msg=output_msg) )		
		}

		## Gatekeeper 6: If there are fewer than N reference alleles in the whole dataset...	
		if(  length(which(gathered_data$Allele == 0)) < min_allele_count ){
			output_msg <- paste("Skip:Less than", min_allele_count, "ref alleles")
			return( list(output_df=tmp, output_msg=output_msg) )		
		}		
		
		## Ensure data.frame encoding is correct
		noms <- colnames(gathered_data)
		if ( any(noms %in% "SampleIdx") ) {
			gathered_data$SampleIdx <- as.factor(gathered_data$SampleIdx)
		}
		if ( any(noms %in% "CyclePosition") ) {
			gathered_data$CyclePosition <- as.integer(gathered_data$CyclePosition)
		}
		if ( any(noms %in% "Mapq") ) {
			gathered_data$Mapq <- as.integer(gathered_data$Mapq)
		}		
		if ( any(noms %in% "Strand") ) {
			gathered_data$Strand <- as.logical(gathered_data$Strand)
		}			
		if ( any(noms %in% "QDiff") ) {
			gathered_data$QDiff <- as.numeric(gathered_data$QDiff)
		}			
		if ( any(noms %in% "Qual") ) {
			gathered_data$Qual <- as.integer(gathered_data$Qual)	
		}			
		if ( any(noms %in% "NumberOfCigarOps") ) {
			gathered_data$NumberOfCigarOps <- as.integer(gathered_data$NumberOfCigarOps)
		}			
		
		## TODO: Test for colinearity (VIF in package(car) is one option) - still problems.
		
		## TODO: Test for outliers (use outlierTest in package(car)); http://www.statmethods.net/stats/rdiagnostics.html
		
		## Add in sample categories (Tumour/Normal) where applicable
		## NO: Results in singularities again - have to test independently
		# if( !(all( is.na(bamSet$isTumour) ) ) ) {
			# tmp <- cbind(tmp, isTumour=as.factor(bamSet$isTumour[tmp$SampleIdx]) ## Match the BAM to its Tumour/Germline annotation
		# }

		return(list(output_df=gathered_data, output_msg='pass'))
	}
}

assessor <- function(z="SampleIdx") {
	## Assign model variables to one of four sets [A-D].
	## vector of strings -> char
	# Args:
	#   z: A comma-delimmited string of model variables.
	# Returns:
	# 		A character [A-D] 
	
	if( all(z %in% "SampleIdx") & length(z) == 1) {
		return("A")
	}else if(length(z) > 1 && z[1] == "SampleIdx"){
		return("B")
	}else if(length(z) > 1 && z[1] != "SampleIdx" && any(z %in% "SampleIdx") ){
		return("C")		
	}else if ( all(!(z %in% "SampleIdx")) ) {
		return("D")	
	}
}

logit2anova <- function(model=formula("Allele~."), x) {
	## Wrapper to perform logistic regression, followed by ANOVA-like test.
	# Args:
	#   model: An object of class formula. Input data (x) is inherited from current frame.
	# Returns:
	# 		A list, with fitted parameters (beta), and call summaries

	new_call <- NA
	new_call_summary <- NA
	glm_out <- glm(model, data=x, family="binomial") # Fit logistic regression model
	walds <- summary(glm_out)$coefficients # coefficients, their standard errors, the z-statistic (sometimes called a Wald z-statistic), and the associated p-values. The logistic regression coefficients give the change in the log odds of the outcome for a one unit increase in the predictor variable.
	output <- anova(glm_out, test="Chisq") # ANOVA-like test.

	output_pval_vector <- paste(rownames(output), format(output[, 'Pr(>Chi)'], digits=3), sep=":") ## Verbose output of actual p-values in a string
	output_pval_vector <- gsub(" ", "", paste(output_pval_vector[2:length(output_pval_vector)], collapse="|") )
	
	## Now test with ANOVA-like test.	
	THRESH <- 0.01	
	test_b <- which(output[ ,"Pr(>Chi)"] <= ( THRESH/ (nrow(output)-1) ) )
	if (length(test_b) > 0) {
		# If significant terms were found in this test, report these.
		tmp <- output[test_b, , drop=FALSE] ## Subset to sig. terms
		tmp <- tmp[order(tmp[, "Pr(>Chi)"]), ] # Re-order by beta (not lowest p-value, as unstable)
		new_call_summary <- paste(rownames(tmp), collapse=",", sep=",")
		new_call <- assessor(z=rownames(tmp)) # Custom function to assign model vars to a set 
	
	}else{
		new_call <- 'E'
		new_call_summary <- "No significant predictors"
	
	}
	return(list(Call=new_call, CallSummary=new_call_summary, p.vals=output_pval_vector))
}

stepLogisticRegression <- function(x=NA, features=c("SampleIdx", "CyclePosition", "Qual", "Mapq", "Strand") ) {
	## Wrapper to perform stepwise logistic regression.
	## (data.frame, vector of strings) -> vector of strings
	# Args:
	#   x: Output from gatherFeaureSet()
	#   features: A vector of features to be included in the models (must be in colnames(x))
	# Returns:
	# 		A named vector of strings containing the call (A-D), and details of the model fit
	
	
	## EXCEPTION [0]: No variants detected as expected, abort here and report.	
	if(length(unique(x$Allele)) == 1){		
		tmp <- paste(paste(features, 'NA', sep=":"), collapse="|")
		output_i <- c(Call=NA, CallSummary="No variants detected", p_vals=tmp)
		return(output_i)		
	}
	## TODO: Issue a warning if the variant allele is found in all samples
	# if( !( any(table(x[, 1:2]) == 0)) ){warning()}
	
	## Subset x to requested features
	features <- features[which(features %in% colnames(x))]
	x <- x[, c("Allele", features)]
	formula_to_test <- formula( paste ("Allele~", paste(features, collapse="+") ) )
	
	## Cross-tabulation of categorical variables
	## "Empty cells or small cells: You should check for empty or small cells by doing a crosstab between categorical predictors and the outcome variable. If a cell has very few cases (a small cell), the model may become unstable or it might not run at all." http://www.ats.ucla.edu/stat/r/dae/logit.htm
	# test_a <- any( xtabs(~Allele + SampleIdx + Strand, data = x) == 0 )
	
	## EXCEPTION [1]: If the variant allele is confined to one strand, abort here and report.
	## Perfect classfiers 'break' logistic regression (same goes for sampleIdx, see 'exception 2')
	## http://www.ats.ucla.edu/stat/mult_pkg/faq/general/complete_separation_logit_models.htm
	## http://stats.stackexchange.com/questions/11109/how-to-deal-with-perfect-separation-in-logistic-regression
	## http://cemsiis.meduniwien.ac.at/fileadmin/msi_akim/CeMSIIS/KB/volltexte/Heinze_Schemper_2002_Statistics_in_Medicine.pdf
	if( any("Strand" == features) ) {
		strand_test <- xtabs(~Allele + Strand, data = x)
		if ( any(strand_test["1", ] == 0) || ncol(strand_test) == 1){
			tmp <- paste(paste(features, 'NA', sep=":"), collapse="|")
			output_i <- c(Call="D", CallSummary="Variant restricted to one strand", p_vals=tmp)
			return(output_i)
		}
	}

	
	# Simple Model without interaction terms, or second order terms (for speed at cost of bias)
	options( warn = 2 ) ## Change warnings to errors, to catch errors due to perfect seperation and non-convergence.
	bicglm1 <- try(bic.glm(Allele ~., data=x, glm.family="binomial", factor.type = TRUE, strict=FALSE), silent = TRUE)	## EV in output is the beta-value for each feature.	
	options( warn = 0 )
	has_bicglm_failed <- class(bicglm1) == "try-error"

	## EXCEPTION [2]: If bicglm has failed (e.g. Dealing with the case when sample ID perfectly correlates with variant (may happen often with clean data, when there are only two samples))
	## NB Could test seperation formally with separation.detection() in library(brglm)
	if ( has_bicglm_failed ){
		## Firth logistic regression, no stepwise feature selection (forward()/backward() too buggy here)
		## Instead, could try basic nested models instead + min(AIC)...	
		# nested_formula_list <- list() ## This is NOT the same as forward/backward selection
		# for(f in 1:length(features)) {
			# nested_formula_list[[f]] <- formula( paste ("Allele~", paste(features[1:f], collapse="+") ) )
		# }
		# fit_logistf_nested = lapply(X=nested_formula_list, FUN=brglm,  family = binomial, data = x, method = "brglm.fit")	
		fit_logistf <- brglm(formula_to_test, family = binomial, data = x, method = "brglm.fit")
		fit_logistf_sum <- summary(fit_logistf)
		walds <- fit_logistf_sum$coefficients
		
		output_pval_vector <- paste(rownames(walds), format(walds[, 'Pr(>|z|)'], digits=3), sep=":") ## Verbose output of actual p-values in a string
		output_pval_vector <- paste(output_pval_vector[2:length(output_pval_vector)], collapse="|")	
		
		new_p <- walds[, "Pr(>|z|)"]
		names(new_p) <- gsub('StrandTRUE', 'Strand', names(new_p))
		new_p2 <- c(
					SampleIdx=min( new_p[ grep('SampleIdx', names(new_p)) ] ), 
					new_p[features] 
					)
		new_p2 <- na.omit(new_p2)
		new_p2[ which(new_p2 > 0.05) ] <- NA ## only return significant p-values
		new_p2 <- sort(new_p2)
		
		if( length(new_p2) >= 1 ) {
			is_A <- names(new_p2)[1] == 'SampleIdx' & length(new_p2) == 1
			is_B <- names(new_p2)[1] == 'SampleIdx' & length(new_p2) > 1
			is_C <- names(new_p2)[1] != 'SampleIdx' & any(names(new_p2) == 'SampleIdx') 
			is_D <- all(names(new_p2) != 'SampleIdx') 
		}else{
			is_A <- is_B <- is_C <- is_D <- FALSE
		}
		
		annot_tmp <- paste(names(na.omit(new_p2)), collapse=",") ## For CallSummary
		
		if( is_A ){
			## SampleIdx is only predictor
			output_i <- c(Call="A", CallSummary=annot_tmp, p_vals=output_pval_vector)
			return(output_i)
		
		}else if( is_B ){
			## SampleIdx is best (but not only) predictor 
			output_i <- c(Call="B", CallSummary=annot_tmp, p_vals=output_pval_vector)
			return(output_i)	
			
		}else if( is_C ){
			## SampleIdx is not best predictor 
			output_i <- c(Call="C", CallSummary=annot_tmp, p_vals=output_pval_vector)
			return(output_i)	
			
		}else if( is_D ){
			## SampleIdx is not a predictor 	
			output_i <- c(Call="D", CallSummary=annot_tmp, p_vals=output_pval_vector)
			return(output_i)			
		
		}else if( all(!(c(is_A, is_B, is_C, is_D))) ) {
			## No predictors
			output_i <- c(Call="E", CallSummary='No significant predictors', p_vals=output_pval_vector)
			return(output_i)			
		
		}
	}		
	
	## If bic.glm worked without errors...
	if ( !(has_bicglm_failed) ) {
		## Check a feature exists with posterior prob. == 100%
		bic_tab <- bicglm1$probne0	 # the posterior probability that each variable is non-zero
		names(bic_tab) <- bicglm1$namesx
		names(bic_tab) <- gsub("TRUE", "", names(bic_tab) )
		bic_pass <- which(bic_tab == 100)
				
		if (length(bic_pass) == 0 & bicglm1$label[1] == "NULL") {
			# No feature could be selected, and none of your predictors appear to relate to 'Allele', i.e. nothing was tested.
			l2a_out <- logit2anova(model=Allele ~., x) # logistic regression/ANOVA combo
			
		}else if (length(bic_pass) == 0  & bicglm1$label[1] != "NULL") {
			# no strong feature could be selected, but a model was suggested...
			tmp <- unlist(strsplit(bicglm1$label[1], ","))	
			tmp <- gsub("TRUE", "", tmp)	
			new_form <- formula( paste ("Allele~", paste(tmp, collapse="+") ) )
			l2a_out <- logit2anova(model=new_form, x) # logistic regression/ANOVA combo			
	
		}else{
			# If an optimal feature[s] was chosen, assess (same as above)
			tmp <- unlist(strsplit(bicglm1$label[1], ","))	
			tmp <- gsub("TRUE", "", tmp)	
			new_form <- formula( paste ("Allele~", paste(tmp, collapse="+") ) )
			l2a_out <- logit2anova(model=new_form, x) # logistic regression/ANOVA combo
			
		}
		### Update output data.frame
		output_i <- c(Call=l2a_out$Call, CallSummary=l2a_out$CallSummary, p_vals=l2a_out$p.vals)
	}

	## Model QC ##
	## [QC] Cross-validation, for model accuracy (a problem in sparse datasets (lowFreqVar))
	if (FALSE) {
		# diagnostic plots
		plot(residuals(glm_out,type="pearson") ~ x$Allele, main="pearson residual vs Allele plot") 
		# "Learning curves" to identify high variance/bias	
	}
	return(output_i)
} 

runVarSLR <- function(x=1, write2disk=FALSE) {
	## Wrapper function
	## (integer, data.frame) -> list
	# Args:
	#   x: Index into (row) loci object 
	# Returns:
	#  list with output results for the specified locus

	## Gather read-derived data
	training_data <- gatherFeaureSet( locus_idx_1=x, d=10000, min_allele_count=4 ) 

	gatherFeaureSet_call <- training_data$output_msg
	training_data <- training_data$output_df
	
	if(write2disk){
		write.table(training_data, file=paste0("/farm/home/salm01/public_html/LIVE/MSalm/raw/", x, ".varslr.in"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
	}
	
	## prepare output object
	output_list <- list()
	for(cc in 1:4){
		output_list[colnames(VARIANTS)[cc]] <- VARIANTS[x, cc]
	}
	output_list['Call'] <- NA
	output_list['CallSummary'] <- NA
	output_list['p_vals'] <- NA
	
	# if ( gatherFeaureSet_call == 'pass' ) {
		# ## Allele data was retrieved.
		
		# if ( VARIANTS[x, 'is_SNV'] ) {
			# feat_to_test=c("SampleIdx", "CyclePosition", "Qual", "Mapq", "Strand", "MMQS") 			
		# }else{
			# ## Note exclusion of mapq (confounder for larger indels; would have to re-calculate taking indel into account).
			# feat_to_test=c("SampleIdx", "CyclePosition", "Qual", "Strand", "QDiff", "MMQS")
		# }

		# outcome <- stepLogisticRegression( x=training_data , features=feat_to_test )
		
		# ## Cross-validation - for design stage only
		# if ( FALSE ) {
			# warning("Cross-validation chunk is active!!!\n")
			# require(boot)
			# ## 10-fold Cross-Validation of the logistic regression model selected
			# ## NB This will obviously not hold if Firth logistic regression was used.
			# model_before <- as.formula( paste("Allele ~ ", paste(feat_to_test, collapse="+") ) )
			# model_after <- as.formula( paste("Allele ~ ", paste(unlist(strsplit(outcome["CallSummary"], ",")), collapse="+") ) )
			# glm_fit_before <- glm(model_before, data=training_data, family="binomial")
			# glm_fit_after <- glm(model_after, data=training_data, family="binomial")
			# cv_error_before <- cv.glm(training_data, glm_fit_before, K=10)$delta[1] ## 10-fold CV
			# cv_error_after <- cv.glm(training_data, glm_fit_after, K=10)$delta[1] ## 10-fold CV
			# ## NB Expect cv_error_before ~= cv_error_after, as this suggests BIC has done a good job in feature selection, retaining the predictive power of the full model.
			# output_list$Call <- cv_error_before
			# output_list$CallSummary <- cv_error_after
			# output_list$p_vals <-  outcome['CallSummary']
			
			
		# }else{
			# ## Add results to output
			# for(cc in 1:length(outcome)) {
				# col_idx <- names(outcome)[cc]
				# if ( !(is.na( col_idx )) ) {
					# output_list[col_idx] <- outcome[cc]
				# }
			
			# }
		# }
		# return(output_list)
		
	# }else{
		# ## No allele data was retrieved.
		# output_list['CallSummary'] <- gatherFeaureSet_call
		# return(output_list)
		
	# }
}

retrieveTable <- function(x) {
	## Gather output data from parallel run and return a user-freindly version
	# Args:
	#   x: Ouput of mclapply()
	# Returns:
	#   a matrix

	output <- vector() # an empty matrix
	for (i in 1:length(x)) {
		if (class(x[[i]]) != "try-error") {
			tmp <- x[[i]][c('chromosome', 'position', 'ref_allele', 'var_allele', 'Call', 'CallSummary', 'p_vals')]
			output <- rbind(output, unlist(tmp))
		}
	}
	return(output)
}


getErrors <- function(x) {
	## Gather output data from parallel run and return all submissions with errors
	# Args:
	#   x: Ouput of mclapply()
	# Returns:
	#   a matrix

	output <- vector() # an empty matrix
	for (i in 1:length(x)) {
		if (class(x[[i]]) == "try-error") {
			tmp <- c(i,
					as.character(attr(x[[i]], "condition")$message),
					as.character(attr(x[[i]], "condition")$call)
					)
			output <- rbind(output, tmp)
		}
	}
	return(output)

}

#########################
### Retired functions ###
#########################
if(FALSE) {

	makeBamSet <- function(germline_bam, tumour_bams) {
		## Convenience wrapper to make BAM list object
		# (string, vector) -> list
		# Args:
		#	germline_bam: Full path to germline BAM [string]
		#   tumour_bams: Full path to tumour BAM file(s) [vector of strings]
		#	open_connection: Open connections to BAMs
		# Returns:
		#   A list
		
		if(length(germline_bam) != 1){
			stop('makeBamSet(): Please specify a single reference BAM file\n')
		}
		if(length(germline_bam) < 1) {
			stop('makeBamSet(): Please specify at least one tumour BAM file\n')
		}
		
		if(any(ls(envir=.GlobalEnv, pattern="BAMSET") %in% "BAMSET")){
			warning("An object named BAMSET already existed in GlobalEnv and has been overwritten!\n")
		}
		
		bams <- c(germline_bam, tumour_bams) ## This ensures that 'germline' is the baseline sample for logistic regression later (ie always sample == 1)
		## Test if BAM/BAI exist
		missing_bam <- !(file.exists(bams))
		if( any( missing_bam ) ){
			cat("Missing file:", bams[missing_bam], "\n")
			stop("makeBamSet(): BAM not found\n")
		}
		
		bais <- paste(bams, ".bai", sep="")
		missing_bai <- !(file.exists( bais ))
		if( any( missing_bai ) ){
			cat("Missing file:", bais[missing_bai], "\n")
			stop("makeBamSet(): BAI not found\n")
		}
		
		## Return list (and assign representation to global env as well)
		BAMSET <- list( bams=bams,
					 bai=bais)
					
		assign("BAMSET", BAMSET, envir=globalenv())
		
		return( BAMSET )
	}

	base2ascii <- function(w=40) { 
		## Convert  phred-scaled base error probability ASCII-encoding (+33).
		## int -> char
		## See, http://datadebrief.blogspot.co.uk/2011/03/ascii-code-table-in-r.html
		# Args:
		#   w: Phred-scaled base error probability 
		# Returns:
		#     A single ASCII character.
		if (length(w) == 1 && is.na(w)) {
			return(w)
		}else if(length(w) != 1) {
			stop("Only one element at a time please!")
		}else{
			return( rawToChar(as.raw( w + 33 )) )
		}
	}

	expandAlign <- function(r=dfr[1, ]) {
		## Convert an alignment object into an alignment matrix.
		## matrix.row -> matrix
		# Args:
		#   r: A single alignment object (see getSNVData2)
		# Returns:
		#   An Mx4 matrix, where each row is a position in the read. Columns refer to BASE, QUAL, Index in Read & Position in Genome, respectively.
		
		## N.B. Strand info (encoded in FLAG) does not alter alignment info as it is encoded relative to the reference (i.e. interpret CIGAR from left to right, and read position as incrementing, irrespective of strand info). 
		
		
		Ci=getCigar(as.character(r$cigar))
		# Expand out CIGAR info
		ci <- rep(as.character(Ci[,"Type"]), Ci[, "Length"]) 
		# Create starting matrix of BASE & QUAL
		si_qi <- cbind( B=unlist(strsplit(as.character(r$seq), "")), Q=unlist(strsplit(as.character(r$qual), "")) )
		si_qi <- data.frame( si_qi, idx_q=c(1:nrow(si_qi)), idx_s=c(as.numeric(r$position): (as.numeric(r$position) + nrow(si_qi) - 1) ) ) 
		
		# Perform CIGAR Ops. SLOW: could be more imaginatively/concisely expressed, but at least it's clear.
		k <- 1
		for (i in ci) {
			if(i == "I" | i == "S"){
				# The subject index is affected by I or S (info downstream shifted along)
				targ <- (k):nrow(si_qi) 
				si_qi$idx_s[targ] <- si_qi$idx_s[targ] - 1
				si_qi$idx_s[k] <- NA
				k <- k + 1
				
			}else if(i == "D"){
				# The query index, B and Q is affected by D, so add in a blank row
				p1 <- 1:(k-1)
				p2 <- k:nrow(si_qi) 
				si_qi <- rbind(si_qi[p1, ],
								rep(NA, 4),
								si_qi[p2, ])
				# Update idx_s
				si_qi$idx_s[k:nrow(si_qi)] <- seq(si_qi$idx_s[k+1], si_qi$idx_s[k+1]+(nrow(si_qi) - k), 1)
				k <- k + 1

			}else if(i == "H"){
				stop("Hard-clipping not expected!\n")
				k <- k + 1				
			}else if(i == "P"){
				stop("Padding not expected!\n")
				k <- k + 1				
			}else if(i == "N"){
				stop("N-padding not expected!\n")
				k <- k + 1				
			}else{
				k <- k + 1
				next()
			}
		}	
		return(si_qi)
	}

	posIdx <- function(r=dfr[1, ], p=pos) {
		## Retrieve base specific features at SNVs
		## matrix.row, int -> list
		# Args:
		#   r: An alignment object (see getSNVData2)
		#	p: A genomic position [int]. N.B. no need to specify chromosome as it is implicit in "r"
		# Returns:
		#   A list, B=BASE, Q=QUAL and I=INDEX values for the specified position
		
		if( as.character(r$nCigar) == 1 )  {
			## No CIGAR string operation, return the base/qual/cycle at position 
			idx <- (p - r$position) + 1 # "Fence-post problem"
			if ( idx < 1 | idx > nchar(as.character(r$seq)) ) {
				#warning("Genome position not covered by read\n")
				return(list(B=NA, Q=NA, I=NA))
			}else{
				B <- unlist(strsplit(as.character(r$seq), ""))[idx]  ## Base
				Q <- unlist(strsplit(as.character(r$qual), ""))[idx] ## Qual score
							
				if (r$revstrand) {
					idx <- nchar(as.character(r$seq)) - idx + 1 # If on reverse strand, flip index...
				}
				return(list(B=B, Q=ascii2base(Q), I=idx))
			}
		}else{
		## CIGAR string operation...SLOW!
		## Could also use sequenceLayer() in Rsamtools.
			ear <- expandAlign(r)
			hit <- which(ear$idx_s == p)
			if (length(hit) > 0) {
				B <- as.character(ear$B[hit])
				Q <- as.character(ear$Q[hit])
				idx <- ear$idx_q[hit]
				# dist_from_indel <- min(abs(which(is.na(ear$idx_s) | is.na(ear$idx_q)) - hit)) ## Closest indel
				if (r$revstrand) {
					idx <- nchar(as.character(r$seq)) - idx + 1 # If on reverse strand, flip index...
				}			
			}else{
				#warning("Genome position not covered by read\n")
				return(list(B=NA, Q=NA, I=NA))
			}

		}
		return(list(B=B, Q=ascii2base(Q), I=idx))
	}



}
