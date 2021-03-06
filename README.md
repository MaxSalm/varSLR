varSLR
======

Accurate mutation calling in clinical tumour samples remains a formidable challenge, confounded by sample complexity, experimental artefacts and algorithmic constraints. In particular, high sequencing error rates (~0.1-1x10-2 per base) entail costly manual review of putative mutations followed by orthogonal validation. Efficient filtering is currently required, given that most mutation callers identify many thousands (in exome sequencing), if not millions (in whole genome sequencing) of candidate mutations per experiment.

To aid in this process, we developed the open-source VarSLR R package to identify somatic nucleotide and insertion-deletion mutations that are likely to be sequencing artefacts. The algorithm incorporates putative confounders of call accuracy (Genome Res. Mar. 2012;22(3):568-576) into stepwise logistic regression models and subsequently classifies variants within a simple, 4-tiered quality schema.

VarSLR is highly scalable and designed to be run in an 'embarrassingly parallel' fashion, thus benefiting from high-performance computing facilities. Moreover, VarSLR performed with high precision when tested with synthetic and experimental data, and has been successfully applied to numerous projects (e.g. Science, Oct. 2014; 346(6206): 251-6)

See <https://www.crick.ac.uk/research/science-technology-platforms/bioinformatics-and-biostatistics/> for more details.

To run the *varSLR* pipeline, execute the following:

    ## Load libraries
    library(parallel)
    source("functions_v3.r")

    ## Load SNV/Indel input file
    # The expected column input order is chr, position, reference allele, alternate allele.
    snv2 <- as.matrix(read.delim("chr5.posMerged.txt", header=FALSE))
    variants <- variantConstructor(x=snv2)
    variants <- filterBlacklist(x=variants)

    ## Define BAM files
    bamSet <- list( bams=paste(c(1:9), ".cleaned.bam", sep=""),
                          bai=paste(c(1:9), ".cleaned.bam.bai", sep=""))
            
    ## Run varSLR
    tst2 <- mclapply(X=test_v, FUN=runVarSLR, mc.cores=CORES) 
    tst3 <- retrieveTable(tst2)

Features
--------

-   Accurate
-   Reproducible
-   Relatively fast (\_8s/variant)
-   Scalable
-   SNVs and Indels
-   Independent of mutation caller

Installation
------------

Contribute
----------

-   Issue Tracker: <https://github.com/MaxSalm/varSLR/issues>
-   Source Code: <https://github.com/MaxSalm/varSLR>

Support
-------

If you are having issues, please let us know via the Issue Tracker and we will endeavour to resolve them as soon as possible.

Known Bugs
----------

License
-------

The project is licensed under the GPL license.

Citation
--------
