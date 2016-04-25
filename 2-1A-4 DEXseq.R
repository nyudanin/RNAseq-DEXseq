################################  PARATHYROIDSE ################################
biocLite(c("GenomicFeatures", "parathyroidSE", "GenomicAlignments"))
library("GenomicFeatures")
library("parathyroidSE")
library("GenomicAlignments")

## Get aligned BAM files
setwd("/Volumes/IBD/Yudanin/RNAseq/2-1A-4 RNAseq/2-1A-4 Aligned Files")
files <- list.files(pattern = "sorted.bam")
files <- files[-grep("bai", files)]

## Download mouse exon annotation database from UCSC
## Get list of annotation tables: supportedUCSCtables() 

mme <- makeTxDbFromUCSC(
  genome="mm10",
  tablename="ensGene",
  transcript_ids=NULL,
  circ_seqs=DEFAULT_CIRC_SEQS,
  url="http://genome.ucsc.edu/cgi-bin/",
  goldenPath_url="http://hgdownload.cse.ucsc.edu/goldenPath",
  taxonomyId=NA,
  miRBaseBuild=NA)

## Count reads by exon
exonicParts <- disjointExons(mme, aggregateGenes=FALSE)
bamlist <- BamFileList( files, index=character(), yieldSize=100000, obeyQname=TRUE )

exonsSE <- summarizeOverlaps( exonicParts, 
                              bamlist, mode="Union", 
                              singleEnd=FALSE,
                              ignore.strand=TRUE, 
                              inter.feature=FALSE, 
                              fragments=TRUE )

################################  DEXSEQ  ################################
biocLite("DEXSeq")
library("DEXSeq")

setwd("/Volumes/IBD/Yudanin/RNAseq/2-1A-4 RNAseq/2-1A-4 DEXseq")
## Assign appropriate sample names to the SE colData --------------------
rownames <- rownames(colData(exonsSE))
rownames <- gsub("2-1A-4 |_sorted.bam","",rownames)
rownames <- gsub(" ","_",rownames)
rownames(colData(exonsSE)) <- rownames

## Assign condition names to the SE colData --------------------
condition <- gsub("R._|D._","",rownames)
colData(exonsSE)$condition <- as.factor(condition)

## Create DEXSeqDataSet from SE --------------------
dxd <- DEXSeqDataSetFromSE(exonsSE, design = ~sample + exon + condition:exon)

## Normalization and dispersion estimates-----------------------------
BPPARAM <- MulticoreParam(workers=4)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, BPPARAM=BPPARAM)
plotDispEsts( dxd )

## Test for differential exon usage and calculate fold change-----------------------------
dxd <- testForDEU( dxd, BPPARAM = BPPARAM)
dxd <- estimateExonFoldChanges( dxd, fitExpToVar="condition", denominator = "MWAT", BPPARAM = BPPARAM)

dxr1 <- DEXSeqResults( dxd )

plotDEXSeq( dxr1, "ENSMUSG00000030669", displayTranscripts=TRUE, splicing=TRUE, names=TRUE, fitExpToVar="condition",
            expression=FALSE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )

DEXSeqHTML( dxr1, FDR=0.1, color=c("#E7A626","#4266F6","#269040","#C7302A"), BPPARAM = BPPARAM )

write.csv (dxr1[dxr1$groupID == "ENSMUSG00000030669",], file="Calca.csv")

?estimateExonFoldChanges
