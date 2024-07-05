# scRNA sequencing analysis

scRNA-seq analysis of T cells from PBMCs stimulated with different antigens/activators.

Memory CD4 T cells were isolated and sequenced using the platform Seq-Well S3 (doi: 10.1007/978-1-0716-2756-3_3.) and the High-Output Kit  for Illumina NextSeq 550.
The demultiplexing and convertion of BCL files to FASTQ were done with bcl2fastq. Then, the QC was done with FastQC and MultiQC. 
The identification of barcodes and UMI was done with umi_tools. The reads were map with STAR and the assigning of reads to genes with featureCounts.

The QC of the cells was done in R using Seurat and auxiliary libraries (file: scRNAseq_seq_well_analysis_QC_cleaning.r) and the cell population analysis with Sparse PCA, Module score and pathways analysis can be found in scRNAseq_seq_well_analysis.r. 
