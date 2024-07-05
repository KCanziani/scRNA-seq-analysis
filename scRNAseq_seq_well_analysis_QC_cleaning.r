library(Seurat)
library(dplyr)
library(biomaRt)
library(tibble)
library(ggplot2)
library(biomaRt)
library(stringr)
library(gridExtra)
library(RColorBrewer)

# Load the functions from the external file
source("R/auxiliary_functions_seq_well.r")

# Read count matrix
# Set the working directory to the location of your data files
setwd("/User/alignment/bam/sorted")
set.seed(533)

# List of count matrix filenames
data <- c("Medium_S1_R2_001_Aligned_sorted.bam_counts.tsv.gz", 
        "CD3_CD28_S2_R2_001_Aligned_sorted.bam_counts.tsv.gz",
        "Allergen_S3_R2_001_Aligned_sorted.bam_counts.tsv.gz",
        "BP_S4_R2_001_Aligned_sorted.bam_counts.tsv.gz")

# Create Seurat objects using the provided data files
seurat_objects <- create_seurat_obj(data)

# Merge the seurat objects
merged_seurat_object <- Merge_Seurat_List(seurat_objects, add.cell.ids = c("Medium", "CD3_CD28", "Allergen", "BP"), project = "PBMC_stimulation" )

# Create features to evaluate quality
merged_seurat_object$percent.mt <- PercentageFeatureSet(merged_seurat_object, pattern = "^MT-")

merged_seurat_object$percent.rb <- PercentageFeatureSet(merged_seurat_object, pattern = "^RP[SL]")

# Percentage hemoglobin genes - includes all genes starting with HB except HBP and Platelet.
merged_seurat_object <- PercentageFeatureSet(merged_seurat_object, "^HB[^(P)]", col.name = "percent_hb")

merged_seurat_object <- PercentageFeatureSet(merged_seurat_object, "PECAM1|PF4", col.name = "percent_plat")

#Rename features
merged_seurat_object$nUMI <- merged_seurat_object$nCount_RNA    

merged_seurat_object$nCount_RNA <- NULL

merged_seurat_object$nGene <- merged_seurat_object$nFeature_RNA  

merged_seurat_object$nFeature_RNA <- NULL

merged_seurat_object$log10GenesPerUMI <- log10(merged_seurat_object$nGene) / log10(merged_seurat_object$nUMI)


###Exploratory analysis to check quality
# Volcano Plot of nGene, nUMI, percent.mt, and percent.rb
VlnPlot(merged_seurat_object, features = c("nGene","nUMI","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

# Feature Scatter Plot: nUMI vs percent.mt
FeatureScatter(merged_seurat_object, feature1 = "nUMI", feature2 = "percent.mt")

# Feature Scatter Plot: nUMI vs nGene
FeatureScatter(merged_seurat_object, feature1 = "nUMI", feature2 = "nGene")

# Feature Scatter Plot: nUMI vs percent.rb
FeatureScatter(merged_seurat_object, feature1 = "nUMI", feature2 = "percent.rb")

# Feature Scatter Plot: percent.rb vs percent.mt
FeatureScatter(merged_seurat_object, feature1 = "percent.rb", feature2 = "percent.mt")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
merged_seurat_object@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.mt)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
merged_seurat_object@meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=percent.rb)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 300) +
  geom_hline(yintercept = 250) +
  facet_wrap(~orig.ident)

# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
merged_seurat_object@meta.data %>%
  ggplot(aes(x=log10GenesPerUMI , color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) +
  xlim(0.75, 1)


# Visualize the relationship between nUMI and cell type
merged_seurat_object@meta.data %>%
  ggplot(aes(x=nUMI , color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 500) +
  scale_x_log10()


# Visualize the relationship between nGene and cell type
merged_seurat_object@meta.data %>%
  ggplot(aes(x=nGene , color = orig.ident, fill=orig.ident)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 250) +
  scale_x_log10()


# Filter bad quality data
merged_seurat_object_clean <- subset(merged_seurat_object, subset= (nUMI>= 500) & 
                                  (nGene >= 250) & 
                                  (log10GenesPerUMI > 0.80) & 
                                  (percent.mt < 20))

## Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = merged_seurat_object_clean, slot = "counts")
# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0
# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10
# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]
# Reassign to filtered Seurat object
filtered_sample <- CreateSeuratObject(filtered_counts, meta.data = merged_seurat_object_clean@meta.data)

saveRDS(filtered_sample, file = "merged_seurat_object.rds")
#filtered_sample <- readRDS("merged_seurat_object.rds")

# Compute the relative expression of each gene per cell Use sparse matrix
# operations, if your dataset is large, doing matrix devisions the regular way
# will take a very long time.
par(mar = c(4, 8, 2, 1))
C <- filtered_sample@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]

boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell",
    col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)


## Visualize the number of cell counts per sample after filtering
filtered_sample@meta.data %>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  scale_fill_discrete(limits = c("Medium", "Allergen", "CD3_CD28", "BP")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Cell number")


## Plot filtered QC
feats <- c("nGene","nUMI","percent_mito","percent_ribo", "percent_hb", "percent_plat")
VlnPlot(filtered_sample, features = feats, ncol = 3, pt.size = 0.1) + 
  theme(plot.title = element_text(size=10)) +
  NoLegend()


## Filter overexpressed genes
dim(filtered_sample)

# Filter MALAT1
filtered_sample <- filtered_sample[!grepl("MALAT1", rownames(filtered_sample)), ]

# Filter Mitocondrial (optional if you do not care about them)
#filtered_sample <- filtered_sample[!grepl("^MT-", rownames(filtered_sample)), ]

# Filter Ribossomal gene (optional if that is a problem on your data) filtered_sample
filtered_sample <- filtered_sample[ ! grepl('^RP[SL]', rownames(filtered_sample)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
#filtered_sample <- filtered_sample[!grepl("^HB[^(P)]", rownames(filtered_sample)), ]


## Cell cycle score
# Before running CellCycleScoring the data need to be normalized and logtransformed.
data.filt <- NormalizeData(filtered_sample)

data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
    s.features = cc.genes$s.genes)

#Check if there are difference in the Scores across samples
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "orig.ident",
    ncol = 4, pt.size = 0.1)


## Predict doublets in each experiment, that is why we split the data into the 4 conditions

split_seurat <- SplitObject(filtered_sample, split.by = "orig.ident")
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2") %>%
            RunPCA(npcs = 30)
})


# Initialize lists to store results
sweep.res <- list()
sweep.stats <- list()
bcmvn <- list()
barplot <- list()
pK <- list()

# Loop through split_seurat
for (i in seq_along(split_seurat)) {
    # Run parameter optimization
    sweep.res[[i]] <- paramSweep_v3(split_seurat[[i]], PCs = 1:10, sct = TRUE)
  
    # Summarize the parameter sweep results
    sweep.stats[[i]] <- summarizeSweep(sweep.res[[i]], GT = FALSE)
  
    # Find the best pK value
    bcmvn[[i]] <- find.pK(sweep.stats[[i]])
  
    pK[[i]] <- bcmvn[[i]] %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
                filter(BCmetric == max(BCmetric)) %>%
                pull(pK) 
    pK[[i]] <- as.numeric(as.character(pK[[i]][[1]]))
}


nExp <- list()
doublets <- list()

# define the expected number of doublet cellscells.
for (i in seq_along(split_seurat)){
    nExp[[i]] <- round(ncol(split_seurat[[i]]) * 0.0237)  # expect 2.37% doublets in seq-well
    doublets[[i]] <- doubletFinder_v3(split_seurat[[i]], pN = 0.25, pK = pK[[i]], 
                                        nExp = nExp, PCs = 1:10, sct = TRUE)
  
}

# Check and save results
for (i in seq_along(split_seurat)){
    plot <- ggplot(bcmvn[[i]], aes(pK, BCmetric, group = 1)) +
                geom_point() +
                geom_line()
    ggsave(paste0(i, "_doublet_finder.png"), plot)
}
#Conclusion: No doublets were found



#Data integration
#Identify integration features
features <- SelectIntegrationFeatures(object.list = split_seurat,
                                  dims = 1:30, nfeatures = 2000) 
#Import immune related genes from ImmPort
immune_genes <- read.csv("ImmPort_immune_genes.csv") 

#Define genes to integrate (include gene of interest)
genes_to_integrate <- unique(append(immune_genes$Symbol, 
                                    c(features, 
                                    transcription_factors, #List of extra genes not included in ImmPort, from source
                                    extra_genes))) #List of extra genes not included in ImmPort, from source

integrated_seurat <- PrepSCTIntegration(object.list = split_seurat, anchor.features = genes_to_integrate)

immune.anchors <- FindIntegrationAnchors(object.list = integrated_seurat , normalization.method = "SCT",
                  anchor.features = features)

integrated_seurat <- IntegrateData(anchorset= immune.anchors, normalization.method = "SCT")

# Run the standard workflow for visualization and clustering
integrated_seurat <- RunPCA(integrated_seurat)
ElbowPlot(integrated_seurat) #Define the ndims

integrated_seurat<- RunUMAP(integrated_seurat, dims=1:10)

integrated_seurat<- FindNeighbors(integrated_seurat, reduction= "pca", dims= 1:10)

integrated_seurat<- FindClusters(integrated_seurat, resolution= c(0.1, 0.2, 0.3, 0.4, 0.5))

saveRDS(integrated_seurat, file = "integrated_seurat.rds")

## Plot the UMAP results and check correct clustering
dimplot1 <- DimPlot(integrated_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle("Clusters") + theme(plot.title = element_text(hjust = 0.5))
dimplot2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "orig.ident") +
  ggtitle("Samples") + theme(plot.title = element_text(hjust = 0.5))

# Arrange the DimPlot objects side by side
grid.arrange(dimplot1, dimplot2, ncol = 2)


DimPlot(integrated_seurat, reduction = "umap", split.by = "orig.ident") +
  ggtitle("Samples") + theme(plot.title = element_text(hjust = 0.5))


###Choosing resolution (number of clusters)
plot1 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.1") +
  ggtitle("RNA_snn_res.0.10") + theme(plot.title = element_text(hjust = 0.5))

plot2 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.2") +
  ggtitle("RNA_snn_res.0.20") + theme(plot.title = element_text(hjust = 0.5))

plot3 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.3") +
  ggtitle("RNA_snn_res.0.30") + theme(plot.title = element_text(hjust = 0.5))

plot4 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.4") +
  ggtitle("RNA_snn_res.0.40") + theme(plot.title = element_text(hjust = 0.5))

plot5 <- DimPlot(integrated_seurat, reduction = "umap", group.by = "integrated_snn_res.0.5") +
  ggtitle("RNA_snn_res.0.50") + theme(plot.title = element_text(hjust = 0.5))

grid.arrange(plot1, plot2, plot3, plot4, plot5, ncol = 2)