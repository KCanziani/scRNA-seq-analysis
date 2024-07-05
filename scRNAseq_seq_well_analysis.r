library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(DESeq2)
library(biomaRt)
library(tibble)
library(scCustomize)
library(ggplot2)
library(biomaRt)
library(stringr)
library(gridExtra)
library(RColorBrewer)
library(PMA)
library(SCPA)
library(msigdbr)
library(purrr)
library(ReactomeGSA)
library(magrittr)

# Load the functions from the external file
source("R/auxiliary_functions_seq_well.r")

dir_plot <- "/User/Seq-well/Plot"

#Load data
integrated_seurat <- readRDS("integrated_seurat.rds")

immune_genes <- read.csv("ImmPort_immune_genes.csv")

genelist <- unique(c(immune_genes$Symbol, transcription_factors, extra_genes))

## Differential expression analysis between conditions
conditions <- c("Medium", "Allergen", "CD3_CD28", "BP")

# Initialize list to store differential gene expression results
list_DGE <- list()

# Perform differential expression analysis between conditions
for (condition_1 in conditions) {
    for (condition_2 in conditions) {
        if (condition_1 != condition_2) {
            # Generate the plot name with file path
            plot_name <- file.path(dir_plot, paste0(condition_1, "_vs_", condition_2, ".png"))
            
            # Perform differential gene expression analysis and save the results
            list_DGE[[paste0(condition_1, "_vs_", condition_2)]] <- differential_gene_expression(
                seurat_object = integrated_seurat, 
                condition_1 = condition_1, 
                condition_2 = condition_2, 
                variable_to_group_by = "orig.ident", 
                genes_to_plot = genelist,
                name_and_path_plot = plot_name
            )
        }
    }
}
        


## Sparse PCA to find Module scores and identify different types of cells
# Identify the most variable genes (default = 2000)
variable_genes <- VariableFeatures(integrated_seurat)
immune_genes <- read.csv("ImmPort_immune_genes.csv")

genelist <- unique(c(immune_genes$Symbol, variable_genes, transcription_factors, extra_genes))

tGenes <- t(as.matrix(GetAssayData(integrated_seurat)[which(rownames(integrated_seurat) %in% genelist),]))

# Remove the lowest-variance genes (1st percentile and below)
# This is to remove immune genes that are not sufficiently detected in our data
variances <- apply(tGenes, 2, var)
tGenes1 <- tGenes[,(variances > quantile(variances, 0.01))]

# Arguments:
# x: A matrix with n rows (samples), p columns (genes). Should be scaled.
# sumabsv: Penalty to achieve sparseness (smaller sumabsv = sparser). Must be between 1 and sqrt(p).
# K: Desired number of components. Must be between 1 and p.
# orth: If TRUE, makes the components pseudo-orthogonal.
# niter: Number of iterations to run per component. Default is 20.
out.orth <- SPC(x = scale(tGenes1), sumabsv=3, K=50, orth=TRUE, niter=20, vpos = TRUE)

out.orth$genenames <- colnames(tGenes1)
# Invert any majority-negative components to be majority-positive
signs <- apply(out.orth$v, 2, function(x) sum(x<0)/sum(x!=0))
out.orth$v.adj <- out.orth$v
out.orth$v.adj[,(signs>0.5)] <- -out.orth$v.adj[,(signs>0.5)]

# Export processed data
scores <- scale(tGenes[,match(out.orth$genenames, colnames(tGenes))]) %*% (out.orth$v.adj)
saveRDS(scores, "Allcells_sparsePCA_scores.rds")
saveRDS(out.orth, "seq_well_sparsePCA.rds")
out.orth<-readRDS("seq_well_sparsePCA.rds")

# Plot results
par(mfrow = c(9,3), mar=c(1.5,1.5,1.5,1.5))
for (i in 1:25){
  features <- which(out.orth$v.adj[,i]!=0)
  features <- features[order(abs(out.orth$v.adj[features,i]))]
  barX <- barplot(out.orth$v.adj[features,i], horiz=TRUE, xaxt='n',
                  space=0.1, border=NA, col="gray80")
  xlimz <- par("usr")[1:2]
  text(x = (xlimz[1] + xlimz[2])/2, y = barX, font=2,
       labels=out.orth$genenames[features], cex=1.5)
  abline(v=0)
  title(paste("Module ", i, sep=''), cex.main=2, font.main=3)
}

par(mfrow = c(9,3), mar=c(1.5,1.5,1.5,1.5))
for (i in 26:50){
  features <- which(out.orth$v.adj[,i]!=0)
  features <- features[order(abs(out.orth$v.adj[features,i]))]
  barX <- barplot(out.orth$v.adj[features,i], horiz=TRUE, xaxt='n',
                  space=0.1, border=NA, col="gray80")
  xlimz <- par("usr")[1:2]
  text(x = (xlimz[1] + xlimz[2])/2, y = barX, font=2,
       labels=out.orth$genenames[features], cex=1.5)
  abline(v=0)
  title(paste("Module ", i, sep=''), cex.main=2, font.main=3)
}

scorematrix = out.orth$u
rownames(scorematrix) = rownames(tGenes1)
colnames(scorematrix) = paste0('MODULE', 1:50)


scorematrix = scorematrix[rownames(integrated_seurat@meta.data),]
integrated_seurat@meta.data = cbind(integrated_seurat@meta.data, scorematrix)

# Defining the information in the seurat object of interest
columns <- c(paste0("MODULE", 1:50),
            "UMAP_1", "UMAP_2")

metaMOD <-integrated_seurat@meta.data

# Extracting this data from the seurat object
module_data <- FetchData(integrated_seurat, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(integrated_seurat, vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
                group_by(ident) %>%
                summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
#options(repr.plot.res = 175, repr.plot.height=20, repr.plot.width=20)
#par(mfrow = c(9,3), mar=c(1.5,1.5,1.5,1.5))
map(paste0("MODULE", 1:25), function(module){
        ggplot(module_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=module), 
                           alpha = 0.7) +
                scale_color_gradient(guide = "none", 
                                     low = "grey90", 
                                     high = "blue") +
                ggtitle(module)
}) %>% 
        plot_grid(plotlist = ., nrow=13)

# Plotting a UMAP plot for each of the PCs
#options(repr.plot.res = 175, repr.plot.height=20, repr.plot.width=20)
#par(mfrow = c(9,3), mar=c(1.5,1.5,1.5,1.5))
map(paste0("MODULE", 26:50), function(module){
        ggplot(module_data, 
               aes(UMAP_1, UMAP_2)) +
                geom_point(aes_string(color=module), 
                           alpha = 0.7) +
                scale_color_gradient(guide = "none", 
                                     low = "grey90", 
                                     high = "blue") +
                ggtitle(module)
}) %>% 
        plot_grid(plotlist = ., nrow=13)

## GATA3 cells
median_17 <- median(integrated_seurat@meta.data$MODULE17)
sd_17 <- sd(integrated_seurat@meta.data$MODULE17)

mod_17_fraction <- integrated_seurat@meta.data %>%
                            group_by(orig.ident) %>%
                            summarize(
                              total_rows = n(),
                              condition_rows = sum(MODULE17 > (median_17+(2*sd_17))),
                              fraction_module_17 = condition_rows / total_rows
                            )

## Tregs
median_2 <- median(integrated_seurat@meta.data$MODULE2)
sd_2 <- sd(integrated_seurat@meta.data$MODULE2)

mod_2_fraction <- integrated_seurat@meta.data %>%
                            group_by(orig.ident) %>%
                            summarize(
                              total_rows = n(),
                              condition_rows = sum(MODULE2 > (median_2+(2*sd_2))),
                              fraction_module_2 = condition_rows / total_rows
                            )

## Cytotoxic T cells
median_10 <- median(integrated_seurat@meta.data$MODULE10)
sd_10 <- sd(integrated_seurat@meta.data$MODULE10)

mod_10_fraction <- integrated_seurat@meta.data %>%
                            group_by(orig.ident) %>%
                            summarize(
                              total_rows = n(),
                              condition_rows = sum(MODULE10 > (median_10+(2*sd_10))),
                              fraction_module_10 = condition_rows / total_rows
                            )

results <- list()

# Loop through columns
for (col_name in paste0('MODULE', 1:50)) {
  column_data <- integrated_seurat@meta.data[[col_name]]
  median_val <- median(column_data)
  sd_val <- sd(column_data)
  
  mod_fraction <- integrated_seurat@meta.data %>%
    group_by(orig.ident) %>%
    summarize(
      total_rows = n(),
      condition_rows = sum(.data[[col_name]] > (median_val + (2 * sd_val))),
      fraction = condition_rows / total_rows
    )
  
  results[[col_name]] <- mod_fraction
}


for (tbl_name in names(results)) {
  cat("Tibble:", tbl_name, "\n")
  print(results[[tbl_name]])
  cat("\n")
}

new_list <- list()

for (module_name in names(results)) {
  
  module_data <- results[[module_name]]
  
  bp_fraction <- module_data[module_data$orig.ident == 'BP',]$fraction
  Allergen_fraction <- module_data[module_data$orig.ident == 'Allergen',]$fraction
  
  if (bp_fraction > Allergen_fraction) {
    new_list[[module_name]] <- module_data
  }
}

## Fraction of cell per cluster
fraction_clusters <- list()

# Loop through columns
for (i in unique(integrated_seurat@meta.data$seurat_clusters)) {
  
  fraction <- integrated_seurat@meta.data %>%
    group_by(orig.ident) %>%
    summarize(
      total_rows = n(),
      condition_rows = sum(seurat_clusters == i),
      fraction = condition_rows / total_rows
    )
  
  fraction_clusters[[i]] <- fraction
}

#Pathway Analysis

pathways <- msigdbr("Homo sapiens", "C7") %>%
   format_pathways()


scpa_out_BP_M <- compare_seurat(integrated_seurat,
                           group1 = "orig.ident", 
                           group1_population = c("BP", "Medium"),
                           pathways = hkr_sets)

scpa_out_CD3_M <- compare_seurat(integrated_seurat,
                           group1 = "orig.ident", 
                           group1_population = c("Medium","CD3_CD28"),
                           pathways = pathways)

scpa_out_Allergen_M <- compare_seurat(integrated_seurat,
                           group1 = "orig.ident", 
                           group1_population = c("Medium","Allergen"),
                           pathways = pathways)

scpa_out_BP_Allergen <- compare_seurat(integrated_seurat,
                           group1 = "orig.ident", 
                           group1_population = c("BP","Allergen"),
                           pathways = pathways)

pathways_comb <- c("hallmark", "kegg", "reactome")
hkr_sets <- msigdbr("Homo sapiens") %>%
  filter(grepl(paste(pathways_comb, collapse = "|"), gs_name, ignore.case = T)) %>%
  format_pathways()


scpa_out_Med_Stim <- compare_seurat(integrated_seurat,
                            group1 = "orig.ident", 
                            group1_population = c(c("Allergen", "BP", "CD3_CD28"),"Medium"),
                            pathways = hkr_sets)

scpa_out_BP <- compare_seurat(integrated_seurat,
                           group1 = "orig.ident", 
                           group1_population = c("BP", "Medium"),
                           pathways = hkr_sets)

#Plot the pathway rank
#Highlight apoptosis pathway
p1 <- plot_rank(scpa_out_BP_M, "Apoptosis", 
                highlight_point_size = 3.5, highlight_point_color = "#60c5f7")

#Highlight IL13 related pathways
p2 <- plot_rank(scpa_out_BP_Allergen, "IL13",
                highlight_point_size = 3.5, highlight_point_color = "#fa815c")

patchwork::wrap_plots(p1, p2)


genes <- hkr_sets[sapply(hkr_sets, function(x) any(c("KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY", "REACTOME_TCR_SIGNALING") %in% x$Pathway))]

overlap_KEGG_REACTOME <- unique(genes[[1]]$Genes[genes[[1]]$Genes %in% genes[[2]]$Genes])

features <- list(
  "pathway" = as.character(overlap_KEGG_REACTOME)
)

DotPlot(object = integrated_seurat, features=features, assay = "RNA",
   col.min = (-1), col.max = 1.5, cols = c("gray", "darkviolet"), group.by = "orig.ident") + 
        theme(axis.text.x = element_text(angle = 90)) + RotatedAxis()


#Calculate Module score of a pathway to identify cells that express those genes
integrated_seurat <- AddModuleScore(integrated_seurat,
               features = list(overlap_KEGG_REACTOME), 
               assay="RNA",
               name="TCR")

FeaturePlot(integrated_seurat,
            features = "TCR1", label = TRUE, repel = TRUE, pt.size=1) +
            scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) + 
            labs(title = "TCR Module score")

TCR_Fraction <- integrated_seurat@meta.data %>% group_by(orig.ident) %>%
  summarize(TCR=sum(TCR1 > 0.1), ntotal=n(), fraction= TCR/ntotal, TCR_median= median(TCR1[TCR1 > 0.1]))

print(TCR_Fraction)


# Test second method: Perform GSVA on Seurat clusters
gsva_result <- analyse_sc_clusters(integrated_seurat, verbose = TRUE)

# Extract pathway expression data from the GSVA result
pathway_expression <- pathways(gsva_result)

# simplify the column names by removing the default dataset identifier
colnames(pathway_expression) <- gsub("\\.Seurat", "", colnames(pathway_expression))

# Find the maximum differentially expressed pathway
max_difference <- do.call(rbind, apply(pathway_expression, 1, function(row) {
    values <- as.numeric(row[2:length(row)])
    return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max_difference$diff <- max_difference$max - max_difference$min

# sort based on the difference
max_difference <- max_difference[order(max_difference$diff, decreasing = T), ]

# Plot heatmap of the pathway results
plot_gsva_heatmap(gsva_result, max_pathways = 100, margins = c(6,10))