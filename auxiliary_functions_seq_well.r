# functions.R


get_gene_symbol_rename_duplicate <- function(data) {

    #' Get Gene Symbols and Rename Duplicates in Seurat Object
    #'
    #' This function takes a Seurat object containing count matrix data, retrieves
    #' gene symbols using Ensembl gene IDs, and handles duplicate gene symbols by
    #' adding suffixes to make them unique. The updated gene symbols are used as
    #' row names in the Seurat object.
    #'
    #' @param data A Seurat object containing count matrix data.
    #'
    #' @return A modified Seurat object with updated gene symbols as row names.
    #'
    #' @details This function connects to the Ensembl BioMart to retrieve gene symbols
    #'          corresponding to the Ensembl gene IDs present in the input Seurat object.
    #'          It handles missing or empty gene symbols by replacing them with the
    #'          corresponding Ensembl gene IDs. Duplicate gene symbols are given unique
    #'          suffixes to ensure uniqueness.
    #'
    #' @examples
    #' \dontrun{
    #' # Create a Seurat object using count matrix data
    #' seurat_obj <- CreateSeuratObject(counts = count_matrix)
    #'
    #' # Get updated Seurat object with modified gene symbols
    #' seurat_obj_updated <- get_gene_symbol_rename_duplicate(seurat_obj)
    #' }
    #'
    #' @import Seurat
    #' @importFrom biomaRt
    #' @export

    # Get the Ensembl gene IDs from the row names of the merged count matrix
    ensembl_gene_ids <- rownames(data)
    
    # Connect to the Ensembl BioMart
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Query the gene symbols using the Ensembl gene IDs
    gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), filters = "ensembl_gene_id", values = ensembl_gene_ids, mart = ensembl)
    
    # Create a mapping from Ensembl gene IDs to gene symbols
    ensembl_to_symbol <- setNames(gene_info$hgnc_symbol, gene_info$ensembl_gene_id)
    
    # Replace missing or empty values in hgnc_symbol with Ensembl gene IDs
    missing_indices <- which(is.na(ensembl_to_symbol) | ensembl_to_symbol == "")
    ensembl_to_symbol[missing_indices] <- names(ensembl_to_symbol)[missing_indices]
    
    # Replace Ensembl gene IDs with gene symbols in the merged Seurat object's RNA assay
    new_row_names <- ensembl_to_symbol[rownames(data)]
    
    # Handle duplicate gene symbols
    duplicated_symbols <- new_row_names[duplicated(new_row_names)]
    
    for (duplicate in duplicated_symbols) {
        duplicate_indices <- which(new_row_names == duplicate)
        
        # Reset the suffix counter to 1 for each set of duplicates
        suffix_counter <- 1
        
        for (index in duplicate_indices) {
            new_suffix <- paste0("_", suffix_counter)
            new_row_names[index] <- paste0(duplicate, new_suffix)
            suffix_counter <- suffix_counter + 1
        }
    }
    
    # Assign updated row names to the data object
    rownames(data) <- new_row_names
    
    return(data)
}


create_seurat_obj <- function(count_matrix_files){

    #' Create Seurat Objects from Count Matrix Files
    #'
    #' This function reads count matrix files, creates Seurat objects for each file,
    #' and returns a list of Seurat objects. Each Seurat object corresponds to a sample's
    #' count matrix data, with the sample name extracted from the filename.
    #'
    #' @param count_matrix_files A character vector of count matrix filenames to be processed.
    #'
    #' @return A list of Seurat objects, where each Seurat object contains the count matrix
    #'         data for a sample. The sample names are used as project names for the Seurat objects.
    #'
    #' @details This function uses regular expressions to extract the sample name from each filename.
    #'          The extracted sample name is used as the project name for the Seurat object.
    #'          The count matrix data from each file is used to create the Seurat object.
    #'
    #' @examples
    #' \dontrun{
    #' # List of count matrix filenames
    #' data <- c(
    #'   "Sample_1",
    #'   "Sample_2"
    #' )
    #' 
    #' # Create Seurat objects using the provided data files
    #' seurat_objects <- create_seurat_obj(data)
    #' }
    #'
    #' @import Seurat
    #' @export

    # Create an empty list to store Seurat objects
    seurat_objs <- list()
    # Loop through each file in the list of count matrix filenames
    for (file in count_matrix_files) {
      # Read the count matrix data from the file
      data <- read.table(file, sep = "\t", header = TRUE, row.names = 1)
  
      data <- get_gene_symbol_rename_duplicate(data)
  
      # Extract sample name from the filename using regular expression
      seurat_ob_name <- sub("^(.*?)_S\\d+_R\\d+_\\d+.*$", "\\1", file)
      
      # Create a Seurat object with the count matrix data
      seurat_ob <- CreateSeuratObject(counts = data, project = seurat_ob_name)
      
      # Add the Seurat object to the list with its name as the key
      seurat_objs[[seurat_ob_name]] <- seurat_ob
  }

  # Return the list of created Seurat objects
  return(seurat_objs)
}


transcription_factors <- c("NFKB1", "JUN", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", 
                            "STAT5B", "STAT6", "IRF1", "IRF2", "IRF3", "IRF4", "IRF5", "IRF6", 
                            "IRF7", "IRF8", "GATA1", "GATA2", "GATA3", "GATA4", "GATA5", "GATA6", 
                            "TBX21", "FOXP1", "FOXP2", "FOXP3", "RORA", "RORB", "RORC", 
                            "MYC", "MYB", "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", 
                            "E2F7", "E2F8", "HIF1A", "HIF1B", "HIF2A", "HIF3A", "SPI1", 
                            "CEBPA", "CEBPB", "CEBPD", "ATF1", "ATF2", "ATF3", "ATF4", "ATF5", 
                            "ATF6", "ATF7", "EGR1", "EGR2", "EGR3", "EGR4", "TCF3", "TCF4", 
                            "PAX5", "EBF1", "IKZF1", "IKZF2", "IKZF3", "ZBTB7A", "ZBTB7B", 
                            "BCL6", "FOXO1", "FOXO3", "FOXO4", "FOXO6", "BACH2", "BATF", "BATF3", 
                            "MAF", "ELF1", "FOS", "FOSB", "FOSL1", "FOSL2", "IRF9", "NFATC1", 
                            "NFATC2", "NFATC3", "NFATC4", "RELA", "RELB", "RXRA", "RXRB", 
                            "RXRG", "SREBF1", "SREBF2", "TCF7", "TGIF1", "TGIF2", "TOX", "REL", "NR4A3")

extra_genes <- c("PTGDR2", "KLRB1", "IL17RB", "TGFBR3", "PLA2G16", "HPGDS", "PPARG", "ACADVL", 
                "ACSL4", "SLC27A2", "LPCAT2", "DGKE", "GK", "CHDH", "ALOX5AP", "BCL6", "SELL", 
                "MAF", "ICOS", "PDCD1", "CD28", "CD40LG", "IL21R", "IL21", "IL2RA","CXCR5", 
                "CCR7", "CXCL13", "REL", "CTLA4", "ICOS", "CXCR5", "NR4A3", "IL2", "NFKB1", 
                "DUSP4", "TNFRSF9", "PDCD1", "TNF", "REL", "CTLA4", "ICOS", "CXCR5", "NR4A3",
                 "IL2", "NFKB1", "DUSP4", "TNFRSF9", "PDCD1", "TNF", "GZMB", "IFNG", "IL4", 
                 "IL5", "IL13","IL9", "IL17A", "IL10", "IL21","RORC", "GATA3", "TBX21", "CCR9", 
                 "CCR7", "CCR6", "CCR5","CXCR6", "CXCR5", "CXCR3", "GPR15", "CD69", "TFRC", 
                 "TNFRSF9", "CD40LG", "CD38", "CD27", "TIGIT", "DPP4", "CD27", "CCL17", 
                 "CCR2", "CCL19", "CLA", "ITGB7")




differential_gene_expression <- function(seurat_object, condition_1, condition_2, variable_to_group_by, genes_to_plot, name_and_path_plot) {
    #' This function performs differential gene expression analysis between two conditions,
    #' identifies significant markers, and generates a volcano plot highlighting these markers.
    #'
    #' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
    #' @param condition_1 The first condition for comparison.
    #' @param condition_2 The second condition for comparison.
    #' @param variable_to_group_by The variable in the metadata to group cells by.
    #' @param genes_to_plot A vector of gene names to include in the plot.
    #' @param name_and_path_plot The file path and name for saving the plot.
    #' @return A data frame of significant markers sorted by log fold change.
    #' @import Seurat
    #' @import ggplot2
    #' @import dplyr

    # Find differentially expressed markers
    markers <- FindMarkers(seurat_object, ident.1 = condition_1, 
                           ident.2 = condition_2, 
                           group.by = variable_to_group_by, 
                           assay = "RNA", logfc.threshold = 0, min.pct = 0)

    # Add gene names to the markers data frame
    markers$gene <- rownames(markers)
    markers <- markers[markers$gene %in% genes_to_plot,]

    # Adjust p-values for plotting
    markers$pval_use <- markers$p_val_adj
    markers$pval_use[markers$pval_use < 1e-60] <- 1e-60

    # Create the volcano plot
    plot <- ggplot(markers, aes(x = avg_log2FC, y = -log10(pval_use))) + 
        geom_point(aes(color = p_val_adj < 0.001 & (avg_log2FC > 0.3 | avg_log2FC < -0.3)), show.legend = FALSE) + 
        geom_vline(xintercept = c(-0.3, 0.3), color = "gray", linetype = "dashed") +
        geom_hline(yintercept = c(-log10(0.001)), color = "gray", linetype = "dashed")

    # Label significant genes
    genes.to.label <- rownames(subset(markers, p_val_adj < 0.001 & (avg_log2FC > 0.3 | avg_log2FC < -0.3)))
    plot <- LabelPoints(plot = plot, points = genes.to.label, repel = TRUE, xnudge = 0, ynudge = 0, max.overlaps = 20)

    # Save the plot
    ggsave(name_and_path_plot, plot)

    # Return significant markers sorted by log fold change
    significant_markers <- markers[markers$p_val_adj < 0.05, ]
    return(significant_markers[order(significant_markers$avg_log2FC),])
}