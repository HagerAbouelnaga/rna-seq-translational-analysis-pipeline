# Title: DESeq2 RNA-seq Pipeline
# Author: Hager Salah Abouelnaga
# Description: End-to-end differential expression analysis

#------------------------- 

# Load libraries 

#------------------------- 

library(DESeq2)  

library(tidyverse)  

library(AnnotationDbi)  

library(org.Hs.eg.db)  

library(apeglm)  

library(EnhancedVolcano)  

library(ComplexHeatmap)  

library(openxlsx)  

# for write.xlsx library(ggrepel) 

#------------------------- 

# Config / paths (EDIT) 

#------------------------- 

#Use relative paths in your repo, e.g. "data/counts/" or "resources/coldata.csv" 

counts_dir <- "path/to/Input_RNA" # folder containing sample_count-files (*.txt)  

coldata_file <- "path/to/resources/coldata.csv" # CSV with sample metadata (sample_name, Sample, ...)  

output_dir <- "path/to/results" # output folder to save tables and plots 

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE) 

#------------------------- 

# Read count files and merge 

#-------------------------  

count_files <- list.files(path = counts_dir, pattern = "\\.txt$", full.names = TRUE) 
if (length(count_files) == 0) {
    stop("No count files found. Check counts_dir path.") 

     }

 

merged_data <- NULL  

for (f in count_files) { 

#get sample name from filename (customize if your naming differs) 

sample_name <- basename(f) %>%  

str_remove_all("\\.txt$") %>%  

str_replace_all("-counts", "")  

tmp <- read_tsv(f, comment = "#", col_names = TRUE) 

 

#pick columns: first column gene id (ENSEMBL) and a counts column. 

#modify below if your count file columns differ. 

#Typical featureCounts output: GeneID, Chr, Start, ... , sample.counts 

#Here we assume the counts are in column 2 or named like sample. 

#Try to detect counts column: 

if (ncol(tmp) < 2) stop(paste("Unexpected file format:", f)) 

 

#assume first column is gene id and last column the counts 

tmp_sel <- tmp %>% select(1, ncol(tmp))  

colnames(tmp_sel) <- c("GeneID", sample_name) 

if (is.null(merged_data)) { 

 merged_data <- tmp_sel  

} else {  

merged_data <- full_join(merged_data, tmp_sel, by = "GeneID")  

 }  

} 

 

#convert to data.frame and set rownames 

merged_df <- as.data.frame(merged_data)  

rownames(merged_df) <- merged_df$GeneID  

merged_df$GeneID <- NULL 

#replace NAs with 0 (if you prefer) 

merged_df[is.na(merged_df)] <- 0 

#------------------------- 

# Add gene SYMBOL column (optional, mapping ENSEMBL -> SYMBOL) 

#------------------------- 

#Only do this if your rownames are ENSG IDs. If they are gene symbols, skip. 

ensembl_ids <- rownames(merged_df) 

symbol_map <- AnnotationDbi::mapIds(org.Hs.eg.db,  

keys = ensembl_ids,  

column = "SYMBOL",  

keytype = "ENSEMBL",  

multiVals = "first") 

#Add as a column if you like (not required by DESeq2) 

merged_df <- merged_df %>%  

rownames_to_column(var = "ENSEMBL") %>%  

mutate(SYMBOL = symbol_map[ENSEMBL]) %>% 

 column_to_rownames(var = "ENSEMBL") 

#------------------------- 

# Read coldata (sample metadata) 

#-------------------------  

coldata <- read.csv(coldata_file, stringsAsFactors = FALSE) 

#Ensure coldata has a column 'sample_name' or adjust accordingly 

#Set rownames to sample names used above 

if (!"sample_name" %in% colnames(coldata)) { stop("coldata.csv must contain a column named 'sample_name' that matches count file names.")  

}  

rownames(coldata) <- coldata$sample_name 


# Check matching between coldata and count matrix

if (!all(rownames(coldata) %in% colnames(merged_df))) {
  stop("Mismatch between count matrix and coldata samples")
}


#Ensure ordering of columns in count matrix matches coldata rows 

count_matrix <- merged_df %>% dplyr::select(all_of(rownames(coldata))) 

#------------------------- 

# Create DESeq2 object 

#------------------------- 

#Factorize condition column - change 'Sample' to your condition column name 

coldata$Sample <- factor(coldata$Sample) 

#set reference level if needed (example: "control1") 

coldata$Sample <- relevel(coldata$Sample, ref = "control1") 

dds <- DESeqDataSetFromMatrix(countData = count_matrix,  

colData = coldata,  

design = ~ Sample) 

#filter low count genes 

keep <- rowSums(counts(dds)) > 1  

dds <- dds[keep, ] 

#------------------------- 

# Run DESeq2 

#-------------------------  

dds <- DESeq(dds) 

 

#optional: run vst for PCA 

vsd <- vst(dds, blind = FALSE)  

pca_plot <- plotPCA(vsd, intgroup = "Sample")  

ggsave(file.path(output_dir, "PCA_plot.png"), plot = pca_plot, width = 6, height = 5) 

 

#------------------------- 

# List results names and perform LFC shrinkage for example contrast 

#-------------------------  

message("Available results names (contrasts):")  

print(resultsNames(dds)) 

 

#Example: replace with the exact coef string you need, 

##e.g. "Sample_G1.8nM_vs_control1" as in your original script. 

contrast_name <- "Sample_G1.8nM_vs_control1"  

if (!(contrast_name %in% resultsNames(dds))) { 
message("Contrast name not found in resultsNames(dds). Use resultsNames(dds) to choose correct coef.")  

} else {  

res_lfc <- lfcShrink(dds, coef = contrast_name, type = "apeglm")  

res_df <- as.data.frame(res_lfc) %>% 

 rownames_to_column(var = "ENSEMBL")  

res_df$SYMBOL <- AnnotationDbi::mapIds(org.Hs.eg.db,  

keys = res_df$ENSEMBL, 

 column = "SYMBOL",  

keytype = "ENSEMBL",  

multiVals = "first") 

 

#save table 

write.xlsx(res_df, file.path(output_dir, paste0("DE_results_", contrast_name, ".xlsx")), overwrite = TRUE) 
} 

 

#------------------------- 

# Volcano plot  

#-------------------------  

 

if (exists("res_df")) {  

EnhancedVolcano(res_df,  

lab = res_df$SYMBOL,  

x = 'log2FoldChange',  

y = 'pvalue',  

pCutoff = 0.05,  

FCcutoff = 1.5,  

title = paste0(contrast_name, " volcano"))  

ggsave(file.path(output_dir, paste0("Volcano_", contrast_name, ".png")), width = 7, height = 6)  

} 

 

#------------------------- 

# Optionally produce heatmap of top genes 

#------------------------- 

 

#Example: top 50 by padj 

if (exists("res_df")) {  

sig <- res_df %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(50) 

 

#extract normalized counts for these genes 

norm_counts <- assay(vsd)[sig$ENSEMBL, , drop = FALSE]  

rownames(norm_counts) <- sig$SYMBOL 

 

#scale rows 

mat <- t(scale(t(norm_counts)))  

Heatmap(mat, name = "zscore")  

} 

 

message("DESeq2 pipeline completed. Check output dir: ", output_dir) 