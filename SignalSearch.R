require('tidyverse')
library(data.table)

#Set input GWAS files
gwas_files <- list(
  AFR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_afr.txt.gz",
  EUR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_eur.txt.gz",
  EAS = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_eas.txt.gz",
  AMR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_amr.txt.gz"
)


#Get Gene Coords
library(httr)
library(jsonlite)
library(dplyr)
library(biomaRt)

#Jagged = JAG1 and JAG2, COX2 = PTGS2, IL17 = IL17A
gene_list <- c("NSD2", "KMT2A", "KDM5B", "KDM6B", "SETDB2", "SETD4", 
               "HDAC11", "HDAC3", "PRMT1", "KDM5C", "STAT3", "DLL4", 
               "NOTCH1", "NOTCH2", "JAG1", "JAG2", "PTGS2", "IL17A", "TGFBR2")

#Get coordinates of each gene
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

gene_coords <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", 
                 "chromosome_name", "start_position", "end_position"),
  filters = "hgnc_symbol",
  values = gene_list,
  mart = ensembl
)

#match GWAS format and add buffer
gene_coords <- gene_coords %>%
  mutate(chromosome_name = paste0("chr", chromosome_name),
         region_start = pmax(1, start_position - 100000),
         region_end = end_position + 100000)


#Function to take in gene data and output corresponding GWAS variants
extract_gene_region <- function(gene, chr, start, end, gwas_file_path, colnames) {
  df <- fread(cmd = paste0("zcat ", gwas_file_path,
                           " | awk -F '\t' '$1==\"", chr, 
                           "\" && $2>=", start, " && $2<=", end, "'"),
              col.names = colnames)
  if (nrow(df) > 0) {
    df$gene <- gene
  }
  return(df)
}

#Function to do Manhattan plot of single search
plot_gene_signal <- function(df, gene_label) {
  ggplot(df, aes(x = POS, y = -log10(p.value))) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", size = 0.8) +
    labs(
      title = paste("Signal near", gene_label),
      x = "Genomic Position",
      y = "-log10(p-value)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid = element_line(color = "gray90")
    )
}

#For all regions
library(patchwork)

all_plots <- list()
base_output_dir <- "/project/damrauer_scratch/Users/saleemsa/MultiRegion_Gene_Signals/"

#Loop through 4 regions
for (region in names(gwas_files)) {
  gwas_file_path <- gwas_files[[region]]
  
  # Create region-specific output directory
  region_output_dir <- file.path(base_output_dir, region)
  dir.create(region_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load column names once per file
  header_line <- system(paste("zcat", gwas_file_path, "| head -n 1"), intern = TRUE)
  colnames <- strsplit(header_line, "\t")[[1]]
  
  for (i in 1:nrow(gene_coords)) {
    gene <- gene_coords$hgnc_symbol[i]
    chr <- gene_coords$chromosome_name[i]
    start <- gene_coords$region_start[i]
    end <- gene_coords$region_end[i]
    
    #Call function that gets GWAS info
    df <- extract_gene_region(gene, chr, start, end, gwas_file_path, colnames)
    
    #Save plot to list
    if (nrow(df) > 0) {
      gene_label <- paste(gene, region, sep = "_")
      p <- plot_gene_signal(df, gene_label)
      all_plots[[gene]][[region]] <- p
      
      # Save individual plot
      ggsave(
        filename = file.path(region_output_dir, paste0(gene_label, "_plot.png")),
        plot = p,
        width = 8, height = 5, dpi = 300
      )
    } else {
      message(paste("No data for", gene, "in", region))
    }
  }
}


#Combined plots
for (gene in names(all_plots)) {
  if (length(all_plots[[gene]]) > 0) {
    combined <- wrap_plots(all_plots[[gene]], ncol = 2) +
      plot_annotation(title = paste("Gene:", gene))
    
    ggsave(
      filename = file.path(base_output_dir, paste0(gene, "_combined_plot.png")),
      plot = combined,
      width = 16, height = 10, dpi = 300
    )
  }
}

#Meta-Analysis

