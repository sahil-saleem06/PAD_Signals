library(tidyverse)
library(GenomicRanges)

gwas_files <- list(
  AFR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_afr.txt.gz",
  EUR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_eur.txt.gz",
  EAS = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_eas.txt.gz",
  AMR = "/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_amr.txt.gz"
)

#Loops through 4 regions
for (region in names(gwas_files)) {
  
  gwas_file_path <- gwas_files[[region]]
  
  temp_raw <- vroom::vroom(gwas_file_path)
  
  temp_tidy <- temp_raw %>% dplyr::select(-c(SNP, BUILD))
  temp_raw <- na.omit(temp_raw)
  
  output_file <- paste0("/project/damrauer_scratch/Users/saleemsa/PAD_Signals/cleaned_sumstats_", region, ".txt")
  
  vroom::vroom_write(temp_raw, output_file)
}

