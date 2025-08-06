library(dplyr)
library(vroom)

# Load raw summary statistics (already in GRCh38)
pad_raw <- fread("/project/voltron/Resources/GWAS-Summary-Statistics/MVP/PAD5c/Apr2025/20250411_mvp_pad5c_eur.txt.gz")

# Fix CHR column
pad_raw[, CHR := gsub("^chr", "", CHR)]
pad_raw[, CHR := as.numeric(CHR)]

# Sort efficiently
pad_sorted <- setorder(pad_raw, CHR, POS)


# Drop any unnecessary columns (example: if 'BUILD' or 'SNP' exists, remove it)
pad_tidy <- pad_sorted %>%
  select(-any_of(c("BUILD", "SNP"))) %>%  # safe way to remove if they exist
  na.omit()                               # drop any rows with NAs

pad_ready <- pad_tidy %>%
  dplyr::mutate(N = (N_case+N_ctrl)) %>%
  select(CHROM = CHR, POS, MARKER = MarkerID, EFFECT_ALLELE = Allele1, OTHER_ALLELE = Allele2, EFF_ALL_FREQ = AF_Allele2, BETA, SE, PVAL = p.value, INFO = imputationInfo, N, N_CASE = N_case, N_CONTROL = N_ctrl)

# Write out the cleaned file
vroom_write(pad_ready, "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/GWASInspectorInput/QC_ready_PAD5c_eur.txt", delim = "\t")
