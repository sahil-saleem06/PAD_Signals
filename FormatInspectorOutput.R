library(data.table)

# Output folder for cleaned versions
output_dir <- "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/MetalInput2/"
dir.create(output_dir, showWarnings = FALSE)

# List of specific files you want to format
input_files <- c(
  "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/MetalInput//QC_QC_ready_PAD5c_afr.txt.gz",
  "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/MetalInput/QC_QC_ready_PAD5c_amr.txt.gz",
  "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/MetalInput/QC_QC_ready_PAD5c_eas.txt.gz",
  "/project/damrauer_scratch/Users/saleemsa/PAD_Signals/MetalInput/QC_QC_ready_PAD5c_eur.txt.gz"
)

for (file in input_files) {
  cat("Processing:", file, "\n")
  
  # Read file
  dt <- fread(file)
  
  # Remove MARKER column if it exists
  if ("highDiffEAF" %in% names(dt)) {
    dt[, highDiffEAF := NULL]
  }
  
  # Write cleaned file to new folder with same name
  out_path <- file.path(output_dir, basename(file))
  fwrite(dt, sub(".gz$", "", out_path), sep = "\t")
  system(paste("gzip -f", shQuote(sub(".gz$", "", out_path))))
}
