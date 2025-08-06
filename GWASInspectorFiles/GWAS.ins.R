#!/usr/bin/Rscript

library(data.table)
library(glue)
library(GWASinspector)
job <- setup_inspector('/project/damrauer_scratch/Users/saleemsa/PAD_Signals/GWASInspectorFiles/config.meta-gr38.ini')
job
job <- run_inspector(job)
result_inspector(job)




