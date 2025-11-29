# ==============================================================================
# TiSMeD: Methylation and Expression Analysis Pipeline
# This file now ONLY defines helper functions. It does NOT run the pipeline
# automatically when sourced.
# ==============================================================================

# 1. Load core functions
source("functions.R", encoding = "UTF-8")

# ----------------------------------------------------------------------
# Helper: auto-detect tissue number from a directory of CSVs
# ----------------------------------------------------------------------
get_tissue_number_from_dir <- function(data_dir, is_methylation = TRUE) {
  filelist <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(filelist) == 0) {
    stop(paste("No CSV files found in", data_dir))
  }
  nameall <- data.table()
  for (file_name in filelist) {
    tmp <- fread(file_name)
    # methylation/expression: first column is ID, others are tissues
    tmp <- tmp[, -1, with = FALSE]
    tmpname <- names(tmp) %>% as.data.table()
    nameall <- rbind(nameall, tmpname)
  }
  uniqueN(nameall$.)
}

# ----------------------------------------------------------------------
# PART A: Methylation pipeline (TSHM / TSUM / HKHM / HKUM)
# ----------------------------------------------------------------------
run_methylation_pipeline <- function(
    base_dir   = getwd(),
    spm_cut    = 0.8,
    run_tshm   = TRUE,
    run_tsum   = TRUE,
    run_hkhm   = TRUE,
    run_hkum   = TRUE,
    meth_data_dir = file.path(base_dir, "data/methylome"),
    spm_jar    = file.path(base_dir, "SPM.jar")
) {
  if (!dir.exists(meth_data_dir)) {
    stop(paste("Methylation data directory does not exist:", meth_data_dir))
  }
  message("--- Starting Methylation Analysis ---")
  message("Calculating total tissue numbers for Methylation...")
  
  meth_tissue_number <- get_tissue_number_from_dir(meth_data_dir, is_methylation = TRUE)
  message(paste("Methylation unique tissues:", meth_tissue_number))
  
  # Output dirs (default layout)
  meth_out_tshm <- file.path(base_dir, "results/tshm")
  meth_out_tsum <- file.path(base_dir, "results/tsum")
  meth_out_hkhm <- file.path(base_dir, "results/hkhm")
  meth_out_hkum <- file.path(base_dir, "results/hkum")
  
  if (run_tshm) {
    identify_tshm(
      input_dir     = meth_data_dir,
      output_dir    = meth_out_tshm,
      tissue_number = meth_tissue_number,
      spm_cut       = spm_cut,
      spm_jar_path  = spm_jar
    )
  }
  
  if (run_tsum) {
    identify_tsum(
      input_dir     = meth_data_dir,
      output_dir    = meth_out_tsum,
      tissue_number = meth_tissue_number,
      spm_cut       = spm_cut,
      spm_jar_path  = spm_jar
    )
  }
  
  if (run_hkhm) {
    identify_hkhm(
      input_dir  = meth_data_dir,
      output_dir = meth_out_hkhm
    )
  }
  
  if (run_hkum) {
    identify_hkum(
      input_dir  = meth_data_dir,
      output_dir = meth_out_hkum
    )
  }
  
  message("Methylation pipeline finished.")
}

# ----------------------------------------------------------------------
# PART B: Expression / Proteome pipeline (TSG/TSP)
# ----------------------------------------------------------------------
run_expression_pipeline <- function(
    base_dir      = getwd(),
    spm_cut       = 0.8,
    trans_data_dir = file.path(base_dir, "data/transcriptome"),
    spm_jar       = file.path(base_dir, "SPM.jar")
) {
  if (!dir.exists(trans_data_dir)) {
    stop(paste("Transcriptome/Proteome data directory does not exist:", trans_data_dir))
  }
  message("--- Starting Transcriptome/Proteome Analysis ---")
  message("Calculating total tissue numbers for Transcriptome/Proteome...")
  
  trans_tissue_number <- get_tissue_number_from_dir(trans_data_dir, is_methylation = FALSE)
  message(paste("Transcriptome unique tissues:", trans_tissue_number))
  
  trans_out_dir <- file.path(base_dir, "results/tsg")
  
  identify_tsg_tsp(
    input_dir     = trans_data_dir,
    output_dir    = trans_out_dir,
    tissue_number = trans_tissue_number,
    spm_cut       = spm_cut,
    spm_jar_path  = spm_jar
  )
  
  message("Transcriptome/Proteome pipeline finished.")
}
