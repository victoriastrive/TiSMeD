# ==============================================================================
# TiSMeD: Methylation and Expression Analysis Pipeline
# Usage: Run this script to execute tissue-specific identification.
# ==============================================================================

# 1. Load configuration and functions
source("functions.R", encoding = "UTF-8")

# 2. Define Directories (USER CONFIGURATION)
BASE_DIR <- getwd()
SPM_JAR  <- file.path(BASE_DIR, "SPM.jar") 

# Methylation Directories
METH_DATA_DIR <- file.path(BASE_DIR, "data/methylome")
METH_OUT_TSHM <- file.path(BASE_DIR, "results/tshm")
METH_OUT_TSUM <- file.path(BASE_DIR, "results/tsum")
METH_OUT_HKHM <- file.path(BASE_DIR, "results/hkhm")
METH_OUT_HKUM <- file.path(BASE_DIR, "results/hkum")

# Transcriptome/Proteome Directories
TRANS_DATA_DIR <- file.path(BASE_DIR, "data/transcriptome")
TRANS_OUT_DIR  <- file.path(BASE_DIR, "results/tsg")


# ==============================================================================
# PART A: Methylation Analysis
# ==============================================================================

if (dir.exists(METH_DATA_DIR)) {
  message("--- Starting Methylation Analysis ---")
  
  # 3. Calculate Tissue Number Dynamically (Methylation)
  message("Calculating total tissue numbers for Methylation...")
  filelist <- list.files(METH_DATA_DIR, pattern = ".csv", full.names = TRUE)
  nameall <- data.table()
  if(length(filelist) > 0) {
    for (file_name in filelist) {
      tmp <- fread(file_name)
      tmp <- tmp[, -1, with = FALSE] 
      tmpname <- names(tmp) %>% as.data.table()
      nameall <- rbind(nameall, tmpname)
    }
    meth_tissue_number <- uniqueN(nameall$.)
    message(paste("Methylation unique tissues:", meth_tissue_number))
    
    # 4. Execute Methylation Modules
    # Note: Function names changed to 'identify_' as requested
    
    identify_tshm(METH_DATA_DIR, METH_OUT_TSHM, meth_tissue_number, 0.8, SPM_JAR)
    identify_tsum(METH_DATA_DIR, METH_OUT_TSUM, meth_tissue_number, 0.8, SPM_JAR)
    identify_hkhm(METH_DATA_DIR, METH_OUT_HKHM)
    identify_hkum(METH_DATA_DIR, METH_OUT_HKUM)
    
  } else {
    message("No CSV files found in Methylation data directory.")
  }
}

# ==============================================================================
# PART B: Transcriptome / Proteome Analysis
# ==============================================================================

if (dir.exists(TRANS_DATA_DIR)) {
  message("--- Starting Transcriptome/Proteome Analysis ---")
  
  # 5. Calculate Tissue Number Dynamically (Transcriptome) [cite: 3]
  message("Calculating total tissue numbers for Transcriptome...")
  filelist <- list.files(TRANS_DATA_DIR, pattern = ".csv", full.names = TRUE)
  nameall <- data.table()
  
  if(length(filelist) > 0) {
    for (file_name in filelist) {
      tmp <- fread(file_name)
      # Remove Gene Symbol column (usually 1st) to count tissues
      tmp <- tmp[, -1, with = FALSE] 
      tmpname <- names(tmp) %>% as.data.table()
      nameall <- rbind(nameall, tmpname)
    }
    trans_tissue_number <- uniqueN(nameall$.)
    message(paste("Transcriptome unique tissues:", trans_tissue_number))
    
    # 6. Execute TSG/TSP Identification
    # NOTE: If your count file is named 'count_mod.txt', change the parameter below.
    identify_tsg_tsp(
      input_dir = TRANS_DATA_DIR, 
      output_dir = TRANS_OUT_DIR, 
      tissue_number = trans_tissue_number, 
      spm_cut = 0.8, 
      spm_jar_path = SPM_JAR
    )
    
  } else {
    message("No CSV files found in Transcriptome data directory.")
  }
}

message("Pipeline completed.")