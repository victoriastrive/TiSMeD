# ==============================================================================
# TiSMeD: Tissue-Specific Methylation and Expression Database
# Core Functions for Identifying Tissue-Specific Markers (Methylation, RNA, Protein)
# ==============================================================================

library(data.table)
library(tidyverse)
library(readxl)
library(tidyr)

# Helper function: Sigmoid transformation
sigmoid <- function(x, a = 1) {
  1 / (1 + exp(-x))
}

# ==============================================================================
# SECTION 1: Methylation Analysis Functions
# ==============================================================================

#' Identify Tissue-Specific hypermethylation sites (TSHMs)
#' @param input_dir Directory containing input dataset expression files (.csv)
#' @param output_dir Directory to save intermediate and final results
#' @param tissue_number Total number of unique tissue types across all datasets
#' @param spm_cut Cutoff value for SPM_adjust value
#' @param spm_jar_path Path to the SPM.jar file

identify_tshm <- function(input_dir, output_dir, tissue_number, spm_cut, spm_jar_path = "SPM.jar") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  count_file <- file.path(input_dir, "count.txt")
  if (!file.exists(count_file)) stop("count.txt not found in input directory.")
  
  count <- fread(count_file)[order(dataset)]
  count$log10sample <- log10(count$sample)
  count$log10tissue <- log10(count$tissue)
  count$sigmoidlog10_tissue <- sapply(count$log10tissue, sigmoid)
  count$sigmoidlog10_tissue_diall <- count$sigmoidlog10_tissue / sigmoid(log10(tissue_number))
  count$sigmoidlog10_sample <- sapply(count$log10sample, sigmoid)
  
  filelist <- list.files(input_dir, pattern = "*.csv")
  TSgeneall_list <- list()
  genetable_list <- list()
  i <- 1
  
  for (file_name in filelist) {
    filepath <- file.path(input_dir, file_name)
    message(paste("Processing (Methylation):", file_name))
    tmp_table <- data.table::fread(filepath)
    names(tmp_table)[1] <- "Site_ID" 
    names(tmp_table) <- gsub(" ", "-", names(tmp_table))
    names(tmp_table) <- gsub(",", "*", names(tmp_table))
    
    temp_csv <- file.path(output_dir, "tmp.csv")
    fwrite(tmp_table, temp_csv)
    
    system(paste("java -jar", spm_jar_path, temp_csv))
    
    spm_output <- paste0(temp_csv, ".spm")
    new_spm_name <- "tmp_spm.csv"
    if(file.exists(spm_output)) file.rename(spm_output, new_spm_name)
    
    test_spm <- fread(new_spm_name)
    if("DPM" %in% names(test_spm)) test_spm[, DPM := NULL]
    test_spm <- pivot_longer(test_spm, !Site_ID, names_to = "tissue", values_to = "spm")
    test_spm$spm <- test_spm$spm * count$sigmoidlog10_tissue_diall[i]
    
    genelist <- test_spm[, 1:2]
    genelist$dataset <- file_name
    genelist$log10sample <- count$sigmoidlog10_sample[i]
    genetable_list[[i]] <- genelist 
    
    TSgene <- test_spm[test_spm$spm >= spm_cut, ] 
    tmp_table_long <- pivot_longer(tmp_table, !Site_ID, names_to = "tissue", values_to = "expre")
    TSgene <- as.data.table(tmp_table_long)[as.data.table(TSgene), on = c("Site_ID", "tissue"), nomatch = 0]
    TSgene$log10sample <- count$sigmoidlog10_sample[i]
    TSgene$dataset <- file_name
    TSgeneall_list[[i]] <- TSgene
    i <- i + 1
  }
  
  genetable <- rbindlist(genetable_list, use.names = TRUE, fill = TRUE)
  TSgeneall <- rbindlist(TSgeneall_list, use.names = TRUE, fill = TRUE)
  
  TSgeneall_summarise <- TSgeneall[, .(
    supportcount = .N,
    tsgsumweight = sum(log10sample, na.rm = TRUE),
    meanspm = mean(spm, na.rm = TRUE)
  ), by = c("Site_ID", "tissue")]
  
  TSgeneall07 <- TSgeneall %>% filter(expre >= 0.7)

  genesumweight <- genetable[, .(
    gene_count = .N,
    genesumweight = sum(log10sample, na.rm = TRUE)
  ), by = c("Site_ID", "tissue")]
  
  tmp <- merge(TSgeneall_summarise, TSgeneall07, by = c("Site_ID", "tissue"), all.y = TRUE)
  tmp <- merge(tmp, genesumweight, by = c("Site_ID", "tissue"), all.x = TRUE)
  
  tmp$Tscore <- tmp$tsgsumweight / tmp$genesumweight
  
  Tscoretable <- tmp %>%
    as.data.table() %>%
    .[, -c("expre", "spm", "log10sample", "dataset","supportcount", "tsgsumweight", "genesumweight")] %>%
    unique()
  Tscoretable$tissue <- gsub("-", " ", Tscoretable$tissue)
  Tscoretable$tissue <- gsub("\\*", ",", Tscoretable$tissue)
  Tscoretable <- Tscoretable %>% filter(gene_count >= 2) %>% filter(Tscore >= 0.6)
  Tscoretable <- Tscoretable[,-c("gene_count")]
  names(Tscoretable)<-c("Site_ID","Tissue","Mean_spm_adjust","Tscore")
  Tscoretable$Mean_spm_adjust<-round(Tscoretable$Mean_spm_adjust,2)
  Tscoretable$Tscore<-round(Tscoretable$Tscore,2)
  fwrite(Tscoretable, file.path(output_dir, "TSHMs.csv"))
  gc()
}

#' Identify Tissue-Specific hypomethylation sites (TSUM)
#' @param input_dir Directory containing input dataset expression files (.csv)
#' @param output_dir Directory to save intermediate and final results
#' @param tissue_number Total number of unique tissue types across all datasets
#' @param spm_cut Cutoff value for SPM_adjust value
#' @param spm_jar_path Path to the SPM.jar file

identify_tsum <- function(input_dir, output_dir, tissue_number, spm_cut, spm_jar_path = "SPM.jar") {
  # (Logic identical to previous calculate_tsum, renamed for consistency)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  count <- fread(file.path(input_dir, "count.txt"))[order(dataset)]
  
  count$log10sample <- log10(count$sample)
  count$log10tissue <- log10(count$tissue)
  count$sigmoidlog10_tissue <- sapply(count$log10tissue, sigmoid)
  count$sigmoidlog10_tissue_diall <- count$sigmoidlog10_tissue / sigmoid(log10(tissue_number))
  count$sigmoidlog10_sample <- sapply(count$log10sample, sigmoid)
  
  filelist <- list.files(input_dir, pattern = "*.csv")
  TSgeneall_list <- list()
  genetable_list <- list()
  i <- 1
  
  for (file_name in filelist) {
    filepath <- file.path(input_dir, file_name)
    message(paste("Processing (Methylation):", file_name))
    tmp_table <- data.table::fread(filepath)
    names(tmp_table)[1] <- "Site_ID"
    names(tmp_table) <- gsub(" ", "-", names(tmp_table))
    names(tmp_table) <- gsub(",", "*", names(tmp_table))
    
    # Invert values
    test <- tmp_table[, 2:ncol(tmp_table)]
    testmax <- (test + 0.0001) / do.call(pmax, test + 0.0001) 
    testcorect <- 1 - testmax
    testcorect <- cbind(tmp_table[, 1], testcorect)
    
    temp_csv <- file.path(output_dir, "tmp.csv")
    fwrite(testcorect, temp_csv)
    
    system(paste("java -jar", spm_jar_path, temp_csv))
    
    spm_output <- paste0(temp_csv, ".spm")
    new_spm_name <- "tmp_spm.csv"
    if(file.exists(spm_output)) file.rename(spm_output, new_spm_name)
    
    test_spm <- fread(new_spm_name)
    if("DPM" %in% names(test_spm)) test_spm[, DPM := NULL]
    test_spm <- pivot_longer(test_spm, !Site_ID, names_to = "tissue", values_to = "spm")
    test_spm$spm <- test_spm$spm * count$sigmoidlog10_tissue_diall[i]
    
    genelist <- test_spm[, 1:2]
    genelist$dataset <- file_name
    genelist$log10sample <- count$sigmoidlog10_sample[i]
    genetable_list[[i]] <- genelist
    
    TSgene <- test_spm[test_spm$spm >= spm_cut, ]
    tmp_table_long <- pivot_longer(tmp_table, !Site_ID, names_to = "tissue", values_to = "expre")
    TSgene <- as.data.table(tmp_table_long)[as.data.table(TSgene), on = c("Site_ID", "tissue"), nomatch = 0]
    TSgene$log10sample <- count$sigmoidlog10_sample[i]
    TSgene$dataset <- file_name
    TSgeneall_list[[i]] <- TSgene
    i <- i + 1
  }
  
  genetable <- rbindlist(genetable_list, use.names = TRUE, fill = TRUE)
  TSgeneall <- rbindlist(TSgeneall_list, use.names = TRUE, fill = TRUE)

  TSgeneall_summarise <- TSgeneall[, .(
    supportcount = .N,
    tsgsumweight = sum(log10sample, na.rm = TRUE),
    meanspm = mean(spm, na.rm = TRUE)
  ), by = c("Site_ID", "tissue")]
  
  TSgeneall03 <- TSgeneall %>% filter(expre <= 0.3)

  genesumweight <- genetable[, .(
    gene_count = .N,
    genesumweight = sum(log10sample, na.rm = TRUE)
  ), by = c("Site_ID", "tissue")]
  
  tmp <- merge(TSgeneall_summarise, TSgeneall03, by = c("Site_ID", "tissue"), all.y = TRUE)
  tmp <- merge(tmp, genesumweight, by = c("Site_ID", "tissue"), all.x = TRUE)

  
  tmp$Tscore <- tmp$tsgsumweight / tmp$genesumweight
  Tscoretable <- tmp %>% as.data.table() %>% .[, -c("expre", "spm", "log10sample", "dataset","supportcount", "tsgsumweight",  "genesumweight")] %>% unique()
  Tscoretable$tissue <- gsub("-", " ", Tscoretable$tissue)
  Tscoretable$tissue <- gsub("\\*", ",", Tscoretable$tissue)
  Tscoretable <- Tscoretable %>% filter(gene_count >= 2) %>% filter(Tscore >= 0.6) 
  Tscoretable <- Tscoretable[,-c("gene_count")]
  names(Tscoretable)<-c("Site_ID","Tissue","Mean_spm_adjust","Tscore")
  Tscoretable$Mean_spm_adjust<-round(Tscoretable$Mean_spm_adjust,2)
  Tscoretable$Tscore<-round(Tscoretable$Tscore,2)
  fwrite(Tscoretable, file.path(output_dir, "TSUMs.csv"))
}

#' Identify Housekeeping hypermethylation sites
#' @param input_dir Directory containing input dataset expression files (.csv)
#' @param output_dir Directory to save intermediate and final results

identify_hkhm <- function(input_dir, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  count <- fread(file.path(input_dir, "count.txt"))[order(dataset)]
  count$log10sample <- log10(count$sample)
  count$sigmoidlog10_sample <- sapply(count$log10sample, sigmoid)
  
  filelist <- list.files(input_dir, pattern = "*.csv")
  pattern_list <- list()
  genetable_list <- list()
  i <- 1
  
  for (file_name in filelist) {
    tmp <- fread(file.path(input_dir, file_name))
    genetable1 <- tmp[, 1, with = FALSE]
    colnames(genetable1)[1] <- "Site_ID"
    genetable1$dataset <- file_name
    genetable1$log10sample <- count$sigmoidlog10_sample[i]
    
    pattern_h <- tmp %>% filter(if_all(2:ncol(tmp), ~ .x >= 0.7))
    if (nrow(pattern_h) > 0) {
      pattern_h <- pattern_h[, 1, with = FALSE]
      colnames(pattern_h)[1] <- "Site_ID"
      pattern_h$dataset <- file_name
      pattern_h$log10sample <- count$sigmoidlog10_sample[i]
      pattern_list[[i]] <- pattern_h
    }
    genetable_list[[i]] <- genetable1
    i <- i + 1
  }
  
  pattern <- rbindlist(pattern_list, use.names = TRUE, fill = TRUE)
  genetable <- rbindlist(genetable_list, use.names = TRUE, fill = TRUE)
  
  TSgeneall_summarise <- pattern[, .(supportcount = .N, tsgsumweight = sum(log10sample, na.rm = TRUE)), by = c("Site_ID")]
  genesumweight <- genetable[, .(gene_count = .N, genesumweight = sum(log10sample, na.rm = TRUE)), by = c("Site_ID")]
  
  tmp <- genesumweight[TSgeneall_summarise, on = c("Site_ID"), nomatch = 0]
  tmp$Tscore <- tmp$tsgsumweight / tmp$genesumweight 
  tmp <- tmp %>% filter(gene_count >= 2) %>% filter(Tscore >= 0.6)
  tmp <- tmp %>% as.data.table() %>% .[, -c("gene_count", "genesumweight", "supportcount", "tsgsumweight")] %>% unique()

  fwrite(tmp, file.path(output_dir, "HKHMs.csv"))
}

#' Identify Housekeeping hypomethylation sites
#' @param input_dir Directory containing input dataset expression files (.csv)
#' @param output_dir Directory to save intermediate and final results

identify_hkum <- function(input_dir, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  count <- fread(file.path(input_dir, "count.txt"))[order(dataset)]
  count$log10sample <- log10(count$sample)
  count$sigmoidlog10_sample <- sapply(count$log10sample, sigmoid)
  
  filelist <- list.files(input_dir, pattern = "*.csv")
  pattern_list <- list()
  genetable_list <- list()
  i <- 1
  
  for (file_name in filelist) {
    tmp <- fread(file.path(input_dir, file_name))
    genetable1 <- tmp[, 1, with = FALSE]
    colnames(genetable1)[1] <- "Site_ID"
    genetable1$dataset <- file_name
    genetable1$log10sample <- count$sigmoidlog10_sample[i]
    
    pattern_h <- tmp %>% filter(if_all(2:ncol(tmp), ~ .x <= 0.3))
    if (nrow(pattern_h) > 0) {
      pattern_h <- pattern_h[, 1, with = FALSE]
      colnames(pattern_h)[1] <- "Site_ID"
      pattern_h$dataset <- file_name
      pattern_h$log10sample <- count$sigmoidlog10_sample[i]
      pattern_list[[i]] <- pattern_h
    }
    genetable_list[[i]] <- genetable1
    i <- i + 1
  }
  
  pattern <- rbindlist(pattern_list, use.names = TRUE, fill = TRUE)
  genetable <- rbindlist(genetable_list, use.names = TRUE, fill = TRUE)
  
  TSgeneall_summarise <- pattern[, .(supportcount = .N, tsgsumweight = sum(log10sample, na.rm = TRUE)), by = c("Site_ID")]
  genesumweight <- genetable[, .(gene_count = .N, genesumweight = sum(log10sample, na.rm = TRUE)), by = c("Site_ID")]
  
  tmp <- genesumweight[TSgeneall_summarise, on = c("Site_ID"), nomatch = 0]
  tmp$Tscore <- tmp$tsgsumweight / tmp$genesumweight 
  tmp <- tmp %>% filter(gene_count >= 2) %>% filter(Tscore >= 0.6)
  tmp <- tmp %>% as.data.table() %>% .[, -c("gene_count", "genesumweight", "supportcount", "tsgsumweight")] %>% unique()
  fwrite(tmp, file.path(output_dir, "HKUMs.csv"))
}

# ==============================================================================
# SECTION 2: Expression Analysis Functions (Abundance Data - TSG/TSP)
# ==============================================================================

#' Identify Tissue-Specific Genes/Proteins (TSG/TSP)
#' @description Identifies tissue-specific expression from transcriptome or proteome data.
#' @param input_dir Directory containing input dataset expression files (.csv)
#' @param output_dir Directory to save intermediate and final results
#' @param tissue_number Total number of unique tissue types across all datasets
#' @param spm_cut Cutoff value for SPM_adjust value
#' @param spm_jar_path Path to the SPM.jar file


identify_tsg_tsp <- function(input_dir, output_dir, tissue_number, spm_cut, spm_jar_path = "SPM.jar") {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Load count file
  count_file <- file.path(input_dir, "count.txt")
  if (!file.exists(count_file)) stop("count.txt not found in input directory.")
  
  count <- fread(count_file)[order(dataset)]
  count$log10sample <- log10(count$sample)
  count$log10tissue <- log10(count$tissue)
  count$sigmoidlog10_tissue <- sapply(count$log10tissue, sigmoid)
  count$sigmoidlog10_tissue_diall <- count$sigmoidlog10_tissue / sigmoid(log10(tissue_number))
  count$sigmoidlog10_sample <- sapply(count$log10sample, sigmoid)
  
  filelist <- list.files(input_dir, pattern = "*.csv")
  TSgeneall_list <- list()
  genetable_list <- list()
  i <- 1
  
  for (file_name in filelist) {
    filepath <- file.path(input_dir, file_name)
    message(paste("Processing (Expression):", file_name))
    
    tmp_table <- data.table::fread(filepath)
    names(tmp_table)[1] <- "gene_symbol"
  
    tmp_table[tmp_table < 1] <- 0 
    
    names(tmp_table) <- gsub(" ", "-", names(tmp_table))
    
    # Prepare for SPM calculation
    temp_csv <- file.path(output_dir, "tmp.csv")
    write_csv(tmp_table, temp_csv)
    
    message(paste("Calculating SPM for", file_name))
    system(paste("java -jar", spm_jar_path, temp_csv))
    
    spm_output <- paste0(temp_csv, ".spm")
    new_spm_name <- "tmp_spm.csv"
    if(file.exists(spm_output)) file.rename(spm_output, new_spm_name)
    
    # Process SPM
    test_spm <- fread(new_spm_name)
    if("DPM" %in% names(test_spm)) test_spm[, DPM := NULL]
    
    test_spm <- pivot_longer(test_spm, !gene_symbol, names_to = "tissue", values_to = "spm")
    test_spm$spm <- test_spm$spm * count$sigmoidlog10_tissue_diall[i]
    
    genelist <- test_spm[, 1:2] 
    genelist$dataset <- file_name
    genelist$log10sample <- count$sigmoidlog10_sample[i]
    genetable_list[[i]] <- genelist
    
    # Filter TSG
    TSgene <- test_spm[test_spm$spm >= spm_cut, ]
    
    # Expression Quantile Logic
    tmp_table_long <- pivot_longer(tmp_table, !gene_symbol, names_to = "tissue", values_to = "expre")
    
    TSgene <- as.data.table(tmp_table_long)[as.data.table(TSgene), on = c("gene_symbol", "tissue"), nomatch = 0]
    TSgene$log10sample <- count$sigmoidlog10_sample[i]
    
    TSgeneall_list[[i]] <- TSgene
    i <- i + 1
  }
  
  genetable <- rbindlist(genetable_list, use.names = TRUE, fill = TRUE)
  TSgeneall <- rbindlist(TSgeneall_list, use.names = TRUE, fill = TRUE)
  
  # Summary Statistics
  TSgeneall_summarise <- TSgeneall[, .(
    supportcount = .N,
    tsgsumweight = sum(log10sample, na.rm = TRUE),
    meanspm = mean(spm, na.rm = TRUE)
  ), by = c("gene_symbol", "tissue")]
  
  # Summary Statistics 
  genesumweight <- genetable[, .(
    gene_count = .N,
    genesumweight = sum(log10sample, na.rm = TRUE)
  ), by = c("gene_symbol", "tissue")]
  
  # Final T-Score Calculation
  tmp <- merge(TSgeneall_summarise, TSgeneall, by = c("gene_symbol", "tissue"), all = TRUE)
  tmp <- merge(tmp, genesumweight, by = c("gene_symbol", "tissue"), all.x = TRUE)
  
  tmp$Tscore <- tmp$tsgsumweight / tmp$genesumweight
  
  Tscoretable <- tmp[, -c("expre", "spm", "log10sample","supportcount", "tsgsumweight", "genesumweight")] %>% unique()
  
  Tscoretable$tissue <- gsub("-", " ", Tscoretable$tissue)
  
  Tscoretable <- Tscoretable %>% filter(gene_count >= 2) %>% filter(Tscore >= 0.6)
  Tscoretable <- Tscoretable[,-c("gene_count")]
  
  names(Tscoretable)<-c("Gene_symbol","Tissue","Mean_spm_adjust","Tscore")
  
  Tscoretable$Mean_spm_adjust<-round(Tscoretable$Mean_spm_adjust,2)
  
  Tscoretable$Tscore<-round(Tscoretable$Tscore,2)
  
  fwrite(Tscoretable, file.path(output_dir, "TSG.csv"))
  message("TSG/TSP Identification Complete.")
  gc()
}