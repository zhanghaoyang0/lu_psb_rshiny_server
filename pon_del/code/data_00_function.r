libs <- c("data.table", "dplyr", "tidyr", "stringi", 'stringr', 
    'GenomicRanges', 'tools')
lapply(libs, require, character.only = TRUE)
options(width = 500)

# Function to load and filter ponp3 data
load_ponp3 <- function(protein_seqs, file_ponp3 = '/srv/shiny-server/db/PON_P3/PONP3_prediction_mean.csv.gz') {
    ponp3 <- fread(file_ponp3)
    ponp3 <- ponp3 %>% filter(NP_id %in% names(protein_seqs))
    return(ponp3)
}

# Get basic info
get_basic_info <- function(df, map_mane, protein_seqs, 
                          AFmapping_file = "/srv/shiny-server/db/mapping/refseq_prot_to_alphafold.csv",
                          ppi_data = ppi,
                          gene_age_file = '/srv/shiny-server/db/all_data_files/gene_age_data.csv') {

  # Introduction               
  cols_to_add <- c("ENST_id", 'gene_symbol', "chr", 'protein_len', "uniprot_id_AF",  names(ppi)[1:6], "gene_age")
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts

  # Add protein length
  tmp <- sapply(protein_seqs, nchar)
  info <- data.frame(NP_id = names(protein_seqs), protein_len = tmp)
  df <- df %>% merge(info)

  # Add uniprot_id_AF mapping
  mapping <- read.csv(AFmapping_file) %>% 
    dplyr::select(refseq_prot, uniprot_accession)
  df <- df %>%
    merge(mapping, by.x = "NP_id", by.y = "refseq_prot", all.x = T) %>%
    dplyr::rename(uniprot_id_AF = uniprot_accession)
  
  # Calculate relative position and closest distance
  df$relative_pos <- df$start/df$protein_len
  df$closest_distance <- pmin(df$start-1, df$protein_len - df$end)
  
  # Get chromosome information
  map <- map_mane %>%
    dplyr::select(NP_id, NC_id) %>%
    mutate(chr = as.integer(str_extract(NC_id, "(?<=NC_0{0,3})[0-9]+")))
  map <- map %>%
    mutate(chr = case_when(
      chr == 23 ~ "chrX",
      chr == 24 ~ "chrY",
      TRUE ~ paste0("chr", chr)
    )) %>% dplyr::select(-NC_id)
  df <- df %>% left_join(map, by='NP_id')
  
  # Add PPI information
  df <- df %>% left_join(ppi_data, by='NP_id')
  
  # Add gene age information
  gene_age <- read.csv(gene_age_file) %>%
    dplyr::rename(NP_id = refseq_ids) %>%
    dplyr::select(NP_id, gene_age)
  df <- df %>% left_join(gene_age, by='NP_id')
  
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get genomic positions, protein positions, NP_id depend on input type
get_pos <- function(df, transvar_output, map_mane, input_type) {
  # process transvar result
  info = read.delim(transvar_output, header=1)
  info <- info %>%
    mutate(
      g_change = str_extract(coordinates.gDNA.cDNA.protein., "(?<=:)g\\.[^/]+"),
      g_start = str_extract(g_change, "\\d+"),
      g_end = as.numeric(sub("g\\.\\d+_(\\d+)del.*", "\\1", g_change)),
      p_change = str_extract(coordinates.gDNA.cDNA.protein., "p\\.[^/]+"),
      start = str_extract(p_change, "\\d+"),
      end = if_else(
        str_detect(p_change, "_\\D*\\d+"),
        str_extract(p_change, "_\\D*(\\d+)") %>% str_extract("\\d+"),
        start
      ),
      transcript = sub(" .*", "", transcript)
    ) %>%
    select(input, transcript, g_start, g_end, start, end) %>%
    distinct() %>% na.omit()
  # filter to mane
  map = map_mane %>% mutate(transcript=sub("\\..*", "", ENST_id))
  info = info %>% left_join(map, by= 'transcript')

  if (input_type == "protein") {
    # Map to ENST_id and create p_change
    map <- map_mane %>%
      dplyr::select(NP_id, ENST_id, gene_symbol)
    df <- df %>% left_join(map, by='NP_id')
    df <- df %>% mutate(input = paste0(ENST_id, ":p.", start, "_", end, 'del'))
    df = df%>%left_join(info%>%select(input, g_start, g_end), by='input')
  }

  if (input_type == "genomic") {
    df = df%>%mutate(input = paste0('chr', chr, ':g.', start, '_', end, 'del'))%>%
      dplyr::rename(g_start=start, g_end=end)
    df = df%>%left_join(info%>%select(input, start, end, NP_id), by='input')
  }

  if (input_type == "transcript") {
    df = df%>%mutate(input = paste0(ENST_id, ':c.', start, '_', end, 'del'))%>%select(-start, -end)
    df = df%>%left_join(info, by='input')
  }
  df = df%>%mutate(start=as.numeric(start), end=as.numeric(end), g_start=as.numeric(g_start), g_end=as.numeric(g_end))

  return(df)
}

# Get del_seq, up5seq and down5seq
get_del_seq <- function(df, protein_seqs) {
    # Introduction               
    cols_to_add <- c("del_seq", "up5seq", "down5seq")
    df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
    print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))

  # Initialize sequence columns
  df <- df %>% mutate(del_seq = NA, up5seq = NA, down5seq = NA)
  
  # Extract sequences
  for (i in seq_len(nrow(df))) {
    if (i %% 100 == 0) {
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    sub_seq <- protein_seqs[df[i, "NP_id"]]
    df[i, "del_seq"] <- substr(sub_seq, df[i, "start"], df[i, "end"])
    df[i, "up5seq"] <- substr(sub_seq, df[i, "start"] - 5, df[i, "start"] - 1)
    df[i, "down5seq"] <- substr(sub_seq, df[i, "end"] + 1, df[i, "end"] + 5)
  }

  # # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get aaindex
# Optimized get_aaindex function with vectorized operations, can address U
get_aaindex <- function(df, aaindex_file = '/srv/shiny-server/db/aaindex/aaindex1.csv') {
  # Read and prepare AAindex data
  aaindex <- col_to_num(read.csv(aaindex_file))
  rownames(aaindex) <- paste0('aaindex_', aaindex[,1])
  aaindex <- as.matrix(aaindex[, -1])
  
  cols_to_add <- rownames(aaindex)
  df <- df[, !names(df) %in% cols_to_add]
  
  # Function to calculate AAindex for a sequence (vectorized)
  calculate_aaindex_vectorized <- function(seqs) {
    # Split all sequences at once
    seq_chars <- strsplit(seqs, '')
    
    # Calculate AAindex for all sequences
    result <- t(sapply(seq_chars, function(seq) {
      if (length(seq) > 1) {
        # Check if all characters exist in aaindex
        valid_chars <- seq[seq %in% colnames(aaindex)]
        if (length(valid_chars) > 0) {
          return(rowMeans(aaindex[, valid_chars, drop = FALSE]))
        } else {
          return(rep(NA, nrow(aaindex)))
        }
      } else if (length(seq) == 1 && seq %in% colnames(aaindex)) {
        return(aaindex[, seq])
      } else {
        return(rep(NA, nrow(aaindex)))
      }
    }))
    
    return(result)
  }
  
  # Calculate AAindex features for all deletion sequences at once
  tryCatch({
    add <- calculate_aaindex_vectorized(df$del_seq)
    
    # Add features to data frame
    df <- cbind(df, add)
  }, error = function(e) {
    print(paste("Error in AAindex calculation:", e$message))
    print("Sample sequences causing issues:")
    print(head(df$del_seq))
    stop(e)
  })
  
  return(df)
}


# Get thermodynamic features (Gw_U, Gs_U, W_U)
get_thermo_feat <- function(df, thermo_dir = '/srv/shiny-server/db/ProtDCal/thermo_dynaminc') {
  # Introduction               
  cols_to_add <- c('Gw_U', 'Gs_U', 'W_U')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding thermodynamic feature columns
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))
  df[,cols_to_add] = NA
  
  # Calculate thermodynamic features
  cat("Calculating thermodynamic features...\n")
  for (i in 1:nrow(df)){
    if (i%%500==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    path_file = sprintf('%s/%s.txt', thermo_dir, df[i, 'NP_id'])
    if (!file.exists(path_file)){next}
    sub = read.table(path_file)
    sub[] = lapply(sub, function(x) if(is.character(x)) gsub(',', '', x) else x) # treat some number like "-1,015.709"
    sub = col_to_num(sub)
    df[i, cols_to_add] = colMeans(sub[df[i, 'start']:df[i, 'end'], 2:4])
  }
  
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get ISA and ECI features
get_isa_eci <- function(df, aa_data_file = '/srv/shiny-server/db/ProtDCal/amino_acid_data.csv') {
  # Introduction               
  cols_to_add <- c('ISA', 'ECI')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding ISA and ECI features
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH del_seq'))

  # Calculate ISA and ECI features
  cat("Calculating ISA and ECI features...\n")
  info <- read.csv(aa_data_file)
  info = info[,c('aa', 'ISA', 'ECI')]
  df[,c('ISA', 'ECI')] = NA
  for (i in seq_len(nrow(df))){
    if (i%%500==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    seq = df[i, 'del_seq']
    seq = strsplit(seq, '')[[1]]
    df[i, c('ISA', 'ECI')] = colMeans(info%>%filter(aa%in%seq)%>%dplyr::select(-aa))
  }
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get Conservation scores
get_conservation <- function(df, cons_dir = '/srv/shiny-server/db/conservation/last2col') {
  # Introduction               
  cols_to_add <- c('Conservation')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding conservation score
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))
  df$Conservation = NA
  for (i in 1:nrow(df)){
    if (i%%500==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    path_file = sprintf('%s/%s.csv', cons_dir, df[i, 'NP_id'])
    if (!file.exists(path_file)){next}
    sub = read.csv(path_file)
    scores = sub[df[i, 'start']:df[i, 'end'], 'score1']
    df[i, 'Conservation'] = mean(scores)
  }
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get Accessibility scores
get_accessibility <- function(df, acc_dir = '/srv/shiny-server/db/6_accessibility_area_freesasa') {
  # Introduction               
  cols_to_add <- c('Accessibility')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding accessibility score
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH uniprot_id_AF, start, end'))
  df$Accessibility = NA
  for (i in 1:nrow(df)){
    if (i%%500==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    path_file = sprintf('%s/%s.txt', acc_dir, df[i, 'uniprot_id_AF'])
    if (!file.exists(path_file)){next}
    sub = read.table(path_file, header=1)
    scores = sub[df[i, 'start']:df[i, 'end'], 'new_access_area']
    df[i, 'Accessibility'] = mean(scores)
  }
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get Secondary Structure features
get_secondary_structure <- function(df, ss_dir = '/srv/shiny-server/db/secondary_structure/process') {
  # Introduction               
  cols_to_add = c('Alphahelix', 'Beta-sheet', 'B_Beta','b_Beta','T/C', 'G_helix', 'pi_helix', 'vlow_conf')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding secondary structure features
  # print(paste('Getting secondary structure features', paste(cols_to_add, collapse = ", "), 'WITH uniprot_id_AF, start, end'))
  df[,cols_to_add] = NA
  for (i in 1:nrow(df)){
    if (i%%500==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    path_file = sprintf('%s/%s.csv', ss_dir, df[i, 'uniprot_id_AF'])
    if (!file.exists(path_file)){next}
    sub = read.csv(path_file)
    tab = table(sub[df[i, 'start']:df[i, 'end'], 'structure'])
    for (j in cols_to_add){
      df[i, j] = ifelse(is.na(tab[j]), 0, 1) # binary
    }
  }
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get ponp3 score
get_ponp3_score <- function(df, ponp3) {
  # Introduction               
  cols_to_add <- c("ponp3_score")
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding ponp3 score
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))
  df$ponp3_score = NA
  # Read and filter PON-P3 data
  ponp3_filtered = ponp3%>%filter(NP_id%in%df$NP_id)
  # Calculate mean probability score for each deletion
  for (i in seq_len(nrow(df))) {
    if (i%%50==0){
      cat(sprintf("Processed %.1f%% (%d/%d rows)\n", 
                  i/nrow(df)*100, i, nrow(df)))
    }
    NP_id <- df$NP_id[i]  # Extract NP_id
    start <- df$start[i]  # Extract start position
    end <- df$end[i]  # Extract end position
    # Filter ponp3 for relevant positions
    sub <- ponp3_filtered %>% 
        filter(NP_id == NP_id & position >= start & position <= end)
    # Compute mean probability score
    score <- mean(sub$meanProb)  # Avoid issues with NA values
    # Store result back in df
    df$ponp3_score[i] <- score
  }
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get g_exon
get_granges_gFeats <- function(df, granges = granges[['g_exon']]) {
  # Introduction               
  cols_to_add <- c("g_exon")
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding g_exon
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))

  # pseudo id
  df$tmp_id = 1:nrow(df)
  # Create GRanges object for non-NA positions
  df$g_exon = 0
  valid_pos <- !is.na(df$g_start) & !is.na(df$g_end)
  if (sum(valid_pos) > 0) {
    gr = with(df[valid_pos,], GRanges(seqnames = chr, ranges = IRanges(start = g_start, end = g_end)))
    hits = findOverlaps(query = gr, subject = granges)
    df$tmp_id[valid_pos][queryHits(hits)] -> valid_hits
    df[df$tmp_id %in% valid_hits, 'g_exon'] = 1
  }
  # remove pseudo id
  df$tmp_id = NULL
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get protein granges features
get_granges_pFeats <- function(df, feats=get('p_', names(granges))) {
  # Introduction               
  cols_to_add <- feats
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding protein granges features
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))

  # Get protein-related features
  gr = with(df, GRanges(seqnames = NP_id, ranges = IRanges(start = start, end = end)))
  
  for (i in feats) {
    df[,i] = 0
    hits = findOverlaps(query = gr, subject = granges[[i]])
    df[queryHits(hits), i] = 1
    # Print distribution information for protein granges features
    # print(paste0(i, ' distribution:'))
    # print(table(df[,i]))
  }
  
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get g_exon
get_granges_gFeats <- function(df, granges) {
  # Introduction               
  cols_to_add <- c("g_exon")
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding g_exon
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH NP_id, start, end'))

  # pseudo id
  df$tmp_id = 1:nrow(df)
  # Create GRanges object for non-NA positions
  df$g_exon = 0
  valid_pos <- !is.na(df$g_start) & !is.na(df$g_end)
  if (sum(valid_pos) > 0) {
    gr = with(df[valid_pos,], GRanges(seqnames = chr, ranges = IRanges(start = g_start, end = g_end)))
    hits = findOverlaps(query = gr, subject = granges)
    df$tmp_id[valid_pos][queryHits(hits)] -> valid_hits
    df[df$tmp_id %in% valid_hits, 'g_exon'] = 1
  }
  # remove pseudo id
  df$tmp_id = NULL
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}

# Get gene group features
get_gene_groups <- function(df, gene_groups) {
  # Introduction
  keep = c('Housekeeping', 'Essential', 'Haploinsufficient', 
    "Redundancy", "Pseudogene", "Duplicate") # oringal is c("Redundancy", "Pseudogene", "Duplicate")
  cols_to_add <- keep
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding gene group features
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH gene_symbol'))
  
  for (i in keep) {
    genes = gene_groups[[i]]
    df[, i] = ifelse(df$gene_symbol %in% genes, 1, 0)
    # Print distribution information for gene group features
    # print(paste0(i, ' distribution:'))
    # print(table(df[,i]))
  }
  
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}






# Get bitScores
calculate_seq_bitScore <- function(seq, bitScores) {
  aas <- unlist(str_split(seq, ""))
  scores <- bitScores$mean_bit[match(aas, bitScores$letter)]
  mean(scores, na.rm = TRUE)
}

get_bitScores <- function(df, seq_cols=c('up5seq', 'down5seq'), bitScores_path='/srv/shiny-server/db/pon_del/data/bitScores.csv') {
  # Introduction               
  cols_to_add <- paste0(seq_cols, '_bitScore')
  df <- df[, !names(df) %in% cols_to_add] # Remove columns to avoid conflicts
  # Print status message about adding bit scores
  # print(paste('Getting', paste(cols_to_add, collapse = ", "), 'WITH', seq_cols))
  
  # Load bitScores data
  bitScores_df <- read.csv(bitScores_path)
  
  for (col in seq_cols) {
    # Get bitScores for this sequence type
    bitScores <- bitScores_df %>% 
      filter(type == col) %>%
      dplyr::select(-type)
    
    # Apply to all sequences
    seqs <- df %>% pull(col)
    seq_bitScores <- sapply(seqs, calculate_seq_bitScore, bitScores = bitScores)
    
    # Add bit score column
    col_name <- paste0(col, '_bitScore')
    df[,col_name] <- seq_bitScores
  }
  
  # Report summary
  # print_column_info(df, cols_to_add)
  return(df)
}



# ======================================================================
# Function in analysis paper
# ======================================================================

aamap = c('CYS', 'C', 'ASP', 'D', 'SER', 'S', 'GLN', 'Q', 'LYS', 'K','ILE', 'I', 'PRO', 'P', 'THR', 'T', 'PHE', 
'F', 'ASN', 'N', 'GLY', 'G', 'HIS', 'H', 'LEU', 'L', 'ARG', 'R', 'TRP', 'W', 'ALA', 'A', 'VAL','V', 'GLU','E', 'TYR', 'Y', 'MET', 'M')
aamap = data.frame(matrix(aamap, ncol=2, byrow=T))
names(aamap) = c('key', 'value')


structures = c('Alpha-helix', 'Beta-sheet', 'B-Beta','b-Beta','T/C', 'G-helix', 'pi-helix', 'vlow conf')
# ============================================================
# info
# ============================================================
# hydrophobic (V, I, L,F, M, W, Y, C), negatively charged (D, E), 
# positively charged (R, K, H), conformational (G, P), polar (N, Q, S), and others (A, T) [15].
amino_groups <- list(
    h = c("V", "I", "L", "F", "M", "W", "Y", "C"), n = c("D", "E"), p = c("R", "K", "H"),
    c = c("G", "P"), l = c("N", "Q", "S"), o = c("A", "T")
)

nt_groups <- list(
    "keto" = c("A", "C"), "amino" = c("T", "G"),
    "strong" = c("G", "C"), "weak" = c("T", "A"),
    "purine" = c("A", "G"), "pyrimidine" = c("T", "C")
)


aamap = c('CYS', 'C', 'ASP', 'D', 'SER', 'S', 'GLN', 'Q', 'LYS', 'K','ILE', 'I', 'PRO', 'P', 'THR', 'T', 'PHE', 
'F', 'ASN', 'N', 'GLY', 'G', 'HIS', 'H', 'LEU', 'L', 'ARG', 'R', 'TRP', 'W', 'ALA', 'A', 'VAL','V', 'GLU','E', 'TYR', 'Y', 'MET', 'M')
aamap = data.frame(matrix(aamap, ncol=2, byrow=T))
names(aamap) = c('key', 'value')


# ======================================================================
# functions to collect features
# ======================================================================
# Dipeptide composition, e.g,m cal_dipeptide('CDS')
cal_dipeptide <- function(sequence) {
    # Generate all possible dipeptides
    dipeptides <- outer(aamap[,2], aamap[,2], paste0)
    # Initialize a named vector to store counts of each dipeptide
    dp_counts <- setNames(rep(0, length(dipeptides)), dipeptides)
    # Count occurrences of each dipeptide in the sequence
    for (i in 1:(nchar(sequence) - 1)) {
        dp <- substr(sequence, i, i + 1)
        if (dp %in% dipeptides) {
            dp_counts[dp] <- dp_counts[dp] + 1
        }
    }
    # Calculate total count of dipeptides
    total_count <- sum(dp_counts)
    # Calculate dipeptide composition
    dp_comp <- dp_counts / total_count
    # Round the composition values to 4 decimal places
    rounded_dp_comp <- round(dp_comp, 4)
    # Convert to a data frame
    dpc_df <- as.data.frame(t(rounded_dp_comp))
    colnames(dpc_df) <- dipeptides
    return(dpc_df)
}



# ======================================================================
# functions
# ======================================================================
## print_column_info
print_column_info <- function(df, cols) {
  # Check if the specified columns exist in the dataframe
  missing_cols <- setdiff(cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste("The following columns are not in the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  # Loop through specified columns and print information
  for (col in cols) {
    na_count <- sum(is.na(df[[col]]))
    total_count <- length(df[[col]])
    na_prop <- round((na_count / total_count) * 100, 2)
    cat("Column:", col, "- NA proportion:", na_prop, "%\n")
    
    # Check if the column is binary (contains only two unique non-NA values)
    unique_vals <- unique(na.omit(df[[col]]))
    if (length(unique_vals) == 2) {
      prop_table <- round(prop.table(table(df[[col]], useNA = "no")) * 100, 2)
      cat("  Class proportions:", paste(names(prop_table), prop_table, "%", collapse = ", "), "\n")
    } else {
      # Print summary statistics for the column
      # print(summary(df[[col]]))
    }
  }
}


## extract element with pattern from vector
# e.g., get('a', c('AS', 'Ba', 'C')) ; get('a', c('AS', 'Ba', 'C'), exact=F)
get <- function(key, vec, exact = F, v = F) {
    vec <- unlist(vec)
    if (exact == T) {
        out <- vec[grep(key, vec, invert = v)]
    } else {
        out <- vec[grep(toupper(key), toupper(vec), invert = v)]
    }
    return(out)
}

is_numeric <- function(x) {
  !any(is.na(suppressWarnings(as.numeric(na.omit(x))))) & is.character(x)
}

# find deletion from fasta
find_del_from_fasta = function(df, seq=seq){
    df = df%>%mutate(ref=NA, alt='-')
    for (i in 1:nrow(df)){
        chr = as.character(df[i, 'chr'])
        start = df[i,] %>% pull(start); end = df[i, ]%>% pull(end)
        sub_seq = seq[chr]
        if (width(sub_seq) <= end | start>end) {
                df[i, 'ref'] = NA 
            } else {
                extracted <- as.character(subseq(sub_seq, start, end))
                df[i, 'ref'] = extracted
            }
        df[i, 'pos'] = df[i, 'start']
    }
    return(df)
}


## find range from
# e.g., find_range(c(1, 3:5, 7))
find_range = function(postions){
    diff_pos <- diff(postions)
    start_idx <- which(diff_pos != 1) + 1
    if (start_idx[1]!=1|length(start_idx)==0) {start_idx=c(1, start_idx)}
    end_idx <- c(start_idx[-1] - 1, length(postions))
    out = c()
    for (i in seq_along(start_idx)) {
        out = c(out, postions[start_idx[i]], postions[end_idx[i]])
    }
    res = data.frame(matrix(out, ncol=2, byrow=T))
    names(res) = c('start', 'end')
    return(res)
}


# conver col to numeric if possible
col_to_num <- function(df) {
    for (i in 1:ncol(df)){
        if (any(is.na(suppressWarnings(as.numeric(df[, i]))))) {next}
        df[,i] = as.numeric(df[, i])
    }
    return(df)
}

# convert p to star
p2star = function(p){
    p <- case_when(p < 0.001 ~ "***",
                p < 0.01 ~ "**",
                p < 0.05 ~ "*",
                TRUE ~ '')
    return(p)
}

# adjust p value in a dataframe
p.adjust_df = function(df){
    p =  p.adjust(unlist(df), method='BH') # extend by col
    df1 = matrix(p, nrow = nrow(df), ncol = ncol(df))
    df1 = as.data.frame(df1)
    rownames(df1) = rownames(df); colnames(df1) = colnames(df)
    return(df1)
}

# convert distribution mat to long formate
get_mat_long <- function(mat, calProp = F, valueName='n') {
    mat_long <- mat %>% mutate(from = as.factor(rownames(mat)))
    mat_long <- melt(setDT(mat_long), id.vars = c("from"), variable.name = "to", value.name = valueName)
    if (calProp == T) {
        mat_long <- mat_long %>% mutate(prop = n / sum(n))
    }
    return(mat_long)
}

get_dismat <- function(col, groups = NULL) {
    from <- substr(col, 1, 1)
    to <- substr(col, nchar(col), nchar(col))
    # if we need to group
    if (!is.null(groups)) {
        for (group in names(groups)) {
            from[from %in% groups[[group]]] <- group
            to[to %in% groups[[group]]] <- group
        }
    }
    # distribute vec and tab
    types <- sort(unique(c(from, to)))
    vec <- c()
    for (i in types) {
        for (j in types) {
            n <- sum(from == i & to == j)
            vec <- c(vec, n)
        }
    }
    mat <- as.data.frame(matrix(vec, ncol = length(types), byrow = T))
    colnames(mat) <- rownames(mat) <- types
    # wide to lone
    mat_long <- get_mat_long(mat, calProp = T, valueName = "n")
    # check
    if (sum(mat) != length(from) | sum(mat) != length(to)) {
        stop("check it!")
    }
    # res
    res <- list()
    res[["mat"]] <- mat
    res[["mat_long"]] <- mat_long
    return(res)
}


# compare two matrix with fisher test
fisher_onMat <- function(mat1, mat2, method) {
    if (!method %in% c("row", "col", "full")) {
        stop("method error")
    }
    types <- colnames(mat1)
    res <- list()
    if (method == "full") {
        mat_or <- mat_p <- mat1 %>% mutate_all(~NA)
        for (i in types) {
            for (j in types) {
                out <- c(mat1[i, j], mat2[i, j], sum(mat1) - mat1[i, j], sum(mat2) - mat2[i, j])
                crosstab <- matrix(out, nrow = 2, byrow = T, dimnames = list(c("yes", "no"), c("ref", "alt")))
                test <- fisher.test(crosstab) # fisher
                mat_or[i, j] <- test$estimate
                mat_p[i, j] <- test$p.value
            }
        }
    } else {
        mat_or <- mat_p <- mat1[,1, drop = FALSE] %>% mutate_all(~NA)
        colnames(mat_or)=colnames(mat_p)=sprintf('BY%s', method)
        mat_prop = data.frame()
        for (i in types) {
            if (method == "row") {
                out <- c(sum(mat1[i, ]), sum(mat2[i, ]))
            } else {
                out <- c(sum(mat1[, i]), sum(mat2[, i]))
            }
            out = c(out, sum(mat1)-out[1], sum(mat2)-out[2])
            crosstab <- matrix(out, nrow = 2, byrow = T, dimnames = list(c("yes", "no"), c("ref", "alt")))
            test <- fisher.test(crosstab) # fisher
            mat_or[i,] <- test$estimate
            mat_p[i,] <- test$p.value
            mat_prop = rbind(mat_prop, crosstab[1,]/colSums(crosstab)) # prop of yes
            names(mat_prop) = c('ref', 'alt')
            res[['prop']] = mat_prop # use for plot
        }
    }
    res[['or']] = mat_or
    res[['p']] = p.adjust_df(mat_p) # adjust multiple testing
    return(res)
}

# for chisqPosthoc_residualTab
chisqPosthoc_residualTab = function(input){
    tab = chisq.posthoc.test(input, method='BH')
    tab[tab=='p values'] = 'p-value'
    # residue tab
    tab1 = tab[seq(1, nrow(tab)-1, 2), 3:ncol(tab)] 
    tab1 = tab1 %>% mutate(across(everything(), ~ round(as.numeric(.), 2)))
    # p tab
    tab2 = tab[seq(2, nrow(tab), 2), 3:ncol(tab)]
    tab2 = tab2 %>% mutate(across(everything(), ~ as.numeric(gsub("\\*", "", .))))
    tab2 <- tab2 %>%mutate(across(everything(), ~ p2star(as.numeric(.))))
    # merge
    tmp <- tab1 %>% mutate(across(everything(), ~ paste0(., "", tab2[[cur_column()]])))
    out = cbind(tab[seq(1, nrow(tab)-1, 2), 1], tmp)
    names(out)[1] = 'var'
    return(out)
}
