library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if task_id is provided
if (length(args) == 0) {
  stop("Please provide a task_id as a command line argument")
}

task_id <- args[1]
input_type <- args[2]

# task_id = "protein_del"; input_type = "protein"
# task_id = "genomic_del"; input_type = "genomic"
# task_id = "transcript_del"; input_type = "transcript"

print(paste0('Input type: ', input_type))

BASE_DIR <- "/srv/shiny-server/pon_del"
TEMP_DIR <- file.path(BASE_DIR, "temp")

print(paste0('Running task id: ', task_id))
print(paste0('Adding pchange ...'))

# # Set working directory
# setwd('/srv/data/website/structure_static/lu_psb_rshiny_server/pon_del/')

# Read and process data
df = read.csv(file.path(TEMP_DIR, paste0(task_id, '.csv')), stringsAsFactors = FALSE)

# deal with different input types
if (input_type == "protein") {
  names(df) <- c('NP_id', 'start', 'end')
  map = read.csv('/srv/shiny-server/db/mapping/mane_dict.csv', stringsAsFactors = FALSE)
  map <- map %>%
    dplyr::select(NP_id, ENST_id)
  df <- df %>% left_join(map, by='NP_id')
  df <- df %>% mutate(p_change = paste0(ENST_id, ":p.", start, "_", end, 'del'))
  out = df$p_change
}

if (input_type == "genomic") {
  names(df) <- c('chr', 'start', 'end')
  df <- df %>% mutate(g_change = paste0("chr", chr, ":g.", start, "_", end, "del"))
  out = df$g_change
}

if (input_type == "transcript") {
  names(df) <- c('ENST_id', 'start', 'end')
  df <- df %>% mutate(c_change = paste0(ENST_id, ":c.", start, "_", end, 'del'))
  out = df$c_change
}

writeLines(out, file.path(TEMP_DIR, paste0(task_id, '_transvarInput.txt')))
print(paste0('Done!'))