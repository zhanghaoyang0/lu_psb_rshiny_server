# Set base directory
BASE_DIR <- "/srv/shiny-server/pon_del"
TEMP_DIR <- file.path(BASE_DIR, "temp")

# Source functions
source(file.path(BASE_DIR, 'code/data_00_function.r'))

# Load feature data
load('/srv/shiny-server/db/feat.rdata')

# Get task ID from command line argument
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide a task ID as a command line argument")
}
task_id <- args[1]
input_type <- args[2]

# task_id = "protein_del"; input_type = "protein"
# task_id = "genomic_del"; input_type = "genomic"
# task_id = "transcript_del"; input_type = "transcript"

print(paste0('Processing task ID: ', task_id))

# Load the corresponding CSV file
input_file <- file.path(TEMP_DIR, paste0(task_id, '.csv'))
if (!file.exists(input_file)) {
  stop(paste("File", input_file, "does not exist"))
}
print(paste0('Reading input file: ', input_file))
df <- read.csv(input_file)


# Add features
print('Adding features ...')

# Get genomic positions, protein positions, NP_id depend on input type
transvar_output <- file.path(TEMP_DIR, paste0(task_id, '_transvarOutput.txt'))
df <- get_pos(df, transvar_output, map_mane, input_type)

# Get basic info
df <- get_basic_info(df, map_mane, protein_seqs)
print('Added basic info')

# Get del_seq, up5seq and down5seqs
df <- get_del_seq(df, protein_seqs)
df <- df[df$del_seq != 'U',] # remove U
print('Added sequence features')

# Get aaindex
df <- get_aaindex(df)

# Get secondary structure features  
df <- get_secondary_structure(df)

# Get conservation
df <- get_conservation(df)

# Get accessibility
df <- get_accessibility(df)

# Get granges protein feats
df <- get_granges_pFeats(df)

# Get gene group features
df <- get_gene_groups(df, gene_groups)

# # Get bitScores
# df <- get_bitScores(df)


# remake protein input 
if (input_type == "protein") {
  df$input <- paste0(df$NP_id, ":p.", df$start, "_", df$end, 'del')
} 


# filter lgbm features
feats <- c('Haploinsufficient', 'Conservation', 'vlow_conf', 'Accessibility', 'p_repeat', 'protein_len', 'aaindex_RACS820103', 'Closeness', 'aaindex_NAKH920103', 'aaindex_NAKH900104', 'p_palindrome', 'aaindex_BIGC670101', 'T/C', 'aaindex_SUEM840102', 'aaindex_LEVM760103', 'aaindex_CHOP780204', 'aaindex_ARGP820102', 'aaindex_GEOR030105', 'Hub_Score', 'p_domain')
feats = c(feats, 'NP_id', 'start', 'end', 'input') # add start information, but not used in prediction

df = df[feats]

# Save the processed data
output_file <- file.path(TEMP_DIR, paste0(task_id, '_withfeat.csv'))
write.csv(df, output_file, row.names = FALSE)
print(paste0('Saved features to: ', output_file))
print(paste0('Output file has ', nrow(df), ' rows'))

print('Done!')