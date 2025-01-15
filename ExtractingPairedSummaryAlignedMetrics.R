# Read the file
SUMMARY_PIC <- readLines(file_path)
print(SUMMARY_PIC)

# Step 1: Obtain headers at line 13
header_line <- SUMMARY_PIC[13] 

# Step 2: Get data--> line 14 forward
data_lines <- SUMMARY_PIC[14:84] 

cat("Header Line:\n", header_line, "\n\n")
cat("Data Lines:\n", head(data_lines, 5), "\n")  # Present 5 lines of data

# Step 3: Split into column names
header <- strsplit(header_line, "\t")[[1]]

# Step 4: Create individual entries
split_data <- strsplit(data_lines, "\t")

# Step 5: Convert to a data frame
data_frame <- do.call(rbind, lapply(split_data, function(x) {
  as.data.frame(t(x), stringsAsFactors = FALSE)
}))

# Step 6: Header --> column names of data frame
colnames(data_frame) <- header

# Step 7: Make columns numeric if needed
data_frame$READ_LENGTH <- as.numeric(data_frame$READ_LENGTH)
data_frame$PAIRED_TOTAL_LENGTH_COUNT <- as.numeric(data_frame$PAIRED_TOTAL_LENGTH_COUNT)
data_frame$PAIRED_ALIGNED_LENGTH_COUNT <- as.numeric(data_frame$PAIRED_ALIGNED_LENGTH_COUNT)

write.csv(data_frame, file = file.path(dir_path, "LEN_PIC_matrix.csv"), row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~

# Checking values in lines 7 to 10 in the content
for (i in 7:10) {
  line_values <- strsplit(SUMMARY_PIC[i], "\t")[[1]]  
  num_values <- length(line_values)  
  cat("Line", i, "has", num_values, "values.\n")
}

###################CREATING MATRIX FOR THE PAIR, FIRST AND SECOND CTAEGORIES IN THE ALIGNED SUMMARY METRICS FOR ALL SAMPLES#################
process_sample_data <- function(main_dir) {
  data_file <- list.files(path = main_pathway, pattern = "*.alignment_summary_metrics", full.names = TRUE)
  if (length(data_file) == 0) {
    stop(paste("No data file found in folder:", main_dir))
  }
  
  SUMMARY_PIC <- readLines(data_file)
  
  headers <- strsplit(SUMMARY_PIC[7], "\t")[[1]][-1]
  headers <- c(headers[35], headers[-35])  
  remove_indices <- c(34, 35)
  headers <- headers[-remove_indices]
  
  data_rows <- list()

  for (line in SUMMARY_PIC[10]) {
    split_line <- strsplit(line, "\t")[[1]]
    remove_indices <- c(34, 35)
    split_line <- split_line[-remove_indices]
    row_label <- split_line[1]
    
    # Convert remaining data values to numeric
    numeric_values <- as.numeric(split_line[-1])
    
    if (length(numeric_values) != 32) {
      stop(paste("Error: Mismatch", row_label))
    }
    
    data_rows[[length(data_rows) + 1]] <- c(row_label, numeric_values)
  }
  
  result_df <- as.data.frame(do.call(rbind, data_rows), stringsAsFactors = FALSE)
  
  colnames(result_df) <- headers
  
  result_df[-1] <- lapply(result_df[-1], as.numeric)
  
  return(result_df)
}


# Get a list of all sample folders
pathway <- list.dirs(path = main_path, full.names = TRUE, recursive = FALSE)
main_path=main_path[grep("/StartSampleName", main_path)].  #Run for negative and positive samples
#main_path=main_dir[grep("/Negative_control", main_dir)]

results_list <- list()

for (main_dir in main_dir) {
  sample_result <- process_sample_data(main_dir)
  
  sample_id <- basename(main_dir)
  
  rownames(sample_result) <- sample_id  
  results_list[[sample_id]] <- sample_result  
}

# Merge sample results in one data frame
qualPIC_matrix_result <- do.call(rbind, results_list)

qualPIC_matrix_result$id=rownames(qualPIC_matrix_result)
rownames(qualPIC_matrix_result)= NULL
colnames(qualPIC_matrix_result)[ncol(qualPIC_matrix_result)] = "id"
qualPIC_matrix_result = qualPIC_matrix_result[, c(ncol(qualPIC_matrix_result), 1:(ncol(qualPIC_matrix_result)-1))]

write.csv(qualPIC_matrix_result, file = file.path(path_here, "Negative_PAIR_qualPIC_matrix_results.csv"), row.names = FALSE)
#write.csv(qualPIC_matrix_result, file = file.path(path_here, "PAIR_qualPIC_matrix_results.csv"), row.names = FALSE)



