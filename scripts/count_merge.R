library(tidyverse)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)
counts_dir <- args[1]
output_file <- args[2]

safe_read <- function(file) {
  tryCatch({
    # Extract library and length class from filename
    filename <- basename(file)
    library <- str_extract(filename, "^[^_]+")
    len_class <- if(str_detect(filename, "short")) "short" else "long"
    
    # Read and process file
    data <- fread(file, header = TRUE) %>%
      as_tibble() %>%
      # Add library and length class
      mutate(
        library = library,
        len_class = len_class
      )
    
    message(sprintf("Successfully read %s with %d rows", filename, nrow(data)))
    return(data)
    
  }, error = function(e) {
    message(sprintf("Error reading %s: %s", basename(file), e$message))
    return(NULL)
  })
}


files <- list.files(counts_dir, pattern = "*.tmp$", full.names = TRUE)
message(sprintf("Found %d files to process", length(files)))


all_counts <- files %>%
  map(safe_read) %>%
  compact() %>%  
  bind_rows() %>%

  filter(
    !is.na(gc),
    !is.nan(gc),
    !is.infinite(gc),
    count > 0
  )


message(sprintf("Final dataset has %d rows after filtering", nrow(all_counts)))


if (nrow(all_counts) == 0) {
  stop("No valid data remains after filtering")
}


write_tsv(all_counts, output_file)
message(sprintf("Successfully wrote output to %s", output_file)) 