library(optparse)
library(tidyverse)

# define options
option_list <- list(
  make_option("--out", type = "character", help = "Output file")
)

parser <- OptionParser(usage = "Usage: %prog --out <output.tsv> [coverage_files...]", 
                       option_list = option_list)

opt <- parse_args(parser, positional_arguments = TRUE)

output_file <- opt$options$out
coverage_files <- opt$args

print(output_file)
print(coverage_files)
# checks
if (is.null(output_file)) {
  stop("Please provide --out")
}

if (length(coverage_files) == 0) {
  stop("Please provide at least one coverage file")
}

# main logic
cov_all <- map_dfr(coverage_files, function(f) {
  df <- read_tsv(f, skip = 1, show_col_types = FALSE)
  colnames(df) <- c("chr", "start", "end", "cnt")
  df %>%
    mutate(file = basename(f)) %>%
    select(chr, start, end, file, cnt)
})

write_tsv(cov_all, output_file)