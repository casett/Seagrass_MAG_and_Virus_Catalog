# Import libraries
library(tidyverse)
library(vroom)

# Combine individual files
cov_mags_ncbi <- "data/mag_coverage/"   # path to the data
cov_mags_ncbi_files <- dir(cov_mags_ncbi, pattern = "*.txt") # get file names

# Load in coverage data
cov_mags_ncbi_data <- data_frame(filename = cov_mags_ncbi_files) %>% # create a data frame
  # holding the file names
  mutate(file_contents = map(filename,          # read files into
                             ~ read_tsv(file.path(cov_mags_ncbi, .))) # a new data column
  )  
cov_mags_ncbi_data

# Unnest data
cov_mags_ncbi_data_un <- unnest(cov_mags_ncbi_data, cols = c(file_contents))

# Fix names
cov_mags_ncbi_data_un <- cov_mags_ncbi_data_un %>% 
  mutate_at("filename", str_replace, "-mean_coverage.txt", "")

# Calculate mean coverage per bin
cov_mags_ncbi_results <- cov_mags_ncbi_data_un %>%
  select(-bin) %>%
  group_by(filename) %>%
  pivot_longer(cols = starts_with('CE'), values_to = "coverage") %>%
  summarize(mean_cov = mean(coverage)) 

# Save results
write_csv(cov_mags_ncbi_results, "results/cov_mags_ncbi_results.csv")
