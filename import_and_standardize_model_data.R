library(tidyverse)

# Function to read and standardize CSV files
read_and_standardize_csv <- function(file_path) {
  df <- read_csv(file_path, col_types = cols(.default = "c")) # Read everything as character
  df
}

# Function to import CSVs and keep only selected columns
import_csv_files <- function(folder_path, columns_to_keep = NULL) {
  file_list <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)

  # Read all files into a list
  data_list <- lapply(file_list, read_and_standardize_csv)

  # Get all unique column names across all files
  all_columns <- unique(unlist(lapply(data_list, colnames)))

  # If no specific columns are provided, keep all columns
  if (is.null(columns_to_keep)) {
    columns_to_keep <- all_columns
  }

  # Standardize columns and filter to selected ones
  standardized_list <- lapply(data_list, function(df) {
    missing_cols <- setdiff(columns_to_keep, colnames(df)) # Identify missing columns
    df <- df %>% mutate(across(everything(), as.character)) # Ensure consistent types

    # Only add missing columns if they exist
    if (length(missing_cols) > 0) {
      df <- bind_cols(df, as_tibble(setNames(replicate(length(missing_cols), NA_character_, simplify = FALSE), missing_cols)))
    }

    # Keep only the desired columns
    df <- df %>% select(all_of(columns_to_keep))
    df
  })

  # Bind all data frames into one
  final_df <- bind_rows(standardized_list)

  return(final_df)
}

# Example usage: Keeping only selected columns
folder_path <- "Data/Originals"
columns_to_keep <- c("Return_Yr", "SPECIES", "Run", "TagDate", "TagReach",
                     "Tag1", "Tag2", "Sex", "FL", "Mark", "ScaleAge",
                     "RecapDate", "RecapTag1", "RecapTag2", "RecapReach")

combined_data <- import_csv_files(folder_path, columns_to_keep)

# View the structure of the final dataset
glimpse(combined_data)

# Check for column mismatches (version 1)
column_mismatches <- map_df(list.files(folder_path, pattern = "\\.csv$", full.names = TRUE), ~{
  df <- read_and_standardize_csv(.x)
  tibble(file = basename(.x), columns = paste(colnames(df), collapse = ", "))
})

print(column_mismatches, width=Inf)


# Check for column mismatches (version 2)
column_sets <- map(file_list, ~ colnames(read_and_standardize_csv(.x)))
column_differences <- map(column_sets, ~ setdiff(.x, column_sets[[1]])) # Identify the set of columns that are not in the first file

mismatch_summary <- tibble( # Combine results into a tibble for easier review
  file = basename(file_list),
  extra_columns = map_chr(column_differences, ~ paste(.x, collapse = ", "))
)

print(mismatch_summary, n = Inf, width = 10000)

