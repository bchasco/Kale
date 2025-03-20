## Function to import CSVs and keep only selected columns
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

  })

  # Bind all data frames into one
  final_df <- bind_rows(standardized_list)

  return(list(file_list = file_list, final_df=final_df))
}
