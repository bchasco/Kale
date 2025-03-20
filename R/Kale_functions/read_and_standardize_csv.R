## Function to read and standardize CSV files
read_and_standardize_csv <- function(file_path) {
  df <- read_csv(file_path, col_types = cols(.default = "c")) # Read everything as character
  df
}
