read_list_from_text <- function(file) {
  lines <- readLines(file)

  list_obj <- list()
  current_name <- NULL
  current_data <- NULL
  is_dataframe <- FALSE
  data_rows <- list()

  for (line in lines) {
    line <- trimws(line)

    # Detect new list item
    if (grepl("^\\$", line)) {
      if (!is.null(current_name)) {
        if (is_dataframe) {
          list_obj[[current_name]] <- as.data.frame(do.call(rbind, data_rows))
        } else {
          list_obj[[current_name]] <- current_data
        }
      }

      # Start new list entry
      current_name <- sub("^\\$", "", line)
      current_data <- NULL
      data_rows <- list()
      is_dataframe <- FALSE
      next
    }

    # Detect numeric or character values
    if (grepl("^\\[.*\\]$", line)) {
      values <- as.numeric(unlist(strsplit(gsub("\\[.*\\]", "", line), "\\s+")))
      current_data <- values
      next
    }

    # Detect data frame format (multiple rows with columns)
    if (grepl("^[0-9]+\\s+", line)) {
      row_values <- unlist(strsplit(line, "\\s+"))
      data_rows[[length(data_rows) + 1]] <- row_values
      is_dataframe <- TRUE
    }
  }

  # Save the last entry
  if (!is.null(current_name)) {
    if (is_dataframe) {
      list_obj[[current_name]] <- as.data.frame(do.call(rbind, data_rows))
    } else {
      list_obj[[current_name]] <- current_data
    }
  }

  return(list_obj)
}

# Load the manually formatted file
x <- read_list_from_text("data/test.txt")
str(x)

