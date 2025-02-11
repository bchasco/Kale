# Parameters for the simulation
n_tag_weeks <- 4  # Number of weeks for tagging
n_rec_weeks <- 4  # Number of weeks for recovery
total_weeks <- n_tag_weeks + n_rec_weeks
gamma <- 0.4  # Probability of capturing a fish each week
omega <- 0.9  # Weekly survival probability

#arrival process
arrival_probs <- rep(0.25,4)
n_tags_week <- rmultinom(1, 400, prob = arrival_probs)


# Simulate data
set.seed(42)

#Create the fish data
fish_data <- list()

for (week in 1:n_weeks_tagging) {
  for (fish_id in 1:n_tags_week[week]) {
    fish_row <- rep(0, total_weeks)  # Initialize capture history
    fish_row[week] <- 1  # Fish is tagged in the current week
    alive <- TRUE

    for (t in (week + 1):total_weeks) {
      if (alive) {
        # Check if fish survives to the next week
        if (runif(1) > survival_prob) {
          alive <- FALSE
        }
        # If alive, determine if it is captured
        else if (runif(1) < capture_prob) {
          fish_row[t] <- 1
        }
      }
    }

    fish_data <- append(fish_data, list(c(paste0("Fish_", week, "_", fish_id), fish_row)))
  }
}

# Convert to data frame
columns <- c("Fish_ID", paste0("Week_", 1:total_weeks))
simulated_data <- do.call(rbind, fish_data)
colnames(simulated_data) <- columns

simulated_data <- as.data.frame(simulated_data)
simulated_data[, -1] <- lapply(simulated_data[, -1], as.numeric)  # Convert capture history to numeric

# Display the first few rows of the simulated data
head(simulated_data)
