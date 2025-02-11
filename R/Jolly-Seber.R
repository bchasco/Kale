set.seed(123)

T <- 10  # Number of sampling occasions
n_fish <- 100  # Initial number of marked individuals

phi <- 0.8  # True survival probability
p <- runif(T, 0.2, 0.6)  # Capture probabilities per occasion

# Generate first capture histories
capture_histories <- matrix(0, n_fish, T)

for (i in 1:n_fish) {
  alive <- TRUE
  for (t in 1:T) {
    if (alive) {
      if (runif(1) < p[t]) {
        capture_histories[i, t] <- 1  # Captured
        break  # Tag is removed, no further recaptures
      }
      if (runif(1) > phi) {
        alive <- FALSE  # Dies or exits study
      }
    }
  }
}

# Convert to RTMB format
capture_data <- as.data.frame(capture_histories)

data <- list(ch = capture_data
             ,T = ncol(capture_data))

# Initial parameter values
parameters <- list(
  logit_phi = qlogis(0.8),
  logit_p = qlogis(runif(T, 0.2, 0.6)),
  log_N0 = log(100)
)

f <- function(parms){

  RTMB::getAll(data,
               parms)

  phi <- plogis(logit_phi)
  p <- plogis(logit_p)
  N0 <- exp(log_N0)

}
#   # Write model in C++
# model_code <- "
# #include <TMB.hpp>
#
# template<class Type>
# Type objective_function<Type>::operator() () {
#     DATA_MATRIX(capture_histories);  // Capture history matrix (individuals x time)
#     int T = capture_histories.cols(); // Number of sampling occasions
#
#     PARAMETER(logit_phi);  // Survival probability (logit scale)
#     PARAMETER_VECTOR(logit_p);  // Capture probabilities at each time step
#     PARAMETER(log_N0);  // Initial population size (log scale)
#
#     // Convert parameters to probability space
#     Type phi = invlogit(logit_phi);
#     vector<Type> p = invlogit(logit_p);
#     Type N0 = exp(log_N0); // Population size
#
#     // Likelihood
#     Type nll = 0;
#
#     // Number of individuals
#     int n_indiv = capture_histories.rows();
#
#     for (int i = 0; i < n_indiv; i++) {
#         bool captured = false;
#         Type prob = 1.0;  // Initialize probability
#
#         for (int t = 0; t < T; t++) {
#             if (capture_histories(i, t) == 1) {
#                 if (!captured) {
#                     // First capture event
#                     prob *= (1 - phi) * p(t);
#                     captured = true;
#                 }
#             } else if (!captured) {
#                 // Survived but not captured
#                 prob *= phi * (1 - p(t));
#             }
#         }
#
#         nll -= log(prob);  // Accumulate negative log-likelihood
#     }
#
#     REPORT(phi);
#     REPORT(p);
#     REPORT(N0);
#
#     return nll;
# }
# "
#
# # Compile the model
# writeLines(model_code, "jolly_seber_removal.cpp")
# compile("jolly_seber_removal.cpp")
# dyn.load(dynlib("jolly_seber_removal"))
#
# # Prepare Data for RTMB
# data_list <- list(
#   capture_histories = as.matrix(capture_data)
# )
#
# # Initial parameter values
# parameters <- list(
#   logit_phi = qlogis(0.8),
#   logit_p = qlogis(runif(T, 0.2, 0.6)),
#   log_N0 = log(100)
# )

obj <- RTMB::MakeADFun(f,
                       parameters,
                       silent = FALSE)
opt <- nlminb(obj$par,
              obj$fn,
              obj$gr)

# # Extract estimates
# rep <- obj$report()
#
# # Print Results
# list(
#   survival = plogis(opt$par["logit_phi"]),
#   capture_probs = plogis(opt$par["logit_p"]),
#   initial_population = exp(opt$par["log_N0"])
# )
