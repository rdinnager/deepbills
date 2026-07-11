# Utility function to validate and process input
process_input <- function(input, is_distance = NULL, verbose = TRUE) {
  if (is.null(is_distance)) {
    is_distance <- inherits(input, "dist") || 
                  (isSymmetric(as.matrix(input)) && all(diag(as.matrix(input)) == 0))
  }
  
  # Convert input to distance matrix
  if (!is_distance) {
    if(verbose) cli::cli_alert_info("Converting input data to distance matrix...")
    D <- as.matrix(Rfast::Dist(input))
  } else {
    if (inherits(input, "dist")) {
      D <- as.matrix(input)
    } else {
      D <- as.matrix(input)
    }
    
    if (!isSymmetric(D)) {
      cli::cli_alert_danger("Error: Distance matrix must be symmetric")
      stop("Distance matrix must be symmetric")
    }
    if (!all(diag(D) == 0)) {
      cli::cli_alert_danger("Error: Distance matrix diagonal must be zero")
      stop("Distance matrix diagonal must be zero")
    }
  }
  return(D)
}

# Function to determine fixed sigma for likelihood calculations
determine_sigma <- function(D, sigma = NULL, sigma_scale = 0.1, verbose = TRUE) {
  if(is.null(sigma)) {
    if(verbose) cli::cli_alert_info("Calculating fixed sigma for likelihood computations...")
    
    # Based on nearest neighbor distances
    nn_dists <- apply(D, 1, function(x) sort(x)[2])
    sigma <- sigma_scale * median(nn_dists)
    
    if(verbose) cli::cli_alert_info(sprintf("Using fixed sigma = %.4f for likelihood calculations", sigma))
  }
  return(sigma)
}

# Helper function for likelihood computation with fixed sigma
compute_likelihood <- function(observed_dist, predicted_dist, fixed_sigma) {
  upper_tri <- upper.tri(observed_dist)
  obs_vec <- observed_dist[upper_tri]
  pred_vec <- predicted_dist[upper_tri]
  sum(dnorm(obs_vec, mean = pred_vec, sd = fixed_sigma, log = TRUE))
}

# Main function
model_based_mds_selection <- function(input, max_dim = 10, 
                                    methods = c("AIC", "BIC", "bootstrap", "crossval", "stress"),
                                    n_bootstrap = 100, 
                                    n_folds = 5,
                                    seed = 123, 
                                    is_distance = NULL,
                                    verbose = TRUE,
                                    sigma = NULL,
                                    sigma_scale = 0.1) {
  require(MASS)
  require(cli)
  require(progress)
  require(Rfast)
  
  set.seed(seed)
  
  # Validate methods
  valid_methods <- c("AIC", "BIC", "bootstrap", "crossval", "stress")
  methods <- match.arg(methods, valid_methods, several.ok = TRUE)
  
  if(verbose) cli::cli_h1("Model-Based MDS Dimension Selection")
  
  # Process input
  D <- process_input(input, is_distance, verbose)
  n <- nrow(D)
  
  # Determine sigma for likelihood calculations
  sigma <- determine_sigma(D, sigma, sigma_scale, verbose)
  
  # Check dimensionality constraints
  max_possible_dim <- n - 1
  if (max_dim > max_possible_dim) {
    cli::cli_alert_warning(sprintf("Reducing max_dim from %d to %d due to data constraints",
                                 max_dim, max_possible_dim))
    max_dim <- max_possible_dim
  }
  
  # Initialize results structure
  results <- initialize_results(methods, max_dim, n, n_bootstrap, n_folds)
  
  # Perform MDS analysis
  results <- compute_mds_dimensions(D, results, methods, max_dim, n_bootstrap, n_folds, 
                                  sigma, verbose)
  
  # Post-process results
  results <- post_process_results(results, methods, max_dim)
  
  # Create visualizations
  if(verbose) create_visualizations(results, methods, max_dim)
  
  # Print summary
  if(verbose) print_summary(results, methods, sigma)
  
  return(results)
}

# Example usage:
# library(MASS)
# 
# # Generate example data
# data_matrix <- matrix(rnorm(100 * 20), nrow = 100)
# dist_matrix <- dist(data_matrix)
# 
# # Run with all methods
# results <- model_based_mds_selection(
#   dist_matrix,
#   max_dim = 10,
#   methods = c("AIC", "BIC", "bootstrap", "crossval"),
#   n_bootstrap = 100,
#   n_folds = 5,
#   verbose = TRUE
# )
# 
# # Run with selected methods only
# results_minimal <- model_based_mds_selection(
#   dist_matrix,
#   methods = c("AIC", "crossval"),
#   verbose = TRUE
# )

# Initialize results structure
initialize_results <- function(methods, max_dim, n, n_bootstrap, n_folds) {
  results <- list()
  
  if(any(c("AIC", "BIC") %in% methods)) {
    results$criteria <- matrix(NA, nrow = max_dim, ncol = 2, 
                             dimnames = list(NULL, c("AIC", "BIC")))
    results$likelihood_vals <- numeric(max_dim)
  }
  
  if("bootstrap" %in% methods) {
    results$bootstrap_stress <- array(NA, dim = c(max_dim, n_bootstrap))
  }
  
  if("crossval" %in% methods) {
    results$cv_error <- matrix(NA, nrow = max_dim, ncol = n_folds)
  }
  
  if("stress" %in% methods) {
    results$stress_vals <- numeric(max_dim)
  }
  
  results$eigenvalues <- matrix(NA, nrow = max_dim, ncol = n)
  
  return(results)
}

# Compute MDS for each dimension
compute_mds_dimensions <- function(D, results, methods, max_dim, n_bootstrap, n_folds, 
                                 sigma, verbose) {
  n <- nrow(D)
  
  # Create progress bar for dimensions
  if(verbose) {
    cli::cli_h2("Computing MDS solutions")
    pb_dim <- progress::progress_bar$new(
      format = "Dimension :dim [:bar] :percent eta: :eta",
      total = max_dim,
      clear = FALSE
    )
  }
  
  # Main computation loop
  for(k in 1:max_dim) {
    if(verbose) {
      pb_dim$tick(tokens = list(dim = sprintf("%d/%d", k, max_dim)))
    }
    
    # Perform classical MDS
    mds_result <- cmdscale(D, k = k, eig = TRUE)
    results$eigenvalues[k,] <- mds_result$eig
    
    # Traditional stress calculation
    if("stress" %in% methods) {
      pred_D <- as.matrix(Rfast::Dist(mds_result$points))
      results$stress_vals[k] <- sqrt(sum((D - pred_D)^2) / sum(D^2))
    }
    
    # Information criteria calculations
    if(any(c("AIC", "BIC") %in% methods)) {
      if(!"stress" %in% methods) {
        pred_D <- as.matrix(Rfast::Dist(mds_result$points))
      }
      log_lik <- compute_likelihood(D, pred_D, sigma)
      results$likelihood_vals[k] <- log_lik
      
      n_params <- n*k - k*(k+1)/2 - k
      
      if("AIC" %in% methods) {
        results$criteria[k, "AIC"] <- -2 * log_lik + 2 * n_params
      }
      if("BIC" %in% methods) {
        results$criteria[k, "BIC"] <- -2 * log_lik + log(n*(n-1)/2) * n_params
      }
    }
    
    # Bootstrap calculations
    if("bootstrap" %in% methods) {
      results <- compute_bootstrap(results, D, k, n_bootstrap, verbose)
    }
    
    # Cross-validation
    if("crossval" %in% methods) {
      results <- compute_crossval(results, D, k, n_folds, verbose)
    }
  }
  
  return(results)
}

# Compute bootstrap iterations
compute_bootstrap <- function(results, D, k, n_bootstrap, verbose) {
  n <- nrow(D)
  
  if(verbose && k == 1) {
    cli::cli_alert_info("Performing bootstrap iterations for dimension {k}...")
    pb_boot <- progress::progress_bar$new(
      format = "Bootstrap [:bar] :percent eta: :eta",
      total = n_bootstrap,
      clear = FALSE
    )
  }
  
  for(b in 1:n_bootstrap) {
    if(verbose) pb_boot$tick()
    
    boot_idx <- sample(n, replace = TRUE)
    boot_D <- D[boot_idx, boot_idx]
    boot_mds <- cmdscale(boot_D, k = k)
    boot_pred_D <- as.matrix(Rfast::Dist(boot_mds))
    results$bootstrap_stress[k, b] <- sqrt(sum((boot_D - boot_pred_D)^2) / sum(boot_D^2))
  }
  
  return(results)
}

# Compute cross-validation
compute_crossval <- function(results, D, k, n_folds, verbose) {
  n <- nrow(D)
  
  if(verbose && k == 1) {
    cli::cli_alert_info("Performing cross-validation for dimension {k}...")
    pb_cv <- progress::progress_bar$new(
      format = "Cross-validation [:bar] :percent eta: :eta",
      total = n_folds,
      clear = FALSE
    )
  }
  
  fold_indices <- split(sample(1:n), rep(1:n_folds, length.out = n))
  
  for(fold in 1:n_folds) {
    if(verbose) pb_cv$tick()
    
    test_idx <- fold_indices[[fold]]
    train_idx <- unlist(fold_indices[-fold])
    
    D_train <- D[train_idx, train_idx]
    D_test <- D[test_idx, train_idx]
    
    mds_train <- cmdscale(D_train, k = k)
    mds_test <- cmdscale(D_test, k = k)
    
    pred_D <- as.matrix(Rfast::Dist(rbind(mds_train, mds_test)))
    test_distances <- pred_D[nrow(mds_train) + 1:length(test_idx), 1:nrow(mds_train)]
    results$cv_error[k, fold] <- mean((D_test - test_distances)^2)
  }
  
  return(results)
}

# Post-process results
post_process_results <- function(results, methods, max_dim) {
  results$optimal_dimensions <- list()
  
  if("AIC" %in% methods) {
    results$optimal_dimensions$AIC <- which.min(results$criteria[, "AIC"])
  }
  if("BIC" %in% methods) {
    results$optimal_dimensions$BIC <- which.min(results$criteria[, "BIC"])
  }
  if("bootstrap" %in% methods) {
    boot_ci <- t(apply(results$bootstrap_stress, 1, quantile, probs = c(0.025, 0.975)))
    results$bootstrap_results <- list(
      mean_stress = rowMeans(results$bootstrap_stress),
      confidence_intervals = boot_ci
    )
    # Find dimension where bootstrap stress stabilizes
    stress_diff <- diff(rowMeans(results$bootstrap_stress))
    results$optimal_dimensions$bootstrap <- which(abs(stress_diff) < mean(abs(stress_diff))/10)[1]
  }
  if("crossval" %in% methods) {
    cv_means <- rowMeans(results$cv_error)
    results$optimal_dimensions$crossval <- which.min(cv_means)
  }
  if("stress" %in% methods) {
    # Find elbow point in stress curve
    stress_diff <- diff(results$stress_vals)
    stress_diff2 <- diff(stress_diff)
    results$optimal_dimensions$stress <- which.max(abs(stress_diff2)) + 1
  }
  
  return(results)
}

# Create visualizations
create_visualizations <- function(results, methods, max_dim) {
  cli::cli_h2("Creating visualizations")
  
  n_plots <- sum(c(
    any(c("AIC", "BIC") %in% methods),
    TRUE,  # eigenvalues always included
    "bootstrap" %in% methods,
    "crossval" %in% methods,
    "stress" %in% methods
  ))
  
  par(mfrow = c(ceiling(n_plots/2), min(2, n_plots)))
  
  # Information Criteria plot
  if(any(c("AIC", "BIC") %in% methods)) {
    plot_information_criteria(results, methods, max_dim)
  }
  
  # Eigenvalue scree plot
  plot_eigenvalues(results)
  
  # Traditional stress scree plot
  if("stress" %in% methods) {
    plot_stress_scree(results, max_dim)
  }
  
  # Bootstrap plot
  if("bootstrap" %in% methods) {
    plot_bootstrap_results(results, max_dim)
  }
  
  # Cross-validation plot
  if("crossval" %in% methods) {
    plot_crossval_results(results, max_dim)
  }
  
  par(mfrow = c(1, 1))
}

# Individual plotting functions
plot_information_criteria <- function(results, methods, max_dim) {
  plot_data <- scale(results$criteria[, methods[methods %in% c("AIC", "BIC")], drop = FALSE])
  matplot(1:max_dim, plot_data, type = "l", lty = 1:2, col = 1:2,
          xlab = "Number of dimensions", ylab = "Standardized criterion value",
          main = "Model Selection Criteria")
  legend("topright", methods[methods %in% c("AIC", "BIC")], 
         lty = 1:2, col = 1:2)
}

plot_eigenvalues <- function(results) {
  n <- ncol(results$eigenvalues)
  plot(1:n, results$eigenvalues[1,], type = "l",
       xlab = "Component", ylab = "Eigenvalue",
       main = "Eigenvalue Scree Plot")
  abline(h = 0, col = "gray", lty = 2)
}

plot_stress_scree <- function(results, max_dim) {
  plot(1:max_dim, results$stress_vals, type = "b",
       xlab = "Number of dimensions", ylab = "Stress",
       main = "Stress Scree Plot")
  abline(v = results$optimal_dimensions$stress, col = "red", lty = 2)
  text(results$optimal_dimensions$stress, mean(results$stress_vals),
       paste("Elbow =", results$optimal_dimensions$stress),
       pos = 4)
}

plot_bootstrap_results <- function(results, max_dim) {
  boot_ci <- results$bootstrap_results$confidence_intervals
  plot(1:max_dim, results$bootstrap_results$mean_stress, type = "l",
       xlab = "Number of dimensions", ylab = "Stress",
       main = "Bootstrap Stress Assessment",
       ylim = range(boot_ci))
  polygon(c(1:max_dim, max_dim:1),
          c(boot_ci[,1], rev(boot_ci[,2])),
          col = rgb(0,0,0,0.2), border = NA)
}

plot_crossval_results <- function(results, max_dim) {
  cv_means <- rowMeans(results$cv_error)
  cv_sd <- apply(results$cv_error, 1, sd)
  plot(1:max_dim, cv_means, type = "l",
       xlab = "Number of dimensions", ylab = "CV Error",
       main = "Cross-validation Error",
       ylim = range(c(cv_means + cv_sd, cv_means - cv_sd)))
  polygon(c(1:max_dim, max_dim:1),
          c(cv_means + cv_sd, rev(cv_means - cv_sd)),
          col = rgb(0,0,0,0.2), border = NA)
}

# Print summary
print_summary <- function(results, methods, sigma) {
  cli::cli_h2("Summary of Results")
  cli::cli_text("")
  cli::cli_text("Fixed sigma used for likelihood: {.val {sprintf('%.4f', sigma)}}")
  cli::cli_text("")
  cli::cli_text("Optimal dimensions by method:")
  for(method in methods) {
    cli::cli_text("{.field {method}}: {.val {results$optimal_dimensions[[method]]}}")
  }
}

#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_dat_all
#' @return
#' @author rdinnager
#' @export
make_phate_analysis <- function(trophic_niche_dat_all, niche_pal) {

  dat <- trophic_niche_dat_all |>
    pull(codes)
  
  init_phate <- phate(dat)
  plot(init_phate$embedding)
  #diff_op <- init_phate$operator$diff_op
  #diff_pot <- init_phate$operator$diff_potential
  
  #test <- model_based_mds_selection(diff_pot, methods = c("AIC", "BIC", "stress"), n_bootstrap = 25, n_folds = 5, verbose = TRUE, max_dim = 25)

  phate_results <- phate(dat, ndim = 5)
  phate_pca <- prcomp(phate_results$embedding)
  phate_pca$x
  
  names(niche_pal) <- levels(as.factor(trophic_niche_dat_all$Trophic.Niche))
  
  rgl::points3d(phate_pca$x, col = niche_pal[trophic_niche_dat_all$Trophic.Niche], 
                size = 10)
  rgl::points3d(phate_pca$x[ , c(1, 2, 4)], col = niche_pal[trophic_niche_dat_all$Trophic.Niche], 
                size = 10)
  rgl::points3d(phate_pca$x[ , c(2, 3, 4)], col = niche_pal[trophic_niche_dat_all$Trophic.Niche], 
                size = 10)
  rgl::points3d(phate_pca$x[ , c(3, 4, 5)], col = niche_pal[trophic_niche_dat_all$Trophic.Niche], 
                size = 10)
  
}
