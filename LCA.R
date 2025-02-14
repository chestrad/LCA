library(poLCA)
library(data.table)

determine_optimal_lca <- function(data, var_names, min_class=2, max_class=5, 
                                  n_rep=10, maxiter=1000, sample_size=10000,
                                  n_cv_folds=5) {
  # Input validation
  if (length(var_names) < 2) stop("At least two variables required for LCA")
  if (!all(var_names %in% names(data))) {
    missing_vars <- var_names[!var_names %in% names(data)]
    stop(sprintf("Variables not found in data: %s", paste(missing_vars, collapse=", ")))
  }
  
  # Convert to data.table and ensure binary coding
  analysis_data <- as.data.table(data[, var_names])
  
  for(col in var_names) {
    if(!all(analysis_data[[col]] %in% c(0,1,NA))) {
      stop(sprintf("Column %s contains non-binary values", col))
    }
    # Add 1 to convert from 0/1 to 1/2 coding required by poLCA
    analysis_data[[col]] <- analysis_data[[col]] + 1
  }
  
  # Create CV folds
  set.seed(42)
  fold_indices <- sample(rep(1:n_cv_folds, length.out=nrow(analysis_data)))
  
  # Initialize results storage
  results <- data.table(
    n_class = rep(min_class:max_class, each=n_cv_folds),
    fold = rep(1:n_cv_folds, times=length(min_class:max_class)),
    BIC = NA_real_,
    AIC = NA_real_,
    entropy = NA_real_,
    convergence_errors = FALSE,
    loglik = NA_real_
  )
  
  # Prepare formula
  formula_str <- paste("cbind(", paste(var_names, collapse=","), ") ~ 1")
  formula <- as.formula(formula_str)
  
  # Process each fold and class size
  for(k in min_class:max_class) {
    cat(sprintf("\nAnalyzing %d classes\n", k))
    
    for(fold in 1:n_cv_folds) {
      # Get training indices for this fold
      train_idx <- which(fold_indices != fold)
      if(length(train_idx) > sample_size) {
        train_idx <- sample(train_idx, sample_size)
      }
      
      # Get training data
      train_data <- analysis_data[train_idx]
      
      # Try multiple starts
      best_bic <- Inf
      best_fit <- NULL
      
      
      for(i in 1:n_rep) {
        cat(sprintf("Fold %d, Start %d/%d\n", fold, i, n_rep))
        
        fit <- tryCatch({
          # Set random seed for reproducibility
          set.seed(42 * k * fold * i)
          
          lca_fit <- poLCA(formula, 
                           train_data, 
                           nclass=k, 
                           maxiter=maxiter,
                           nrep=1,
                           verbose=FALSE)
           
          if(!is.null(lca_fit) && !lca_fit$eflag && !is.na(lca_fit$bic) && 
             lca_fit$bic < best_bic) {
            best_bic <- lca_fit$bic
            best_fit <- lca_fit
            cat("  Converged! BIC:", lca_fit$bic, "\n")
          }
          lca_fit
        }, error = function(e) {
          cat("  Error:", conditionMessage(e), "\n")
          NULL
        })
        
        gc()
      }
       
      
      # Store results from best model
      if(!is.null(best_fit)) {
        idx <- which(results$n_class == k & results$fold == fold)
        results[idx, `:=`(
          BIC = best_fit$bic,
          AIC = best_fit$aic,
          entropy = -sum(best_fit$posterior * log(best_fit$posterior + 1e-10), 
                         na.rm=TRUE) / nrow(best_fit$posterior),
          convergence_errors = best_fit$eflag,
          loglik = best_fit$llik
        )]
        
        cat(sprintf("Class %d, Fold %d Results:\n", k, fold))
        print(results[idx])
      }
      
      rm(train_data, best_fit)
      gc()
    }
  }
  
  # Calculate average metrics across folds
  cv_results <- results[, .(
    mean_BIC = mean(BIC, na.rm=TRUE),
    sd_BIC = sd(BIC, na.rm=TRUE),
    error_rate = mean(convergence_errors, na.rm=TRUE)
  ), by=n_class]
  
  # Select optimal model
  if(any(!is.na(cv_results$mean_BIC))) {
    optimal_k <- cv_results[which.min(mean_BIC)]$n_class
    cat(sprintf("\nOptimal number of classes based on CV-BIC: %d\n", optimal_k))
  } else {
    warning("No valid models found. Using minimum number of classes.")
    optimal_k <- min_class
  }
  
  return(list(
    optimal_k = optimal_k,
    results = results,
    cv_results = cv_results
  ))
}

apply_lca <- function(data, var_names, k, maxiter=1000, n_rep=10) {
  # Convert to data.table and ensure binary coding
  analysis_data <- as.data.table(data[, var_names])
  for(col in var_names) {
    if(!all(analysis_data[[col]] %in% c(0,1,NA))) {
      stop(sprintf("Column %s contains non-binary values", col))
    }
    analysis_data[[col]] <- analysis_data[[col]] + 1
  }
  
  # Prepare formula
  formula_str <- paste("cbind(", paste(var_names, collapse=","), ") ~ 1")
  formula <- as.formula(formula_str)
  
  # Try multiple starts
  best_bic <- Inf
  best_fit <- NULL
  
  for(i in 1:n_rep) {
    tryCatch({
      cat(sprintf("Attempt %d/%d\n", i, n_rep))
      
      lca_fit <- poLCA(formula, 
                       analysis_data, 
                       nclass=k, 
                       maxiter=maxiter,
                       nrep=1,
                       verbose=FALSE)
      
      if(!is.null(lca_fit) && !lca_fit$eflag && !is.na(lca_fit$bic) && 
         lca_fit$bic < best_bic) {
        best_bic <- lca_fit$bic
        best_fit <- lca_fit
      }
      
    }, error = function(e) {
      cat(sprintf("Warning in attempt %d: %s\n", i, e$message))
    })
    
    gc()
  }
  
  if(is.null(best_fit)) {
    stop("Failed to fit valid model after all attempts")
  }
   
  pattern_freq <- analysis_data[,
                                .(Pattern = {
                                  binary_values <- lapply(.SD, function(x) x-1)
                                  pattern <- do.call(paste0, binary_values)
                                  sprintf("%06d", as.numeric(pattern))
                                }),
                                .SDcols = var_names
  ][, .N, by = Pattern]
   
  pattern_summary <- data.frame( 
    Pattern = format(pattern_freq$Pattern, justify="right", width=6),
    Frequency = pattern_freq$N,
    Proportion = pattern_freq$N / nrow(analysis_data),
    stringsAsFactors = FALSE   
  )
  
  results <- list(
    class_assignments = best_fit$predclass,
    class_sizes = table(best_fit$predclass),
    model_fit = data.table(
      metric = c("BIC", "AIC", "LogLik", "df", "convergence_errors", "n_iter"),
      value = c(best_fit$bic, best_fit$aic, best_fit$llik, 
                best_fit$resid.df, best_fit$eflag, best_fit$numiter)
    ),
    class_probs = setNames(
      lapply(best_fit$probs, function(x) x[,2]),
      var_names
    ),
    posterior_probs = best_fit$posterior,
    pattern_summary = pattern_summary  #REVISED 
  )
  
  return(results)
}


##### USAGE #####
### Cross-validation to obtain optimal k for LCA
results <- determine_optimal_lca(data, 
                                 var_names=c("S1", "S2", "S3", "S4", "S5", "S6"),
                                 sample_size=34041,   
                                 n_cv_folds=5,      
                                 n_rep=10,          
                                 maxiter=2000)     
write.csv(data.frame(optimal_k = results$optimal_k), "optimal_k.csv", row.names = FALSE) 
write.csv(results$results, "cv_results.csv", row.names = FALSE) 
write.csv(results$cv_results, "cv_results_summary.csv", row.names = FALSE)

### Fit final model with optimal number of classes
final_model <- apply_lca(data, 
                         var_names=c("S1", "S2", "S3", "S4", "S5", "S6"),
                         k=results$optimal_k,
                         maxiter=2000,  
                         n_rep=10) 

# Save class sizes
write.csv(data.frame(class = names(final_model$class_sizes),
                     size = as.numeric(final_model$class_sizes)), 
          "class_sizes.csv", row.names = FALSE)

# Save model fit statistics
write.csv(final_model$model_fit, "model_fit.csv", row.names = FALSE)

# Save class probabilities
n_classes <- length(final_model$class_sizes)  # Get number of classes from class_sizes
class_probs_df <- data.frame(
  variable = rep(names(final_model$class_probs), each = n_classes),
  class = rep(paste0("class", 1:n_classes), length(final_model$class_probs)),
  probability = unlist(final_model$class_probs)
)
write.csv(class_probs_df, "class_probabilities.csv", row.names = FALSE)
 
# class assignment, then perform survival analyses using Class variable
data$Class <- factor(final_model$class_assignments) 

# Save co-occurrence pattern summary 
final_model$pattern_summary$Pattern <- paste0("=\"", final_model$pattern_summary$Pattern, "\"")
write.csv(final_model$pattern_summary, "pattern_frequency.csv", row.names = FALSE)
 
