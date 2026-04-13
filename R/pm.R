pm <- function(y, d, method = "multinom", ...) {
  
  y_fac <- factor(y)
  y_levels <- levels(y_fac)
  y_num <- as.numeric(y_fac)
  target_levels <- as.character(seq_len(length(y_levels)))
  d <- data.matrix(d)
  
  if (method == "multinom") {
    fit <- nnet::multinom(y_fac ~ d, ...)
    predict.test.probs <- predict(fit, type = "probs")
    predict.test.df <- data.frame(predict.test.probs, check.names = FALSE)
    
    if (length(target_levels) == 2 && ncol(predict.test.df) == 1) {
      pp <- data.frame(
        "1" = 1 - predict.test.df[, 1],
        "2" = predict.test.df[, 1],
        check.names = FALSE
      )
    } else {
      old_names <- colnames(predict.test.df)
      idx <- match(old_names, y_levels)
      if (any(is.na(idx))) {
        stop("Unexpected class names returned by multinom prediction.")
      }
      colnames(predict.test.df) <- as.character(idx)
      pp <- .complete_prob_df(predict.test.df, target_levels)
    }
    
    return(pp)
  }
  
  if (method == "tree") {
    fit <- rpart::rpart(y_fac ~ d, ...)
    predict.test.probs <- predict(fit, type = "prob")
    predict.test.df <- data.frame(predict.test.probs, check.names = FALSE)
    
    old_names <- colnames(predict.test.df)
    idx <- match(old_names, y_levels)
    if (any(is.na(idx))) {
      stop("Unexpected class names returned by tree prediction.")
    }
    colnames(predict.test.df) <- as.character(idx)
    pp <- .complete_prob_df(predict.test.df, target_levels)
    
    return(pp)
  }
  
  if (method == "svm") {
    fit <- e1071::svm(y_fac ~ d, ..., probability = TRUE)
    predict.test <- predict(fit, d, probability = TRUE)
    predict.test <- attr(predict.test, "probabilities")
    predict.test.df <- data.frame(predict.test, check.names = FALSE)
    
    old_names <- colnames(predict.test.df)
    idx <- match(old_names, y_levels)
    if (any(is.na(idx))) {
      stop("Unexpected class names returned by svm prediction.")
    }
    colnames(predict.test.df) <- as.character(idx)
    pp <- .complete_prob_df(predict.test.df, target_levels)
    
    return(pp)
  }
  
  if (method == "lda") {
    fit <- MASS::lda(y_fac ~ d, ...)
    predict.test <- predict(fit)$posterior
    predict.test.df <- data.frame(predict.test, check.names = FALSE)
    
    old_names <- colnames(predict.test.df)
    idx <- match(old_names, y_levels)
    if (any(is.na(idx))) {
      stop("Unexpected class names returned by lda prediction.")
    }
    colnames(predict.test.df) <- as.character(idx)
    pp <- .complete_prob_df(predict.test.df, target_levels)
    
    return(pp)
  }
  
  if (method == "prob") {
    pp_sum <- rowSums(d)
    if (any(pp_sum < 0.999 | pp_sum > 1.001)) {
      cat("ERROR: The input value \"d\" should be a probability matrix.")
      return(NULL)
    }
    pp <- d
    return(pp)
  }
  
  stop("Unsupported method.")
}




.complete_prob_df <- function(prob, target_levels) {
  prob <- as.data.frame(prob, check.names = FALSE)
  
  out <- as.data.frame(
    matrix(0, nrow = nrow(prob), ncol = length(target_levels)),
    check.names = FALSE
  )
  colnames(out) <- target_levels
  
  common <- intersect(colnames(prob), target_levels)
  if (length(common) > 0) {
    out[, common] <- prob[, common, drop = FALSE]
  }
  
  out
}