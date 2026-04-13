pdi=function(y,d,method="multinom", withTies=TRUE,
             Var=FALSE, level=0.95, ...){

  y=factor(y)
  y_levels=levels(y)
  y=as.numeric(y)
  d=data.matrix(d)
  num=length(unique(y))

  pp=pm(y=y,d=d,method=method,...)
  pp=data.matrix(pp)
  
  #########################################
  if(Var){
    data_fast <- data.frame(
      outcome = factor(y, levels = seq_len(num))
    )
    data_fast <- cbind(data_fast, as.data.frame(pp))
    colnames(data_fast) <- c("outcome", paste0("p", seq_len(num)))
    
    var_res <- .pdi_var_efficient(
      data = data_fast,
      withTies = withTies
    )
    
    z_value <- qnorm(1 - (1 - level) / 2)
    
    pdi_overall <- var_res$pdi
    se_overall <- sqrt(var_res$sigma2 / nrow(data_fast))
    ci_overall <- c(
      max(0, pdi_overall - z_value * se_overall),
      min(1, pdi_overall + z_value * se_overall)
    )
    
    pdi_components <- var_res$pdi_categ
    se_components <- sqrt(var_res$sigma2_categ / nrow(data_fast))
    lower_components <- pmax(0, pdi_components - z_value * se_components)
    upper_components <- pmin(1, pdi_components + z_value * se_components)
    
    df <- data.frame(
      CATEGORIES = y_levels,
      VALUES = pdi_components,
      SE = se_components,
      LOWER_CI = lower_components,
      UPPER_CI = upper_components
    )
    
    result <- list(
      call = match.call(),
      measure = pdi_overall,
      se = se_overall,
      ci = ci_overall,
      level = level,
      table = df
    )
    
    class(result) <- "mcca.pdi.var"  
    return(result)
  }

  #########################################
  if(num==3 && !withTies){
    n1=which(y==1) #return the label
    n2=which(y==2)
    n3=which(y==3)
    pp1=pp[n1,]
    pp2=pp[n2,]
    pp3=pp[n3,]
    pdi1=0
    pdi2=0
    pdi3=0
    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pp1[i,1]>pp2[,1])*sum(pp1[i,1]>pp3[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pp2[i,2]>pp1[,2])*sum(pp2[i,2]>pp3[,2])
    }
    for(i in 1:length(n3)){
      pdi3=pdi3+sum(pp3[i,3]>pp1[,3])*sum(pp3[i,3]>pp2[,3])
    }
    pdi=(pdi1+pdi2+pdi3)/(3*length(n1)*length(n2)*length(n3))

    pdis=c(pdi1,pdi2,pdi3)/(length(n1)*length(n2)*length(n3))
    df=data.frame(CATEGORIES=sapply(1:num, function(i) y_levels[i]),VALUES=pdis)
    result=list(call=match.call(),measure=pdi,table=df)
    class(result)="mcca.pdi"
    return(result)

  #########################################
  }else if(num==4 && !withTies){
    n1=which(y==1) #return the label
    n2=which(y==2)
    n3=which(y==3)
    n4=which(y==4)
    pp1=pp[n1,]
    pp2=pp[n2,]
    pp3=pp[n3,]
    pp4=pp[n4,]
    pdi1=0
    pdi2=0
    pdi3=0
    pdi4=0
    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pp1[i,1]>pp2[,1])*sum(pp1[i,1]>pp3[,1])*sum(pp1[i,1]>pp4[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pp2[i,2]>pp1[,2])*sum(pp2[i,2]>pp3[,2])*sum(pp2[i,2]>pp4[,2])
    }
    for(i in 1:length(n3)){
      pdi3=pdi3+sum(pp3[i,3]>pp1[,3])*sum(pp3[i,3]>pp2[,3])*sum(pp3[i,3]>pp4[,3])
    }
    for(i in 1:length(n4)){
      pdi4=pdi4+sum(pp4[i,4]>pp1[,4])*sum(pp4[i,4]>pp2[,4])*sum(pp4[i,4]>pp3[,4])
    }
    pdi=(pdi1+pdi2+pdi3+pdi4)/(4*length(n1)*length(n2)*length(n3)*length(n4))

    pdis=c(pdi1,pdi2,pdi3,pdi4)/(length(n1)*length(n2)*length(n3)*length(n4))

    df=data.frame(CATEGORIES=sapply(1:num, function(i) y_levels[i]),VALUES=pdis)
    result=list(call=match.call(),measure=pdi,table=df)
    class(result)="mcca.pdi"
    return(result)

    #########################################

  }else if(num==2 && !withTies){
    n1=which(y==1) #return the label
    n2=which(y==2)
    pp1=pp[n1,]
    pp2=pp[n2,]

    pdi1=0
    pdi2=0

    for(i in 1:length(n1)){
      pdi1=pdi1+sum(pp1[i,1]>pp2[,1])+sum(pp1[i,1]==pp2[,1])
    }
    for(i in 1:length(n2)){
      pdi2=pdi2+sum(pp2[i,2]>pp1[,2])
    }

    pdi=(pdi1+pdi2)/(2*length(n1)*length(n2))

    pdis=c(pdi1,pdi2)/(length(n1)*length(n2))
    df=data.frame(CATEGORIES=sapply(1:num, function(i) y_levels[i]),VALUES=pdis)
    result=list(call=match.call(),measure=pdi,table=df)
    class(result)="mcca.pdi"
    return(result)
    
    #########################################
  }else{
    data_dover <- data.frame(outcome = y, pp)
    colnames(data_dover) <- c("outcome", paste0("p", 1:num))
    
    est <- .pdiest_dover(data_dover)
    
    pdi_overall <- as.numeric(est[1, 1])
    pdi_components <- as.numeric(est[-1, 1])
    
    df <- data.frame(
      CATEGORIES = sapply(1:num, function(i) y_levels[i]),
      VALUES = pdi_components
    )
    
    result <- list(
      call = match.call(),
      measure = pdi_overall,
      table = df
    )
    class(result) <- "mcca.pdi"
    return(result)
  }
  
}




###########################################################
# PDI Estimation in general cases using Dover et al. (2021) method
# Reference: Dover DC, Savu A, Engel B.
# Statistics in Medicine, 2021.
###########################################################
# Original implementation by:
# Anamaria Savu, PhD
# Canadian VIGOUR Centre
# University of Alberta
# savu@ualberta.ca
# November 2020
###########################################################

# Internal function: Estimates PDI and its components from probability matrix
# data: data.frame with $outcome and probability columns ($p1, $p2, ...)
#
# Usage notes:
# data must have one of its components named $outcome for storing
# the outcome values. The values of the outcome must cover all integer
# values between a lower and an upper bound. The remaining components
# of the data frame ($p1, $p2, ...) must record the estimated probabilities
# in a specific order: $p1 probabilities for the lowest value of the outcome,
# $p2 probabilities for the second lowest value of the outcome, and so on.
###########################################################
.pdiest_dover <- function(data) {
  y_raw <- data$outcome
  ymin <- min(y_raw)
  ymax <- max(y_raw)
  y_levels <- as.character(ymin:ymax)
  y <- factor(y_raw, levels = ymin:ymax)
  
  noutcome <- ymax - ymin
  p <- prod(table(y))
  pdi <- NULL
  
  for (i in seq_len(noutcome + 1)) {
    predprob <- data[[i + 1]]
    
    t0 <- table(predprob, y)
    t0 <- as.matrix(t0)
    
    if (!all(y_levels %in% colnames(t0))) {
      missing_cols <- setdiff(y_levels, colnames(t0))
      add_cols <- matrix(
        0,
        nrow = nrow(t0),
        ncol = length(missing_cols),
        dimnames = list(rownames(t0), missing_cols)
      )
      t0 <- cbind(t0, add_cols)
    }
    
    t0 <- t0[, y_levels, drop = FALSE]
    
    t <- cbind(
      t0[, i, drop = FALSE],
      t0[, -i, drop = FALSE]
    )
    
    restrictt <- if (noutcome == 1) {
      matrix(t[, 2, drop = TRUE], ncol = 1)
    } else {
      t[, 2:(noutcome + 1), drop = FALSE]
    }
    
    cmat <- apply(restrictt, 2, cumsum)
    if (is.null(dim(cmat))) {
      cmat <- matrix(cmat, ncol = noutcome)
    }
    
    if (nrow(cmat) == 1) {
      cnew <- matrix(0, nrow = 1, ncol = noutcome)
    } else {
      cnew <- rbind(
        matrix(0, nrow = 1, ncol = noutcome),
        cmat[seq_len(nrow(cmat) - 1), , drop = FALSE]
      )
    }
    
    mat <- NULL
    for (j in seq_len(noutcome)) {
      mat0 <- cbind(mat, 0)
      mat1 <- cbind(mat, 1)
      mat <- rbind(mat0, mat1)
    }
    
    r <- 0
    for (k in seq_len(nrow(mat))) {
      dt <- t(apply(restrictt, 1, function(x) mat[k, ] * x))
      dcnew <- t(apply(cnew, 1, function(x) (1 - mat[k, ]) * x))
      
      if (is.null(dim(dt))) {
        dt <- matrix(dt, ncol = noutcome)
      }
      if (is.null(dim(dcnew))) {
        dcnew <- matrix(dcnew, ncol = noutcome)
      }
      
      if (noutcome == 1) {
        dfinal <- cbind(
          t[, 1, drop = FALSE],
          matrix(dt + dcnew, ncol = 1)
        )
      } else {
        dfinal <- cbind(
          t[, 1, drop = FALSE],
          dt + dcnew
        )
      }
      
      r <- r + sum(apply(dfinal, 1, prod)) / (1 + sum(mat[k, ]))
    }
    
    r <- r / p
    pdi <- rbind(pdi, r)
  }
  
  pdi <- rbind(mean(pdi), pdi)
  rownames(pdi) <- NULL
  pdi
}


###########################################################
# PDI with efficient variance estimation
# Reference: Feng Q, Jia T, Liu P, Van Calster B, Li J.
# Statistics in Medicine, 2026.
###########################################################

# Internal function: Estimates PDI, its variance, standard deviation, and 
# confidence interval from probability matrix data: data.frame with $outcome 
# and probability columns ($p1, $p2, ...)
#
# Usage notes:
# data must have one of its components named $outcome for storing
# the outcome values. The values of the outcome must cover all integer
# values between a lower and an upper bound. The remaining components
# of the data frame ($p1, $p2, ...) must record the estimated probabilities
# in a specific order: $p1 probabilities for the lowest value of the outcome,
# $p2 probabilities for the second lowest value of the outcome, and so on.
###########################################################

.pdi_var_efficient <- function(data, withTies = TRUE) {
  
  y <- data$outcome
  y_char <- as.character(y)
  n <- length(y)
  categ_table <- table(y)
  categ <- names(categ_table)
  M <- length(categ)
  sz_prod <- prod(categ_table)
  
  if (withTies) {
    
    matC_init <- c()
    for (j in 1:(M - 1)) {
      matC0 <- cbind(matC_init, 0)
      matC1 <- cbind(matC_init, 1)
      matC_init <- rbind(matC0, matC1)
    }
    matC_init <- cbind(1, matC_init)
    
    pdi_categ <- c()
    pdi_sub <- c()
    
    for (a in 1:M) {
      
      predprob <- data[, a + 1]
      nn <- table(predprob, y)
      nn <- as.matrix(nn)
      
      NN_temp <- apply(nn, 2, cumsum)
      
      if (is.null(dim(NN_temp))) {
        NN_temp <- matrix(NN_temp, nrow = nrow(nn), ncol = ncol(nn))
        colnames(NN_temp) <- colnames(nn)
        rownames(NN_temp) <- rownames(nn)
      }
      
      if (nrow(NN_temp) == 1) {
        NN <- matrix(0, nrow = 1, ncol = M)
        colnames(NN) <- colnames(NN_temp)
        rownames(NN) <- rownames(NN_temp)
      } else {
        NN <- rbind(
          matrix(0, nrow = 1, ncol = M),
          NN_temp[1:(nrow(NN_temp) - 1), , drop = FALSE]
        )
        colnames(NN) <- colnames(NN_temp)
        rownames(NN) <- rownames(NN_temp)
      }
      
      aid <- which(y == levels(y)[a])
      S_a <- unique(predprob[aid])
      S_card <- length(S_a)
      
      nn_a <- nn[as.character(S_a), , drop = FALSE]
      NN_a <- NN[as.character(S_a), , drop = FALSE]
      
      matC <- matC_init
      matC[, 1] <- matC[, a]
      matC[, a] <- 1
      
      r <- 0
      A <- matrix(0, nrow = S_card, ncol = n)
      
      for (k in 1:nrow(matC)) {
        
        cc <- matC[k, ]
        names(cc) <- categ
        
        cc_mat <- matrix(cc, nrow = nrow(nn_a), ncol = ncol(nn_a), byrow = TRUE)
        dfinal <- cc_mat * nn_a + (1 - cc_mat) * NN_a
        
        row_prod <- apply(dfinal, 1, prod)
        zero_count <- rowSums(dfinal == 0)
        dfinal_all <- matrix(0, nrow(dfinal), ncol(dfinal))
        
        idx0 <- which(zero_count == 0)
        if (length(idx0) > 0) {
          dfinal_all[idx0, ] <- row_prod[idx0] / dfinal[idx0, ]
        }
        
        idx1 <- which(zero_count == 1)
        if (length(idx1) > 0) {
          for (i in idx1) {
            zero_idx <- which(dfinal[i, ] == 0)
            dfinal_all[i, zero_idx] <- prod(dfinal[i, -zero_idx])
          }
        }
        
        colnames(dfinal_all) <- categ
        dfinal_all <- dfinal_all[, y_char, drop = FALSE]
        
        matC_all <- cc[y_char]
        numerator <- t(sapply(S_a, function(pp) {
          matC_all * (predprob == pp) + (1 - matC_all) * (predprob < pp)
        }))
        
        A <- A + numerator * dfinal_all / sum(cc)
        r <- r + sum(apply(dfinal, 1, prod)) / sum(cc)
      }
      
      AA <- apply(A, 2, sum)
      pdi_sub <- cbind(pdi_sub, AA / sz_prod)
      pdi_categ <- c(pdi_categ, r / sz_prod)
    }
    
  } else {
    
    pdi_categ <- c()
    pdi_sub <- c()
    
    for (a in 1:M) {
      
      predprob <- data[, a + 1]
      nn <- table(predprob, y)
      nn <- as.matrix(nn)
      
      NN_temp <- apply(nn, 2, cumsum)
      
      if (is.null(dim(NN_temp))) {
        NN_temp <- matrix(NN_temp, nrow = nrow(nn), ncol = ncol(nn))
        colnames(NN_temp) <- colnames(nn)
        rownames(NN_temp) <- rownames(nn)
      }
      
      if (nrow(NN_temp) == 1) {
        NN <- matrix(0, nrow = 1, ncol = M)
        colnames(NN) <- colnames(NN_temp)
        rownames(NN) <- rownames(NN_temp)
      } else {
        NN <- rbind(
          matrix(0, nrow = 1, ncol = M),
          NN_temp[1:(nrow(NN_temp) - 1), , drop = FALSE]
        )
        colnames(NN) <- colnames(NN_temp)
        rownames(NN) <- rownames(NN_temp)
      }
      
      aid <- which(y == levels(y)[a])
      S_a <- unique(predprob[aid])
      
      nn_a <- nn[as.character(S_a), , drop = FALSE]
      NN_a <- NN[as.character(S_a), , drop = FALSE]
      
      dfinal <- NN_a
      dfinal[, a] <- 1
      
      row_prod <- apply(dfinal, 1, prod)
      zero_count <- rowSums(dfinal == 0)
      dfinal_all <- matrix(0, nrow(dfinal), ncol(dfinal))
      
      idx0 <- which(zero_count == 0)
      if (length(idx0) > 0) {
        dfinal_all[idx0, ] <- row_prod[idx0] / dfinal[idx0, ]
      }
      
      idx1 <- which(zero_count == 1)
      if (length(idx1) > 0) {
        for (i in idx1) {
          zero_idx <- which(dfinal[i, ] == 0)
          dfinal_all[i, zero_idx] <- prod(dfinal[i, -zero_idx])
        }
      }
      
      colnames(dfinal_all) <- categ
      dfinal_all <- dfinal_all[, y_char, drop = FALSE]
      
      numerator <- t(sapply(S_a, function(pp) {
        predprob < pp
      }))
      
      A <- numerator * dfinal_all
      r <- sum(apply(dfinal, 1, prod))
      
      AA <- apply(A, 2, sum)
      AA_a <- dfinal_all[, aid[1]][match(predprob[aid], S_a)]
      AA[aid] <- AA_a
      
      pdi_sub <- cbind(pdi_sub, AA / sz_prod)
      pdi_categ <- c(pdi_categ, r / sz_prod)
    }
  }
  
  names(pdi_categ) <- categ
  colnames(pdi_sub) <- categ
  
  pdi <- mean(pdi_categ)
  
  sigma2_term1 <- n * sum(apply(pdi_sub, 1, function(x) (mean(x))^2))
  sigma2_term2 <- n * pdi^2 * sum(1 / categ_table)
  sigma2 <- sigma2_term1 - sigma2_term2
  
  # category-specific asymptotic variance
  sigma2_categ <- sapply(1:M, function(a) {
    n * sum(pdi_sub[, a]^2) - n * pdi_categ[a]^2 * sum(1 / categ_table)
  })
  
  sigma2 <- max(as.numeric(sigma2), 0)
  sigma2_categ <- pmax(as.numeric(sigma2_categ), 0)
  
  list(
    pdi = as.numeric(pdi),
    sigma2 = sigma2,
    pdi_categ = as.numeric(pdi_categ),
    sigma2_categ = sigma2_categ,
    categories = categ
  )
}

