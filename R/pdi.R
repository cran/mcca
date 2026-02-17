pdi=function(y,d,method="multinom",...){

  y=factor(y)
  y_levels=levels(y)
  y=as.numeric(y)
  d=data.matrix(d)
  num=length(unique(y))

  pp=pm(y=y,d=d,method=method,...)

  #########################################
  if(num==3){
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
  }else if(num==4){
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

  }else if(num==2){
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
  }
}


###########################################################
# PDI Variance Estimation using Dover et al. (2021) method
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
  y <- data$outcome
  ymin <- min(y)
  ymax <- max(y)
  noutcome <- ymax - ymin
  p <- prod(table(y))
  pdi <- c()

  for (i in 1:(noutcome + 1)) {
    predprob <- data[, (i + 1)]
    #READ predicted probabilities for level i
    t0 <- table(predprob, y)
    #CALCULATE frequencies of predicted probabilities for level i by outcome

    dim1 <- dim(t0)[1]
    dim2 <- dim(t0)[2]
    t <- cbind(t0[, i], t0[, -i])
    #REORDER columns
    restrictt <- if (noutcome == 1) {
      matrix(t[, 2:(noutcome + 1)], ncol = 1)
    } else {
      t[, 2:(noutcome + 1)]
    }
    #REMOVE first column of t

    c <- apply(restrictt, 2, cumsum)
    #CALCULATE cumulative frequencies of predicted probabilities for level i by outcome
    cnew <- if (noutcome == 1) {
      rbind(rep(0, noutcome), matrix(c[1:(dim(c)[1] - 1), ], ncol = 1))
    } else {
      rbind(rep(0, noutcome), c[1:(dim(c)[1] - 1), ])
    }
    #INTRODUCE a row of zeros at the begining of c

    mat <- c()
    #MATRIX of 0s and 1s of dimension 2^(noutcome) x noutcome
    for (j in 1:noutcome) {
      mat0 <- cbind(mat, 0)
      mat1 <- cbind(mat, 1)
      mat <- rbind(mat0, mat1)
    }

    r <- 0
    for (k in 1:dim(mat)[1]) {
      dt <- t(apply(restrictt, 1, function(x) mat[k, ] * x))
      dcnew <- t(apply(cnew, 1, function(x) (1 - mat[k, ]) * x))
      dfinal <- if (noutcome == 1) {
        cbind(t[, 1], t(dt + dcnew))
      } else {
        cbind(t[, 1], dt + dcnew)
      }
      #TAKE all combinations of frequencies and cumulative frequencies
      r <- r + sum(apply(dfinal, 1, prod)) / (1 + sum(mat[k, ]))
      #MULTIPLY across rows
    }

    r <- r / p
    #PDI component for outcome i
    pdi <- rbind(pdi, r)
  }
  pdi <- rbind(mean(pdi), pdi)
  pdi
}

# PDI with variance estimation using Dover et al. (2021) method
# y: outcome vector (factor or numeric)
# d: probability matrix or marker data
# method: classification method
# B: number of bootstrap samples (default 250)
# level: confidence level (default 0.95)
# ...: additional arguments passed to pm()
pdi_var <- function(y, d, method = "multinom", B = 250, level = 0.95, ...) {

  y <- factor(y)
  y_levels <- levels(y)
  y_numeric <- as.numeric(y)
  d <- data.matrix(d)
  num <- length(unique(y_numeric))

  pp <- pm(y = y_numeric, d = d, method = method, ...)
  pp <- data.matrix(pp)

  data_dover <- data.frame(outcome = y_numeric, pp)
  colnames(data_dover) <- c("outcome", paste0("p", 1:num))

  estimate <- .pdiest_dover(data_dover)

  samplesize <- nrow(data_dover)
  pdibs <- NULL

  for (i in 1:B) {
    vec <- sample.int(samplesize, size = samplesize, replace = TRUE)
    mydatabs <- data_dover[vec, ]

    if (length(unique(mydatabs$outcome)) < num) {
      next
    }

    bs_result <- .pdiest_dover(mydatabs)
    if (is.null(pdibs)) {
      pdibs <- bs_result
    } else {
      pdibs <- cbind(pdibs, bs_result)
    }
  }

  stderr <- sqrt(apply(pdibs, 1, var))
  z_value <- qnorm(1 - (1 - level) / 2)
  lowerci <- pmax(0, estimate - z_value * stderr)
  upperci <- pmin(1, estimate + z_value * stderr)

  pdi_overall <- estimate[1]
  pdi_components <- estimate[-1]
  se_overall <- stderr[1]
  se_components <- stderr[-1]
  ci_overall <- c(lowerci[1], upperci[1])
  ci_components <- cbind(lowerci[-1], upperci[-1])

  df <- data.frame(
    CATEGORIES = y_levels,
    VALUES = as.numeric(pdi_components),
    SE = as.numeric(se_components),
    LOWER_CI = as.numeric(ci_components[, 1]),
    UPPER_CI = as.numeric(ci_components[, 2])
  )

  result <- list(
    call = match.call(),
    measure = as.numeric(pdi_overall),
    se = as.numeric(se_overall),
    ci = ci_overall,
    level = level,
    B = B,
    table = df
  )

  class(result) <- "mcca.pdi.var"
  return(result)
}
