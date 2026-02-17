\name{pdi_var}
\alias{pdi_var}
\title{Calculate PDI Value with Variance Estimation}
\usage{
pdi_var(y, d, method="multinom", B=250, level=0.95, \dots)
}
\description{
Compute the Polytomous Discrimination Index (PDI) value with bootstrap variance estimation and confidence intervals for two or three or four categories classifiers. This function implements the variance estimation method described in Dover et al. (2021).
}
\arguments{
  \item{y}{The multinomial response vector with two, three or four categories. It can be factor or integer-valued.}
  \item{d}{The set of candidate markers, including one or more columns. Can be a data frame or a matrix; if the method is "prob", then d should be the probability matrix.}
  \item{method}{Specifies what method is used to construct the classifier based on the marker set in d. Available option includes the following methods:"multinom": Multinomial Logistic Regression which is the default method, requiring R package nnet;"tree": Classification Tree method, requiring R package rpart;"svm": Support Vector Machine (C-classification and radial basis as default), requiring R package e1071;"lda": Linear Discriminant Analysis, requiring R package lda;"prob": d is a risk matrix resulted from any external classification algorithm obtained by the user.}
  \item{B}{Number of bootstrap samples. Default is 250.}
  \item{level}{The confidence level. Default value is 0.95.}
  \item{\dots}{Additional arguments in the chosen method's function.}
}
\details{
The function returns the PDI value with standard error and confidence interval based on bootstrap resampling. This implements the variance estimation method from Dover et al. (2021) for the polytomous discrimination index.

The PDI measures the discrimination ability of a classification model by calculating the probability that a randomly selected subject from each outcome category is correctly ranked by the model's predicted probabilities.
}
\value{
Returns an object of class "mcca.pdi.var". The PDI value with variance estimation using a particular learning method on a set of marker(s).

An object of class "mcca.pdi.var" is a list containing at least the following components:
  \item{call}{the matched call.}
  \item{measure}{the overall PDI value.}
  \item{se}{the bootstrap standard error of the overall PDI.}
  \item{ci}{the confidence interval for the overall PDI.}
  \item{level}{the confidence level used.}
  \item{B}{the number of bootstrap samples used.}
  \item{table}{a data frame with category-specific PDI values, standard errors, and confidence intervals.}
}
\references{
Dover, D.C., Savu, A., Engel, B. (2021). Polytomous discrimination index: Estimation and inference. Statistics in Medicine.

Li, J., Gao, M., D'Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.

Van Calster B, Vergouwe Y, Looman CWN, Van Belle V, Timmerman D and Steyerberg EW. Assessing the discriminative ability of risk models for more than two outcome categories. European Journal of Epidemiology 2012; 27: 761-770.

Li, J., Feng, Q., Fine, J.P., Pencina, M.J., Van Calster, B. (2018). Nonparametric estimation and inference for polytomous discrimination index. Statistical Methods in Medical Research. 27(10): 3092-3103.
}
\author{
Ming Gao: gaoming@umich.edu

Jialiang Li: stalj@nus.edu.sg

Variance estimation method by Anamaria Savu, PhD (Canadian VIGOUR Centre, University of Alberta)
}
\note{
Users are advised to change the operating settings of various classifiers since it is well known that machine learning methods require extensive tuning. Currently only some common and intuitive options are set as default and they are by no means the optimal parameterization for a particular data analysis. Users can put machine learning methods' parameters after tuning. A more flexible evaluation is to consider "method=prob" in which case the input d should be a matrix of membership probabilities with k columns and each row of d should sum to one.

The number of bootstrap samples (B) controls the precision of the variance estimate. Larger values give more stable estimates but take longer to compute.
}
\seealso{
\code{\link{pdi}}, \code{\link{ests}}
}
\examples{
str(iris)
data <- iris[, 1:4]
label <- iris[, 5]

# PDI with variance estimation using 50 bootstrap samples
pdi_var(y = label, d = data, method = "multinom", B = 50, trace = FALSE)

## Call:
## pdi_var(y = label, d = data, method = "multinom", B = 50, trace = FALSE)

## Overall Polytomous Discrimination Index:
##  0.9848

## Standard Error:
##  0.0089

## 95\% Confidence Interval:
##  [0.9674, 1.0000]

## Bootstrap Samples: 50

## Category-specific Polytomous Discrimination Index:
##  CATEGORIES VALUES     SE LOWER_CI UPPER_CI
##      setosa 1.0000 0.0000   1.0000   1.0000
##  versicolor 0.9772 0.0133   0.9511   1.0000
##   virginica 0.9772 0.0133   0.9511   1.0000


# Using tree method
pdi_var(y = label, d = data, method = "tree", B = 50)

# Using probability matrix directly
require(nnet)
fit <- multinom(label ~ ., data = data.frame(label = label, data),
                maxit = 1000, MaxNWts = 2000, trace = FALSE)
predict.probs <- predict(fit, type = "probs")
pp <- data.frame(predict.probs)
pdi_var(y = label, d = pp, method = "prob", B = 50)

}
\keyword{PDI}
\keyword{variance}
