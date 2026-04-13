\name{pdi}
\alias{pdi}
\title{Calculate PDI Value}
\usage{
pdi(y, d, method="multinom", withTies = TRUE, Var = FALSE, level = 0.95, \dots)
}
\description{
compute the Polytomous Discrimination Index (PDI) value for a classifier, with optional variance estimation, using either a specified built-in classification method or user-supplied class probabilities.
}
\arguments{
  \item{y}{A multinomial response vector with two or more categories. It can be a factor or an integer-valued vector.}
  \item{d}{A set of candidate markers with one or more columns. It can be a data frame or a matrix. If \code{method = "prob"}, then \code{d} should be a matrix of predicted class probabilities.}
  \item{method}{The method used to construct the classifier based on the marker set in \code{d}. Available options are \code{"multinom"} (multinomial logistic regression; default; requires package \pkg{nnet}), \code{"tree"} (classification tree; requires package \pkg{rpart}), \code{"svm"} (support vector machine classification; requires package \pkg{e1071}), \code{"lda"} (linear discriminant analysis; requires package \pkg{MASS}), and \code{"prob"} (user-supplied predicted class probabilities).}
  \item{withTies}{A logical value indicating whether tied predicted probabilities may occur in the PDI calculation. This argument is only used when \code{Var = TRUE}. The default is \code{TRUE}.}
  \item{Var}{A logical value indicating whether to compute the variance and confidence interval of the PDI. The default is \code{FALSE}.}
  \item{level}{The confidence level for confidence interval construction. The default is \code{0.95}.}
  \item{\dots}{Additional arguments passed to the function corresponding to the chosen classification method.}
}
\details{
The function returns the PDI value for predictive markers based on a user-chosen machine learning method. Currently available methods include logistic regression (default), tree, lda, svm and user-computed risk values. This function is general since we can evaluate the accuracy for marker combinations resulted from complicated classification algorithms.


}
\value{
Returns an object of class \code{"mcca.pdi"}. The PDI value of the classifier obtained using a particular learning method on a set of marker(s).

An object of class \code{"mcca.pdi"} is a list containing at least the following components:
  \item{call}{the matched call.}
  \item{measure}{the value of the PDI.}
  \item{table}{the category-specific values of the PDI.}
If \code{Var = TRUE}, the list also contains the following:
  \item{se}{the calculated standard error of the overall PDI.}
  \item{ci}{the confidence interval for the overall PDI.}
  \item{level}{the confidence level used.}
and component \code{table} also contains the standard errors and confidence intervals of the category-specific PDI values.
}
\references{
Feng, Q., Jia, T., Liu, P., Van Calster, B., Li J. (2026). Efficient variance estimation for the polytomous discrimination index. Statistics in Medicine.

Dover, D.C., Savu, A., Engel, B. (2021). Polytomous discrimination index: Estimation and inference. Statistics in Medicine.

Li, J., Gao, M., D’Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.

Van Calster, B., Vergouwe, Y., Looman, C.W.N., Van Belle, V., Timmerman, D. and Steyerberg, E.W. Assessing the discriminative ability of risk models for more than two outcome categories. European Journal of Epidemiology 2012; 27: 761-770.

Li, J., Feng, Q., Fine, J.P., Pencina, M.J., Van Calster, B. (2018). Nonparametric estimation and inference for polytomous discrimination index. Statistical Methods in Medical Research. 27(10): 3092—3103.
}
\author{
Ming Gao: gaoming@umich.edu

Pan Liu: pan.liu@nus.edu.sg

Jialiang Li: stalj@nus.edu.sg
}
\note{
Users are advised to change the operating settings of various classifiers since it is well known that machine learning methods require extensive tuning. Currently only some common and intuitive options are set as default and they are by no means the optimal parameterization for a particular data analysis. Users can put machine learning methods' parameters after tuning. A more flexible evaluation is to consider "method=prob" in which case the input d should be a matrix of membership probabilities with k columns and each row of d should sum to one.
}
\seealso{
\code{\link{hum}}
}
\examples{
str(iris)
data <- iris[, 3]
label <- iris[, 5]
pdi(y = label, d = data,method = "multinom")

## Call:
## pdi(y = label, d = data, method = "multinom")

## Overall Polytomous Discrimination Index:
##  0.9881333

## Category-specific Polytomous Discrimination Index:
##    CATEGORIES VALUES
##  1     setosa 1.0000
##  2 versicolor 0.9822
##  3  virginica 0.9822

# PDI with variance estimation
pdi(y = label, d = data, method = "multinom", Var = TRUE)

## Call:
## pdi(y = label, d = data, method = "multinom", Var = TRUE)

## Overall Polytomous Discrimination Index:
##  0.9881333

## Standard Error:
##  0.006440373

## 95\% Confidence Interval:
##  [0.9755104, 1.0000000]

## Category-specific Polytomous Discrimination Index:
##  CATEGORIES VALUES     SE LOWER_CI UPPER_CI
##      setosa 1.0000 0.0000   1.0000   1.0000
##  versicolor 0.9822 0.0097   0.9633   1.0000
##   virginica 0.9822 0.0097   0.9633   1.0000

pdi(y = label, d = data,method = "tree")
pdi(y = label, d = data,method = "tree",control = rpart::rpart.control(minsplit = 200))

data <- data.matrix(iris[, 3])
label <- as.numeric(iris[, 5])
# multinomial
require(nnet)
# model
fit <- multinom(label ~ data, maxit = 1000, MaxNWts = 2000)
predict.probs <- predict(fit, type = "probs")
pp<- data.frame(predict.probs)
# extract the probablity assessment vector
head(pp)
pdi(y = label, d = pp, method = "prob")
}
\keyword{ PDI }
