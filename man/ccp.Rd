\name{ccp}
\alias{ccp}
\title{Calculate CCP Value}
\usage{
ccp(y, d, method="multinom", \dots)
}
\description{
compute the Correct Classification Percentage (CCP) value of two or three or four categories classifiers with an option to define the specific model or user-defined model.
}
\arguments{
  \item{y}{The multinomial response vector with two, three or four categories.}
  \item{d}{The set of candidate markers, including one or more columns. Can be a data frame or a matrix; if the method is "label", then d should be the label vector.}
  \item{method}{Specifies what method is used to construct the classifier based on the marker set in d. Available option includes the following methods:"multinom": Multinomial Logistic Regression which is the default method, requiring R package nnet;"tree": Classification Tree method, requiring R package rpart; "svm": Support Vector Machine (C-classification and radial basis as default), requiring R package e1071; "lda": Linear Discriminant Analysis, requiring R package lda; "label": d is a label vector resulted from any external classification algorithm obtained by the user, should be encoded from 1; "prob": d is a probability matrix resulted from any external classification algorithm obtained by the user.}
  \item{\dots}{Additional arguments in the chosen method's function.}
}
\details{
The function returns the CCP value for predictive markers based on a user-chosen machine learning method. Currently available methods include logistic regression (default), tree, lda, svm and user-computed risk values. This function is general since we can evaluate the accuracy for marker combinations resulted from complicated classification algorithms.

}
\value{
Returns an object of class "mcca.ccp". The CCP value of the classification using a particular learning method on a set of marker(s).

An object of class "mcca.ccp" is a list containing at least the following components:
  \item{call}{the matched call.}
  \item{measure}{the value of measure.}
  \item{table}{the category-specific value of measure.}
}
\references{
Li, J., Gao, M., D’Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.

Li, J., Jiang, B. and Fine, J. P. (2013). Multicategory reclassification statistics for assessing Improvements in diagnostic accuracy. Biostatistics. 14(2): 382-394.

Li, J., Jiang, B., and Fine, J. P. (2013). Letter to Editor: Response. Biostatistics. 14(4): 809-810.
}
\author{
Ming Gao: gaoming@umich.edu

Pan Liu: pan.liu@nus.edu.sg

Jialiang Li: stalj@nus.edu.sg
}
\note{
Users are advised to change the operating settings of various classifiers since it is well known that machine learning methods require extensive tuning. Currently only some common and intuitive options are set as default and they are by no means the optimal parameterization for a particular data analysis. Users can put machine learning methods' parameters after tuning. A more flexible evaluation is to consider "method=label" in which case the input d should be a label vector.
}
\seealso{
\code{\link{pdi}}
}
\examples{
str(iris)
data <- iris[, 1:4]
label <- iris[, 5]
ccp(y = label, d = data, method = "multinom",maxit = 1000,MaxNWts = 2000,trace=FALSE)

## Call:
## ccp(y = label, d = data, method = "multinom", maxit = 1000,
##     MaxNWts = 2000, trace = FALSE)

## Overall Correct Classification Probability:
##  0.9866667

## Category-specific Correct Classification Probability:
##   CATEGORIES VALUES PREVALENCE
## 1     setosa   1.00  0.3333333
## 2 versicolor   0.98  0.3333333
## 3  virginica   0.98  0.3333333

ccp(y = label, d = data, method = "multinom")
ccp(y = label, d = data, method = "svm")
ccp(y = label, d = data, method = "svm",kernel="sigmoid",cost=4,scale=TRUE,coef0=0.5)
ccp(y = label, d = data, method = "tree")

p = as.numeric(label)
ccp(y = label, d = p, method = "label")


table(mtcars$carb)
for (i in (1:length(mtcars$carb))) {
  if (mtcars$carb[i] == 3 | mtcars$carb[i] == 6 | mtcars$carb[i] == 8) {
    mtcars$carb_new[i] = 9
  }else{
    mtcars$carb_new[i] = mtcars$carb[i]
  }
}
data <- data.matrix(mtcars[, c(1)])
mtcars$carb_new <- factor(mtcars$carb_new)
label <- mtcars$carb_new
str(mtcars)
ccp(y = as.numeric(label), d = data, method = "svm",kernel="radial",cost=1,scale=TRUE)

## Call:
## ccp(y = as.numeric(label), d = data, method = "svm", kernel = "radial", cost = 1, scale = TRUE)

## Overall Correct Classification Probability:
##  0.4375

## Category-specific Correct Classification Probability:
##   CATEGORIES    VALUES PREVALENCE
## 1          1 0.5714286    0.21875
## 2          2 0.2000000    0.31250
## 3          3 0.8000000    0.31250
## 4          4 0.0000000    0.15625

}
\keyword{ CCP }
