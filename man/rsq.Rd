\name{rsq}
\alias{rsq}
\title{Calculate RSQ Value}
\usage{
rsq(y, d, method="multinom", \dots)
}
\description{
compute the R-squared (RSQ) value of two or three or four categories classifiers with an option to define the specific model or user-defined model.
}
\arguments{
  \item{y}{The multinomial response vector with two, three or four categories. It can be factor or integer-valued.}
  \item{d}{The set of candidate markers, including one or more columns. Can be a data frame or a matrix; if the method is "prob", then d should be the probablity matrix.}
  \item{method}{Specifies what method is used to construct the classifier based on the marker set in d. Available option includes the following methods:"multinom": Multinomial Logistic Regression which is the default method, requiring R package nnet;"tree": Classification Tree method, requiring R package rpart;"svm": Support Vector Machine (C-classification and radial basis as default), requiring R package e1071;"lda": Linear Discriminant Analysis, requiring R package lda;"prob": d is a risk matrix resulted from any external classification algorithm obtained by the user.}
  \item{\dots}{Additional arguments in the chosen method's function.}
}
\details{
The function returns the RSQ value for predictive markers based on a user-chosen machine learning method. Currently available methods include logistic regression (default), tree, lda, svm and user-computed risk values. This function is general since we can evaluate the accuracy for marker combinations resulted from complicated classification algorithms.

}
\value{
Returns an object of class "mcca.rsq". The RSQ value of the classification using a particular learning method on a set of marker(s).

An object of class "mcca.rsq" is a list containing at least the following components:
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
Users are advised to change the operating settings of various classifiers since it is well known that machine learning methods require extensive tuning. Currently only some common and intuitive options are set as default and they are by no means the optimal parameterization for a particular data analysis. Users can put machine learning methods' parameters after tuning. A more flexible evaluation is to consider "method=prob" in which case the input d should be a matrix of membership probabilities with k columns and each row of d should sum to one.
}
\seealso{
\code{\link{ccp}}
}
\examples{
str(iris)
data <- iris[, 1:4]
label <- iris[, 5]
rsq(y = label, d = data, method="multinom")

## Call:
## rsq(y = label, d = data, method = "multinom")

## Overall R-squared value:
##  0.9637932

## Category-specific R-squared value:
##   CATEGORIES    VALUES
## 1     setosa 0.9999824
## 2 versicolor 0.9456770
## 3  virginica 0.9457203

rsq(y = label, d = data, method = "tree")


data <- data.matrix(iris[, 1:4])
label <- as.numeric(iris[, 5])
# multinomial
require(nnet)
# model
fit <- multinom(label ~ data, maxit = 1000, MaxNWts = 2000)
predict.probs <- predict(fit, type = "probs")
pp<- data.frame(predict.probs)
# extract the probablity assessment vector
head(pp)
rsq(y = label, d = pp, method = "prob")


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

rsq(y = label, d = data, method="tree")
rsq(y = label, d = data, method="lda")
rsq(y = label, d = data, method="lda",prior = c(100,1,1,1)/103)

}
\keyword{ RSQ }
