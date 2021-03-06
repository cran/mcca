\name{idi}
\alias{idi}
\title{Calculate IDI Value}
\usage{
idi(y, m1, m2, method="multinom", \dots)
}
\description{
compute the integrated discrimination improvement (IDI) value of two or three or four categories classifiers with an option to define the specific model or user-defined model.
}
\arguments{
  \item{y}{The multinomial response vector with two, three or four categories. It can be factor or integer-valued.}
  \item{m1}{The set of marker(s) included in the baseline model, can be a data frame or a matrix; if the method is "prob", then m1 should be the prediction probablity matrix of the baseline model.}
  \item{m2}{The set of additional marker(s) included in the improved model, can be a data frame or a matrix; if the method is "prob", then m2 should be the prediction probablity matrix of the improved model.}
  \item{method}{Specifies what method is used to construct the classifier based on the marker set in m1 & m2. Available option includes the following methods:"multinom": Multinomial Logistic Regression which is the default method, requiring R package nnet;"tree": Classification Tree method, requiring R package rpart;"svm": Support Vector Machine (C-classification and radial basis as default), requiring R package e1071;"lda": Linear Discriminant Analysis, requiring R package lda;"prob": m1 & m2 are risk matrices resulted from any external classification algorithm obtained by the user.}
  \item{\dots}{Additional arguments in the chosen method's function.}
}
\details{
The function returns the IDI value for predictive markers based on a user-chosen machine learning method. Currently available methods include logistic regression (default), tree, lda, svm and user-computed risk values. This function is general since we can evaluate the accuracy for marker combinations resulted from complicated classification algorithms.

}
\value{
The IDI value of the classification using a particular learning method on a set of marker(s).
}
\references{
Li, J., Gao, M., D’Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.

Li, J., Jiang, B. and Fine, J. P. (2013). Multicategory reclassification statistics for assessing Improvements in diagnostic accuracy. Biostatistics. 14(2): 382—394.

Li, J., Jiang, B., and Fine, J. P. (2013). Letter to Editor: Response. Biostatistics. 14(4): 809-810.
}
\author{
Ming Gao: gaoming@umich.edu

Jialiang Li: stalj@nus.edu.sg
}
\note{
Users are advised to change the operating settings of various classifiers since it is well known that machine learning methods require extensive tuning. Currently only some common and intuitive options are set as default and they are by no means the optimal parameterization for a particular data analysis. Users can put machine learning methods' parameters after tuning. A more flexible evaluation is to consider "method=prob" in which case the input m1 & m2 should be a matrix of membership probabilities with k columns and each row of m1 & m2 should sum to one.
}
\seealso{
\code{\link{nri}}
}
\examples{
table(mtcars$carb)
for (i in (1:length(mtcars$carb))) {
  if (mtcars$carb[i] == 3 | mtcars$carb[i] == 6 | mtcars$carb[i] == 8) {
    mtcars$carb_new[i] = 9
  }else{
    mtcars$carb_new[i] = mtcars$carb[i]
  }
}
data <- data.matrix(mtcars[, c(1,5)])
mtcars$carb_new <- factor(mtcars$carb_new)
label <- mtcars$carb_new
str(mtcars)

idi(y = label, m1 = data[, 1], m2 = data[, 2], "tree")
## [1] 0.09979413
idi(y = label, m1 = data[, 1], m2 = data[, 2], "tree",control=rpart::rpart.control(minsplit=4))
## [1] 0.2216707

}
\keyword{ IDI }
