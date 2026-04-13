\name{print.mcca.pdi.var}
\alias{print.mcca.pdi.var}
\title{Print Method for mcca pdi.var class}
\usage{
\method{print}{mcca.pdi.var}(x, ...)
}
\description{
An S3 method for the print generic. It is designed for a quick look at PDI values with variance estimation.
}
\arguments{
  \item{x}{object of class `mcca.pdi.var`.}
  \item{\dots}{further arguments to `print.default`.}
}

\value{
An object of class "mcca.pdi.var" is a list containing at least the following components:
  \item{call}{the matched call.}
  \item{measure}{the overall PDI value.}
  \item{se}{the calculated standard error of the overall PDI.}
  \item{ci}{the confidence interval for the overall PDI.}
  \item{level}{the confidence level used.}
  \item{table}{a data frame with category-specific PDI values, standard errors, and confidence intervals.}
}
\references{
Feng, Q., Jia, T., Liu, P., Van Calster, B., Li J. (2026). Efficient variance estimation for the polytomous discrimination index. Statistics in Medicine.

Dover, D.C., Savu, A., Engel, B. (2021). Polytomous discrimination index: Estimation and inference. Statistics in Medicine.

Li, J., Gao, M., D'Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.
}
\author{
Ming Gao: gaoming@umich.edu

Pan Liu: pan.liu@nus.edu.sg

Jialiang Li: stalj@nus.edu.sg
}

\examples{
data = iris[, 1:4]
label = iris[, 5]
pdi_var_object = pdi(y = label, d = data, method = "multinom", Var = TRUE)
print(pdi_var_object)

}
\keyword{PDI}
