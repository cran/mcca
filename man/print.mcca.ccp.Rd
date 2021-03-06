\name{print.mcca.ccp}
\alias{print.mcca.ccp}
\title{Print Method for mcca ccp class}
\usage{
\method{print}{mcca.ccp}(x, ...)
}
\description{
An S3 method for the print generic. It is designed for a quick look at CCP values.
}
\arguments{
  \item{x}{object of class `mcca.ccp`.}
  \item{\dots}{further arguments to `print.default`.}
}

\value{
An object of class "mcca.ccp" is a list containing at least the following components:
  \item{call}{the matched call.}
  \item{measure}{the value of measure.}
  \item{table}{the category-specific value of measure.}
}
\references{
Li, J., Jiang, B. and Fine, J. P. (2013). Multicategory reclassification statistics for assessing Improvements in diagnostic accuracy. Biostatistics. 14(2): 382-394.

Li, J., Jiang, B., and Fine, J. P. (2013). Letter to Editor: Response. Biostatistics. 14(4): 809-810.
}
\author{
Ming Gao: gaoming@umich.edu

Jialiang Li: stalj@nus.edu.sg
}

\examples{
data = iris[, 1:4]
label = iris[, 5]
ccp_object=ccp(y = label, d = data, method = "multinom",maxit = 1000,MaxNWts = 2000,trace=FALSE)
print(ccp_object)

}
\keyword{ CCP }
