\name{print.mcca.hum}
\alias{print.mcca.hum}
\title{Print Method for mcca hum class}
\usage{
\method{print}{mcca.hum}(x, ...)
}
\description{
An S3 method for the print generic. It is designed for a quick look at hum values.
}
\arguments{
  \item{x}{object of class `mcca.hum`.}
  \item{\dots}{further arguments to `print.default`.}
}

\value{
An object of class "mcca.hum" is a list containing at least the following components:
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

Pan Liu: pan.liu@nus.edu.sg

Jialiang Li: stalj@nus.edu.sg
}

\examples{
data = iris[, 1:4]
label = iris[, 5]
hum_object=hum(y = label, d = data, method = "multinom",maxit = 1000,MaxNWts = 2000,trace=FALSE)
print(hum_object)

}
\keyword{ hum }
