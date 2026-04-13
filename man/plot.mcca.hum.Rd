\name{plot.mcca.hum}
\alias{plot.mcca.hum}
\title{Plot 3D ROC surface}
\usage{
\method{plot}{mcca.hum}(x,labs=levels(x$y),coords=1:3,nticks=5,filename='fig.png',cex=0.7, \dots)
}
\description{
plot the 3D ROC surface for a three-category classifier using the 3-dimensional point coordinates, computed by obj which is a mcca.hum class.
}
\arguments{
  \item{x}{An mcca.hum class object, containing probability matrix and labels.}
  \item{labs}{The label names of three coordinates. Default is `levels(x$y)`.}
  \item{coords}{The coordinates markers. Default is `c(1,2,3)`, which means labs[1] is the x-axis (class 1), labs[2] is the z-axis (class 3) and labs[3] is the y-axis (class 2).}
  \item{nticks}{Suggested number of ticks.}
  \item{filename}{Filename to save snapshot.}
  \item{cex}{Size for text.}
  \item{\dots}{further arguments to `plot.default`.}
}
\details{
This function is to plot the 3D ROC surface according to the correct classification probabilities for the three categories, resulted from any statistical or machine learning methods. This function complements the HUM package which can only plot 3D ROC surface for a single diagnostic test or biomarker for three classes.

}
\value{
The function doesn't return any value.
}
\references{
Li, J., and Zhou, X. H. (2009). Nonparametric and semiparametric estimation of the three way receiver operating characteristic surface. Journal of Statistical Planning and Inference. 139: 4133—4142.

Li, J., Gao, M., D’Agostino, R. (2019). Evaluating Classification Accuracy for Modern Learning Approaches. Statistics in Medicine (Tutorials in Biostatistics). 38(13): 2477-2503.
}
\author{
Ming Gao: gaoming@umich.edu

Pan Liu: pan.liu@nus.edu.sg

Jialiang Li: stalj@nus.edu.sg
}

\examples{
data <- iris[, 1]
label <- iris[, 5]
a=hum(y = label, d = data,method = "multinom")
#plot(a,filename='fig.png')

}
\keyword{ Plot }
