% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rcode.r
\name{FOS_ad}
\alias{FOS_ad}
\title{The adaptive FOS method}
\usage{
FOS_ad(x, alpha = 0.0027, pn = 0.9)
}
\arguments{
\item{x}{Input data vector}

\item{alpha}{False alarm rate}

\item{pn}{Nominal coverage probability}
}
\value{
the lower and upper confidence bounds
\itemize{
 \item{ LCL:}{ The lower confidence limit}
 \item{ UCL:}{ The upper confidence limit}
 }
}
\description{
The adaptive FOS method
}
\examples{
FOS_ad(x6, alpha = 0.0027, pn = 0.3)
}
