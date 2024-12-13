% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Rcode.r
\name{FOS_3terms}
\alias{FOS_3terms}
\title{The 3-term FOS method}
\usage{
FOS_3terms(x, alpha = 0.0027, pn = 0.9)
}
\arguments{
\item{x}{Input data vector}

\item{alpha}{False alarm rate}

\item{pn}{Nominal coverage probability}
}
\value{
The lower and upper confidence bounds
\itemize{
 \item{ LCL:}{ The lower confidence limit}
 \item{ UCL:}{ The upper confidence limit}
 }
}
\description{
The 3-term FOS method
}
\examples{
FOS_3terms(x6, alpha = 0.0027, pn = 0.3)
}
