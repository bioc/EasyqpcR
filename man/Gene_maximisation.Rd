\name{Gene_maximisation}
\alias{Gene_maximisation}
\docType{data}
\title{
A data frame containing all the data from four different qPCR runs.
}
\description{
A data frame containing all the data from four different qPCR runs where samples
 and genes are spread accross runs.
}
\usage{data(Gene_maximisation)}
\format{
  A data frame with 192 observations on the following 4 variables.
  \describe{
    \item{\code{Samples}}{a factor with levels \code{} \code{IRC 1.1} 
 \code{IRC 1.2} \code{IRC 1.3} \code{IRC 1.4} \code{IRC 1.5} \code{IRC 1.6} 
 \code{IRC 1.7} \code{IRC 1.8} \code{IRC 2.1} \code{IRC 2.2} \code{IRC 2.3} 
 \code{IRC 2.4} \code{IRC 2.5} \code{IRC 2.6} \code{IRC 2.7} \code{IRC 2.8} 
 \code{IRC 3.1} \code{IRC 3.2} \code{IRC 3.3} \code{IRC 3.4} \code{IRC 3.5} 
 \code{IRC 3.6} \code{IRC 3.7} \code{IRC 3.8} \code{NTC 1} \code{NTC 2} 
 \code{NTC 3} \code{NTC 4} \code{NTC 5} \code{NTC 6} \code{NTC 7} \code{NTC 8} 
 \code{Sample 1} \code{Sample 10} \code{Sample 11} \code{Sample 12} 
 \code{Sample 13} \code{Sample 14} \code{Sample 15} \code{Sample 16} 
 \code{Sample 17} \code{Sample 18} \code{Sample 19} \code{Sample 2} 
 \code{Sample 20} \code{Sample 21} \code{Sample 22} \code{Sample 23} 
 \code{Sample 24} \code{Sample 25} \code{Sample 26} \code{Sample 27} 
 \code{Sample 28} \code{Sample 29} \code{Sample 3} \code{Sample 30} 
 \code{Sample 31} \code{Sample 32} \code{Sample 4} \code{Sample 5} 
 \code{Sample 6} \code{Sample 7} \code{Sample 8} \code{Sample 9}}
    \item{\code{RG1}}{a numeric vector}
    \item{\code{RG2}}{a numeric vector}
    \item{\code{TG}}{a numeric vector}
  }
}
\details{
This data frame is composed by 64 samples including IRCs (3 per gene), NTC (1 
 per gene). The control samples are Samples 1 to 16, and the test samples are
 Samples 17 to 32. There are 2 reference genes (RG1, RG2) and one target gene  
 (TG).
}
\source{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\references{
S. Le Pape, IRTOMIT-INSERM U1082.
}
\examples{
data(Gene_maximisation)
}
\keyword{datasets}
