\name{totData}
\alias{totData}

\title{
Aggregation of qPCR biological replicates and data transformation
}
\description{
This function aggregates qPCR biological replicates and calculates the main
 parameters : mean (arithmetic or geometric), the standard deviation and the
 standard error from your biological replicates of your experience. This
 function has an algorithm published by Willems et al. (2008) which performs a
 standardization procedure that can be applied to data sets that display high
 variation between biological replicates. This enables proper statistical
 analysis and drawing relevant conclusions. The procedure is not new, it has
 been used in microarray data analysis and is based on a series of sequential
 corrections, including log transformation, mean centering, and autoscaling.
}
\usage{
totData(data, r, base, nSpl, nCTL, na.rm = na.rm)
}

\arguments{
  \item{data}{
data.frame containing row datas (genes in columns, samples in rows, Cq values).
}
  \item{r}{
numeric, number of qPCR replicates.
}
  \item{base}{
numeric, the logarithmic base (2 or 10).
}
  \item{nSpl}{
numeric, the number of samples.
}
  \item{nCTL}{
numeric, number of samples forming your control group.
}
\item{na.rm}{
logical, indicating whether NA values should be stripped before the computation
 proceeds.
}
}
\details{
The standardization procedure used in this function (if TRUE for the
 transformation argument) is based on the article of Willems et al. (2008).
 This function perform successively thEerror operations : log-transformation of
 your raw data, mean of your log-transformed data for each biological replicate,
 mean centering for each biological replicate, standard deviation of each
 mean-centered biological replicate, autoscaling of your data, i.e., your
 mean-centered data for each biological replicate will be divided by the
 standard deviation of the mean-centered biological replicate and then
 multiplicated by the mean of the standard deviation of all the biological
 replicates.

For more information for the way to use this function, please see the vignette.
}
\value{

 \item{Standardized data according to Willems procedure (2008)}{The 
 standardized data according to the Willems procedure (2008). These data must 
 be used for statistical analysis.}
 \item{Mean}{The mean of the biological replicates for each sample.}
 \item{Lower Confidence Interval 95\%}{The Lower Confidence Interval 95\% of 
 the biological replicates for each sample.}
 \item{Upper Confidence Interval 95\%}{The Upper Confidence Interval 95\% of 
 the biological replicates for each sample.}


}
\references{
Erik Willems Luc Leyns, Jo Vandesompele. Standardization of real-time PCR gene
 expression data from independent biological replicates. Analytical Biochemistry
 379 (2008) 127-129 (doi:10.1016/j.ab.2008.04.036).
 <url:http://www.sciencedirect.com/science/article/pii/S0003269708002649>
}
\author{
Sylvain Le pape <sylvain.le.pape@univ-poitiers.fr>
}

\examples{
data(qPCR_run1,qPCR_run2,qPCR_run3)

nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
        Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)


## Isolate the calibrator NRQ values of the first biological replicate

a <- nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]] 

## Isolate the calibrator NRQ values of the first biological replicate

b <- nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]]
 
## Isolate the calibrator NRQ values of the first biological replicate

c <- nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=c(1, 1, 1, 1), CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[[3]]

## Regrouping the calibrator NRQ values of all the biological replicates

d <- rbind(a, b, c) 

## Calibration factor calculation

e <- calData(d) 


## Attenuation of inter-run variation thanks to the calibration factor for the
## first biological replicate

nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)

## Attenuation of inter-run variation thanks to the calibration factor for the
## second biological replicate

nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE) 

## Attenuation of inter-run variation thanks to the calibration factor for the
## third biological replicate

nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE) 

## Isolate the NRQs scaled to control of the first biological replicate

a1 <- nrmData(data = qPCR_run1 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1] 

## Isolate the NRQs scaled to control of the second biological replicate 

b1 <- nrmData(data = qPCR_run2 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1] 

## Isolate the NRQs scaled to control of the third biological replicate

c1 <- nrmData(data = qPCR_run3 , r=3, E=c(2, 2, 2, 2), 
	      Eerror=c(0.02, 0.02, 0.02, 0.02), nSpl=5, 
	      nbRef=2, Refposcol=1:2, nCTL=2, 
	      CF=e, CalPos=5, trace=TRUE, geo=TRUE, na.rm=TRUE)[1]

## Data frame transformation

a2 <- as.data.frame(a1)
b2 <- as.data.frame(b1)
c2 <- as.data.frame(c1)

## Aggregation of the three biological replicates

d2 <- rbind(a2, b2, c2)

totData(data=d2, r=3, base=2, nSpl=5, nCTL = 2,
	na.rm=TRUE)

}

\keyword{Biological replicates}
\keyword{Standardization procedure}
