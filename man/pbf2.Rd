\name{pbf2}
\alias{pbf2}
\title{prob being in response function}
\description{Calculate the probability of being in response function (PBRF) and the variances}
\usage{
  pbf2(y1,y2,d1,d2,times=y2[order(y2)])
}
\arguments{
  \item{y1}{ a numeric vector of event times denoting the minimum of event times \eqn{T_1}, \eqn{T_2} and censoring time \eqn{C}, where \eqn{T_2} 
    corresponds to the event time, \eqn{T_1} corresponds to the response time.}
  \item{y2}{ a numeric vector of event times denoting the minimum of event time \eqn{T_2} and censoring time \eqn{C}. Clearly, y2 is not smaller than y1.}
  \item{d1}{ a numeric vector of event indicators with 1 denoting the response is observed and 0 denoting otherwise.}
  \item{d2}{ a numeric vector of event indicators with 1 denoting the event is observed and 0 denoting otherwise.}
  \item{times}{ a numeric vector of timepoints at which we want to estimate the PBRF}
}
\details{There three methods to estimate PBRF: the subtraction, the division and the Semi-Markov methods are presented in Tsai et al. (2017).
There are two sub-methods for division and the Semi-Markov methods when the censoring distrbution is estimated in two different ways (looking at y1=min(\eqn{T_1},\eqn{T_2},\eqn{C}) and at y2=min(\eqn{T_2},\eqn{C})). So there are 5 methods in total reported. Method 1: division and based on y1; Method 2: division and based on y2; Method 3: Semi-Markov and based on y1; Method 4: Semi-Markov and based on y2; Method 5: subtraction. 
The methods based on y2 perform better than the corresponding ones based on y1.  
}
\value{
  \item{pbrf}{The estimates at each timepoints (row) and by methods 1-5 (column)}
  \item{vpbrf}{The variance estimates at each timepoints (row) and by methods 1-5 (column)}
}
\references{
Tsai W.Y., Luo X., Crowley J. (2017) The Probability of Being in Response Function and Its Applications. In: Matsui S., Crowley J. (eds) Frontiers of Biostatistical Methods and Applications in Clinical Oncology. Springer, Singapore. <doi: 10.1007/978-981-10-0126-0_10>.
}
\author{Xiaodong Luo}
\examples{
  n<-300
  rho<-0.5
  lambda10<-0.1;lambda20<-0.08;lambdac0<-0.09
  lam1<-rep(0,n);lam2<-rep(0,n);lamc<-rep(0,n)
  z<-rep(0,n)
  z[1:(n/2)]<-1
  lam1<-lambda10
  lam2<-lambda20
  lamc<-lambdac0
  tem<-matrix(0,ncol=3,nrow=n)
  
  y2y<-matrix(0,nrow=n,ncol=3)
  y2y[,1]<-rnorm(n);y2y[,3]<-rnorm(n)
  y2y[,2]<-rho*y2y[,1]+sqrt(1-rho^2)*y2y[,3]
  tem[,1]<--log(1-pnorm(y2y[,1]))/lam1
  tem[,2]<--log(1-pnorm(y2y[,2]))/lam2
  tem[,3]<--log(1-runif(n))/lamc
  
  y1<-apply(tem,1,min)
  y2<-apply(tem[,2:3],1,min)
  d1<-as.numeric(tem[,1]<=y1)
  d2<-as.numeric(tem[,2]<=y2)
  btemp<-pbf2(y1,y2,d1,d2,times=c(1,3,5))
  btemp
}
\keyword{composite endpoints}
\keyword{survival analysis}
\keyword{multiple events}





