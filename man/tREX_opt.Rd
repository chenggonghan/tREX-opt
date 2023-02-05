\name{tREXopt}
\alias{tREXopt}
\title{

}
\description{
tREX (Park, J. and Lin, S. (2016)) is a model-based method to construct 3D structure from Hi-C matrix which takes dependency of Hi-C data into consideration and addresses the excessive zero issues. However, its computational cost is very expensive.Here we propose tREX-opt, a 3D structure recapitulation method to combine the advantages from both modeling-based and optimization-based approaches. Specifically, we form an objective function based on the tREX model that has greater fidelity to real data. 
}
\usage{
trex_opt(P.int=NULL,N.int=3,X,max.iter=50,tol=1e-4,method.type='tRex',parallel=FALSE,cpu=4)
}
\arguments{
  \item{P.int}{
  A N*3 initial coordinate matrix. If not specified, a randomly generated coordinate matrix will be used. The default value is NULL.
}
\item{N.int}{
	The number of the columns of the initial P.int. The default value is 3.
}
\item{X}{
   The Hi-C matrix. Only the lower triangular values are needed.
}
\item{max.iter}{
  The maximum number of iterations..
}
\item{tol}{
	 The relative difference in log likelihood between two adjacent iterations
}
  \item{method.type}{
	A string value. The method options. method.type should be one of the following: tRex, tPAM, PRAM, PAM.
}
  \item{parallel}{
	 A logical value enabling the use of parallel function. The default is FALSE. 
}
  \item{cpu}{
 The number of CPUs of parallelization when parallelization is turned on.
}

\value{
It returns a list containing the following dimensions:
  \item{struct}{The output 3D coordinate matrix}
  \item{neg_loglik}{negative loglikelihood function value}
  \item{parameter}{Model parameters}
}

\references{
Park, J. and Lin, S. (2016). Impact of data resolution on three-dimensional structure inference methods. BMC bioinformatics, 17(1), 70.
}

\examples{
library(trexopt)
load('sampledata.rda')
X<-matrix(0,nrow=61,ncol = 61)
X[lower.tri(X)] <- sample_data
result<-trex_opt(X,max.iter=100)
}

