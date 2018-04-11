\name{pETM}
\alias{pETM}
\title{Penalized Exponential Tilt Model}
\description{
    Fit a penalized exponential tilt model (ETM) to identify differentially
    methylated loci between cases and controls. ETM is able to detect any
    differences in means only, in variances only or in both means and 
    variances. \cr
    \cr
    A penalized exponential tilt model using combined lasso and Laplacian penalties
    is applied to high-dimensional DNA methylation data with case-control 
    association studies. When CpG sites are correlated with each other within
    the same gene or the same genetic region, Laplacian matrix can be imposed
    into the penalty function to encourage grouping effects among linked CpG
    sites. The selection probability of an individual CpG site is computed based
    on a finite number of resamplings.
    }
\usage{pETM(x,y,cx=NULL,alpha=0.1,maxit=100000,thre=1e-6,group=NULL,lambda=NULL,
      type=c("ring","fcon"),etm=c("none","normal","beta"),psub=0.5,nlam=10,
      kb=10,K=100)
    }
\arguments{
    \item{x}{Observed DNA methylation beta values consisting of \eqn{n} samples 
        and \eqn{p} CpG sites. It should be (\eqn{n} x \eqn{p}) design matrix
        without an intercept.}
    \item{y}{The phenotype outcome coded as 1 for cases and 0 for the 
        controls.}
    \item{cx}{The covariates such as age and gender. It should be 
        (\eqn{n} x \eqn{m}) matrix, where \eqn{m} is the number of the 
         covariates.}
    \item{alpha}{The penalty mixing parameter with \eqn{0\le\alpha\le 1} and
        default is 0.1. See details.}
    \item{maxit}{Maximum number of passes over the data for all regularization
        values, and default is 10^5. For fast computation, use a smaller value
        than the default value.}
    \item{thre}{Convergence threshold for coordinate descent algorithm.
        The default value is \code{1E-6}. For fast computation, use a larger
        value than the default value.}
    \item{group}{The integer vector describing the size of genes or genetic 
        regions. The length of \code{group} should be equivalent to the total 
        number of genes or genetic regions, and the sum of \code{group} should
        be the same as the total number of CpG sites. If no group information
        is available, i.e., not specified, the \code{pETM} performs an 
        elastic-net regularization procedure of a logistic regression.
        See details.}
    \item{lambda}{A sequence of regularization tuning parameter can be 
        specified. Typical usage is to have the program compute its own
        \code{lambda} sequence based on \code{nlam} and \code{kb}}.     
    \item{type}{A type of network within each group when \code{group} is 
        specified. "\code{ring}" and "\code{fcon}" represent a ring and fully
        connected network, respectively. Default is "\code{ring}". See details.}     
    \item{etm}{A type of an exponential tilt model. \code{none} does not perform
        an exponential tilt model, instead an ordinary penalized logistic 
        regression model is applied. \code{normal} performs a penalized 
        exponential tilt model based on a Gaussian distribution, and \code{beta}
        performs a penalized exponential tilt model based on a Beta 
        distribution. See details.}
    \item{psub}{The proportion of subsamples used for resamplings, and 
        \code{psub}\eqn{\in[0.5,1)}. The default is 0.5.}
    \item{nlam}{The number of \code{lambda} values used for resamplings,
        and default is 10. For fast computation, use a smaller value than the 
        default value.}
    \item{kb}{The number of burn-out replications before resamplings to properly
        adjust a sequence of \code{lambda} values and default is 10.}
    \item{K}{The number of resamplings, and default is 100.}    
        }
\details{
       The exponential tilt model based on a logistic regression is defined as
       \deqn{\log\frac{p(x_i)}{1-p(x_i)} = 
       \beta_0+h_1(x_i)^{T}\beta_1+h_2(x_i)^{T}\beta_2,}
       where \eqn{h_1(\cdot)} and \eqn{h_2(\cdot)} are pre-specified functions.
       For example \eqn{h_1(x)=x} and \eqn{h_2(x)=x^2} if \code{etm} is 
       \code{normal} and \eqn{h_1(x)=-\log(x)} and \eqn{h_2(x)=-\log(1-x)} if
       \code{etm} is \code{beta}. \cr
       \cr
       The penalty function of \code{pETM} is defined as
       \deqn{\alpha||\beta||_1+(1-\alpha)(\beta^{T}L\beta)/2,}
       where \eqn{L} is a Laplacian matrix describing a group structure of
       CpG sites. This penalty is equivalent to the Lasso penalty if
       \code{alpha=1}. When \code{group} is not defined, \eqn{L} is replaced by
       an identity matrix. In this case, \code{pETM} performs an elastic-net
       regularization procedure since the second term of the penalty simply
       reduces to the squared \eqn{l_2} norm of \eqn{\beta}. \cr
       \cr
       If group sizes of CpG sites are listed in \code{group}, it is assumed
       that CpG sites within the same genes are linked with each other like
       a ring or a fully connected network. In this case, the Laplacian matrix
       forms a block-wise diagonal matrix. The ring network assumes only 
       adjacent CpG sites within the same genes are linked with each other,
       while every CpG sites within the same genes are linked with each other
       for fully connected network. For a big gene, ring network is recommended
       for computational speed-up. \cr
       \cr
       The selection result is summarized as the selection probability of 
       individual CpG sites. The \code{psub} portions of \code{n} samples are
       randomly selected without replacement \code{K} times. For each 
       subsample of (\code{x},\code{cx},\code{y}), \code{pETM} is applied to
       find non-zero coefficients of CpG sites along with \code{nlam} lambda
       values. The selection probability of each CpG site is then computed
       based on the maximum proportion of non-zero regression coefficients
       among \code{K} replications.
       }
\value{
    \item{selprob}{The selection probabilities of \eqn{p} CpG sites}
    \item{topsp}{The selection probability of each CpG site is listed in
        descending order along with the name of CpG sites.}
    \item{lambda}{The actual sequence of \code{lambda} values used}
    \item{valid.K}{The actual number of resamplings used}
     }
\references{
    H. Sun and S. Wang (2012)
    \emph{Penalized Logistic Regression for High-dimensional DNA Methylation
    Data with Case-Control Studies}, Bioinformatics 28(10), 1368--1375 \cr
    \cr
    H. Sun and S. Wang (2013)
    \emph{Network-based Regularization for Matched Case-Control
    Analysis of High-dimensional DNA Methylation Data}, Statistics in Medicine
    32(12), 2127--2139 \cr
    \cr
    H. Sun and S. Wang (2016)
    \emph{Penalized Exponential Tilt Model for Analysis of High-dimensional
    DNA Methylation Data}, Manuscript \cr
    }
\author{Hokeun Sun <hsun@pusan.ac.kr>}
\examples{
    n <- 100
    p <- 500
    x <- matrix(rnorm(n*p), n, p)
    y <- rep(0:1, c(50,50))
     
    # a total of 200 genes each of which consists of 1, 2, or 5 CpG sites 
    gr <- rep(c(1,2,5), c(50,100,50))
    
    # ordinary penalized logistic regression   
    g1 <- pETM(x, y, group=gr, K=10) 
    
    # penalized exponential tilt model based on Gaussian distribution 
    g2 <- pETM(x, y, group=gr, etm = "normal", K=10) 
    
    # penalized exponential tilt model based on Beta distribution
    x2 <- matrix(runif(n*p), n, p) 
    g3 <- pETM(x2, y, group=gr, etm = "beta", K=10) 
}    



