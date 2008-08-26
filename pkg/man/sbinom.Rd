\name{sbinom}
\alias{sbinom}
\alias{spois}
\alias{snbinom}
\alias{sbeta}
\title{Abstract definitions of distributions}
\description{
  Functions returning values for summary statistics
  (mean, median, etc.) of distributions
}
\usage{
sbinom(size,prob)
snbinom(size, prob, mu)
spois(lambda)
sbeta(shape1,shape2)
}
\arguments{
  \item{prob}{probability as defined for \code{dbinom} or \code{dnbinom}}
  \item{size}{size/variance parameter as defined for \code{dnbinom}}
  \item{mu}{mean parameter as defined for \code{dnbinom}}
 \item{shape1}{shape parameter for \code{dbeta}}
 \item{shape2}{shape parameter for \code{dbeta}}
  \item{lambda}{rate parameter as defined for \code{dpois}}
}
\value{
  \item{title}{name of the distribution}
  \item{[parameters]}{input parameters for the distribution}
  \item{mean}{theoretical mean of the distribution}
  \item{median}{theoretical median of the distribution}
  \item{mode}{theoretical mode of the distribution}
  \item{variance}{theoretical variance of the distribution}
  \item{sd}{theoretical standard deviation of the distribution}
}
\author{Ben Bolker}
\seealso{\code{\link{dbinom}},\code{\link{dpois}}, \code{\link{dnbinom}}}
\examples{
  sbinom(prob=0.2,size=10)
  snbinom(mu=2,size=1.2)
}
\note{these definitions are tentative, subject to change
as I figure this out better.  Perhaps construct functions
that return functions? Strip down results? Do more
automatically?}

\keyword{misc}
