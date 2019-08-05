\name{Discriminator.EComparator}
\alias{EComparator}
\alias{Discriminator.EComparator}
\alias{Discriminator.EComparator.CONTAINS}
\alias{Discriminator.EComparator.EQUAL}
\alias{Discriminator.EComparator.GREATER_THAN}
\alias{Discriminator.EComparator.GREATER_THAN_EQ}
\alias{Discriminator.EComparator.LESS_THAN}
\alias{Discriminator.EComparator.LESS_THAN_EQ}
\alias{Discriminator.EComparator.MATCHES}
\alias{Discriminator.EComparator.NOT_CONTAINS}
\alias{Discriminator.EComparator.NOT_EQUAL}
\alias{Discriminator.EComparator.NOT_MATCHES}
\docType{data}
\title{
Comparators used with Discriminators
}
\description{
These are the \code{op} used with \code{\link{Discriminator.setComparator}}
}
\value{
\item{Discriminator.EComparator.CONTAINS }{The value is a comma-delimited list of values, of which the Engine value must match one.\cr The wildcard character (*) is allowed.}
\item{Discriminator.EComparator.EQUAL }{The parsed numerical values are equivalent.}
\item{Discriminator.EComparator.GREATER_THAN }{The parsed numerical Engine value must be greater than the given parsed value.}
\item{Discriminator.EComparator.GREATER_THAN_EQ }{The parsed numerical Engine value must be greater than or equal to the given parsed value.}
\item{Discriminator.EComparator.LESS_THAN }{The parsed numerical Engine value must be less than the given parsed value.}
\item{Discriminator.EComparator.LESS_THAN_EQ }{The parsed numerical Engine value must be greater than or equal to the given parsed value.}
\item{Discriminator.EComparator.MATCHES }{The value is a Java regular expression.\cr The Engine value must match the expression.}
\item{Discriminator.EComparator.NOT_CONTAINS }{The value is a comma-delimited list of values, none of which The Engine value matches.\cr The wildcard character (*) is allowed.}
\item{Discriminator.EComparator.NOT_EQUAL }{The parsed numerical values are not equivalent.}
\item{Discriminator.EComparator.NOT_MATCHES }{The value is a Java regular expression.\cr The Engine value must not match the expression.}
}
\seealso{
\code{\link{EngineProperties}}, \code{\link{Discriminator.setComparator}}
}
\keyword{datasets}
