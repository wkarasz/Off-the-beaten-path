\name{Discriminator}
\alias{Discriminator}
\alias{Discriminator.setComparator}
\alias{BasicEngineDiscriminator.setComparator}
\alias{PropertyDiscriminator.setComparator}
\alias{PropertyScoringDiscriminator.setComparator}
\title{Discriminator}
\description{
\code{Discriminator} defines the interface for customer-defined job and task-level discriminators.
}
\usage{
  Discriminator.setComparator(self, name, op, value, nullCompare = TRUE)
  BasicEngineDiscriminator.setComparator(self, name, op, value, nullCompare = TRUE)
  PropertyDiscriminator.setComparator(self, name, op, value, nullCompare = TRUE)
  PropertyScoringDiscriminator.setComparator(self, name, op, value, nullCompare = TRUE)
}
\details{
  \code{Discriminator} defines the interface for customer-defined Service and Task-level discriminators.\cr
  Such a discriminator must define the allow method, which enables discrimination based on properties of the Engine. The properties that are given on any Engine are in \code{\link{EngineProperties}}. Any properties defined in the Engine tracker, Engine properties page, and \code{EngineSession}, are also passed in.
  \bold{Any customer-defined discriminator classes must be placed in the Manager's lib directory}

  \describe{
	\item{Discriminator.setComparator(self, name, op, value, nullCompare = TRUE)}{
	  Discriminates based on given properties and comparators.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}        \tab An instance of \code{Discriminator}\cr
			\code{name}        \tab the property name. See \code{\link{EngineProperties}}\cr
			\code{op}          \tab the comparison operation See \code{\link{Discriminator.EComparator}}\cr
			\code{value}       \tab the value against which comparisons will be made\cr
			\code{nullCompare} \tab result when the value is null (default TRUE)
		  }
		}
	  }
	}
  }
}
\seealso{
code{\link{EngineProperties}}, \code{\link{Discriminator.EComparator}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
