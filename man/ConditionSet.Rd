\name{ConditionSet}
\alias{ConditionSet}
\alias{ConditionSet.add}
\title{ConditionSet}
\description{
An marker interface for classes that represents a \code{Condition} that can be satisfied.
}
\usage{
ConditionSet.add(self, condition)
}
\details{
  An marker interface for classes that represents a \code{Condition} that can be satisfied.\cr
  All methods that require an instance of \code{ConditionSet} object can use the \code{$} method to invoke the desired method.\cr
  For example:\cr
  \code{
    obj$add(SchedulingConditionFactory.createQueueJumpCondition())
  }

  \describe{
	\item{ConditionSet.add(self, condition)}{
	  Adds a condition.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}      \tab An instance of \code{ConditionSet}\cr
			\code{condition} \tab An instance of \code{Condition}
		  }
		}
	  }
	}
  }

}
\seealso{
\code{\link{ConditionSetFactory}},
\code{\link{SchedulingConditionFactory}}
}
\examples{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
