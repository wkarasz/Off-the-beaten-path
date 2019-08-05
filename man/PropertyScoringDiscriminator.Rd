\name{PropertyScoringDiscriminator}
\alias{PropertyScoringDiscriminator}
\alias{PropertyScoringDiscriminator.new}
\alias{PropertyScoringDiscriminator.getMaxTasks}
\alias{PropertyScoringDiscriminator.setMaxTasks}
\title{PropertyScoringDiscriminator}
\description{
  Deprecated: Use \code{\link{SchedulingConditionFactory}}.
}
\usage{
  PropertyScoringDiscriminator()
  PropertyScoringDiscriminator.getMaxTasks(self)
  PropertyScoringDiscriminator.setMaxTasks(self, max)
}
\details{
  Deprecated: Use \code{\link{SchedulingConditionFactory}}.\cr
  Calculates a score which is the number of true property comparisons, based on a given comparison Properties set. All comparison properties must have a name, a comparator, and a value. See Discriminator for defined comparators and values.\cr
  All methods that require an instance of \code{ConditionSet} object can use the \code{$} method to invoke the desired method.\cr
  For example:\cr
  \code{
    max <- obj$getMaxTasks()
  }

  \describe{
	\item{PropertyScoringDiscriminator()}{
	  Creates an instance of \code{PropertyScoringDiscriminator} object.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			an instance.
		  }
		}
	  }
	}
  }

  \describe{
	\item{PropertyScoringDiscriminator.getMaxTasks(self)}{
	  The maximum number of Tasks to score.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}       \tab An instance of \code{PropertyScoringDiscriminator}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the number. Default 20
		  }
		}
	  }
	}
  }

  \describe{
	\item{PropertyScoringDiscriminator.setMaxTasks(self, max)}{
	  The maximum number of Tasks to score.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}       \tab An instance of \code{PropertyScoringDiscriminator}\cr
			\code{max}        \tab the number. Default 20
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the number. Default 20
		  }
		}
	  }
	}
  }

}
\seealso{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
