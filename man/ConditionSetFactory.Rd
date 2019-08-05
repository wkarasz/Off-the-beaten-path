\name{ConditionSetFactory}
\alias{ConditionSetFactory}
\alias{ConditionSetFactory.createANDSet}
\alias{ConditionSetFactory.createORSet}
\alias{ConditionSetFactory.createServiceSet}
\title{ConditionSetFactory}
\description{
  A Factory for creating ConditionSets.
}
\usage{
ConditionSetFactory.createANDSet()
ConditionSetFactory.createORSet()
ConditionSetFactory.createServiceSet()
}
%- maybe also 'usage' for other objects documented here.
\details{
  A Factory for creating ConditionSets.

  \describe{
	\item{ConditionSetFactory.createANDSet()}{
	  Creates a set for which all Conditions must be satisfied for the set to be satisfied.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			the set
		  }
		}
	  }
	}
  }

  \describe{
	\item{ConditionSetFactory.createORSet()}{
	  Creates a set for which any \code{Condition} can be satisfied for the set to be satisfied.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			the set
		  }
		}
	  }
	}
  }

  \describe{
	\item{ConditionSetFactory.createServiceSet()}{
	  Creates a base set for adding different \code{Condition} types to an Invocation or Session.\cr
	  You can only add one of each type. Currently there are four types: Discriminator/Dependency, Affinity, QueueJump, and Descriptor.
	  So you can add one Discriminator or Dependency or set of such, one Affinity or set of such, one QueueJump, and one Descriptor (sets of Descriptors are not allowed).
	  The set can only be added to \code{ServiceFactory.createService(...)}, \code{Service.execute(...)}, and \code{Service.submit(...)}
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			the set
		  }
		}
	  }
	}
  }

}
\seealso{
\code{\link{ServiceFactory.createService}}, \code{\link{Service.execute}}, \code{\link{Service.submit}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
