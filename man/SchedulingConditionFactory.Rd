\name{SchedulingConditionFactory}
\alias{SchedulingConditionFactory}
\alias{SchedulingConditionFactory.createCustomAffinity}
\alias{SchedulingConditionFactory.createCustomDiscriminator}
\alias{SchedulingConditionFactory.createPropertyAffinity}
\alias{SchedulingConditionFactory.createPropertyDiscriminator}
\alias{SchedulingConditionFactory.createQueueJumpCondition}
\title{SchedulingConditionFactory}
\description{
Factory for creating Conditions that affect how Service Sessions and Invocation are scheduled to to Engines.
}
\usage{
SchedulingConditionFactory.createCustomAffinity(className, props)
SchedulingConditionFactory.createCustomDiscriminator(className, props)
SchedulingConditionFactory.createPropertyAffinity(name, comparision, value, score)
SchedulingConditionFactory.createPropertyDiscriminator(name, comparision, value, nullCompare)
SchedulingConditionFactory.createQueueJumpCondition()
}

\details{
These Conditions perform both affinity and discrimination for \bold{Service} Session and Invocation requests, based on Engine properties. They can be Session- and/or Invocation-based.\cr
\bold{Session} \code{\link{Condition}} When a schedule event is performed, all free Engines are first filtered by any \bold{Discriminator} Conditions set on the session. if \emph{Consider Task Discriminators For Affinity} is set, it will also look up to the \emph{Task Discrimination Depth} and any Engines does not satisfy any tasks will also be filtered out. Then each is assigned an Affinity score to the Session based on the intrinsic affinity (homed, state, gridlibs) plus any Affinity Condition set on the session. The list is sorted, and Sessions are assigned their allocated amount of Engines in the order of highest to lowest affinity score. If at any point a session has received its allocated amount of Engines, any further matches are ignored.\cr
\bold{Invocation} Condition Use Affinity on an invocation basis carefully, as the Scheduler is optimized for Session-based affinity. Note that you also must set \code{Options.AFFINITY_DEPTH} to > 0 for it to have effect. When all Engines have been allocated to a Session during a schedule event, a score between each request and Engine is calculated based on the affinity for each allowed match. If a Discriminator Condition is set, Engines will be filtered out based on this. The score is zero if no affinity is set. If \code{Options.AFFINITY_DEPTH} is greater than 0, It will search up to that value or the number of Engines allocated to the Session, whichever is greater, and then assign tasks to Engines in order of affinity.\cr

  \describe{
	\item{SchedulingConditionFactory.createCustomAffinity(className, props)}{
	  Creates a \code{\link{Condition}} for Affinity from Engines to a Session or an Invocation request, using a \code{CustomAffinity} implementation.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{className}    \tab 	The name of the class. It must implement \code{CustomAffinity}, and the class must be deployed to the Broker in its lib directory.\cr
			\code{props}        \tab \code{\link{StringProperties}} that can be retrieved by the Affinity object. None of the properties can include non-ASCII or XML reserved characters.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The \code{\link{Condition}}. It will not be an object of specified class, as that is only available on the Broker, but \code{getProperties} is valid.
		  }
		}
	  }
	}
  }

  \describe{
	\item{SchedulingConditionFactory.createCustomDiscriminator(className, props)}{
	  Creates a \code{\link{Condition}} for Discrimination from Engines to a Session or an Invocation request, using a \code{CustomDiscriminator} implementation.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{className}    \tab 	The name of the class. It must implement \code{CustomDiscriminator}, and the class must be deployed to the Broker in its lib directory.\cr
			\code{props}        \tab \code{\link{StringProperties}} that can be retrieved by the Affinity object. None of the properties can include non-ASCII or XML reserved characters.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The \code{\link{Condition}}. It will not be an object of specified class, as that is only available on the Broker, but getProperties is valid.
		  }
		}
	  }
	}
  }

  \describe{
	\item{SchedulingConditionFactory.createPropertyAffinity(name, comparision, value, score)}{
	  Creates a new \code{PropertyAffinity} from Engines to a Session or an Invocation request.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{name}         \tab The name of the \code{Engine} property. It cannot include non-ASCII or XML reserved characters. Use property names defined in \code{\link{EngineProperties}}\cr
			\code{comparision}  \tab  One of the defined comparisons in \code{\link{Discriminator.EComparator}}, which compares the Engine's value against this.\cr
			\code{value}        \tab The property value. It cannot include non-ASCII or XML reserved characters.\cr
			\code{score}        \tab The score to add of the operation is true
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Affinity condition
		  }
		}
	  }
	}
  }

  \describe{
	\item{SchedulingConditionFactory.createPropertyDiscriminator(name, comparision, value, nullCompare)}{
	  Creates a new \code{PropertyDiscriminator} from Engines to a Session or an Invocation request.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{name}         \tab The name of the \code{Engine} property. It cannot include non-ASCII or XML reserved characters. Use property names defined in \code{\link{EngineProperties}}\cr
			\code{comparision}  \tab  One of the defined comparisons in \code{\link{Discriminator.EComparator}}, which compares the Engine's value against this.\cr
			\code{value}        \tab The property value. It cannot include non-ASCII or XML reserved characters.\cr
			\code{nullCompare}  \tab The result if the property does not exist.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Affinity condition
		  }
		}
	  }
	}
  }

  \describe{
	\item{SchedulingConditionFactory.createQueueJumpCondition()}{
	  Creates a Queue Jump Condition, which indicates that the submitted Task should be placed at the front of the queue of waiting tasks.\cr
	  This condition can only be added to a submit or execute call. The condition can be passed as that method's condition, or as a part of a set created by \code{ConditionSetFactory.createServiceSet()}
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			the queue jump condition
		  }
		}
      }
	}
  }

}
\examples{
  dualAffinity <- SchedulingConditionFactory.createPropertyAffinity(EngineProperties.CPU_NO, Discriminator.EComparator.EQUAL, "2", 4.0)
  preferFast <- SchedulingConditionFactory.createPropertyAffinity(EngineProperties.CURRENT_CPU, Discriminator.EComparator.GREATER_THAN, "200", 5.0)
  preferSlow <- SchedulingConditionFactory.createPropertyAffinity(EngineProperties.CURRENT_CPU, Discriminator.EComparator.LESS_THAN, "100", 5.0)
  isLinux <- SchedulingConditionFactory.createPropertyDiscriminator(EngineProperties.OS, Discriminator.EComparator.CONTAINS, "Linux", FALSE)

  # ---------------- Service ---------------------------------

  # First create a Service with dual affinity
  cs <- ServiceFactory.createService("MyService", list(), StringProperties(), StringProperties(), dualAffinity);
  rs <- ResultSet()

  # Now submit a Task that must run on a linux machine, and prefers > 200 mflops
  set <- ConditionSetFactory.createANDSet();
  set$add(isLinux);
  set$add(preferFast);

  # Title it
  set$add(DescriptorFactory.create("C++ Driver Linux Build"));
  id <- cs$submit(rs, "buildLinuxCPP", data, set);

  # Now submit a Task that must run on a linux machine, and prefers < 100 mflops
  set <- condSetFactory.createANDSet();
  set$add(preferSlow);
  discriminatorSet <- ConditionSetFactory.createANDSet();

  # this Task must also run on a linux machine
  # and it is dependent on the previous Task completing.
  # because it's a multiple of type discriminator/dependency, need a new set
  discriminatorSet$add(isLinux);
  discriminatorSet$add(DependencyFactory.createDependency(s$getId(), id, true));
  set$add(discriminatorSet);

  # Title it
  set$add(DescriptorFactory.create("C++ Engine Linux Build"))
  s$submit(rs, "buildLinuxEngine", data, set)

}
\seealso{
\code{\link{Service}}, 
\code{\link{DescriptorFactory}}, 
\code{\link{ConditionSetFactory}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
