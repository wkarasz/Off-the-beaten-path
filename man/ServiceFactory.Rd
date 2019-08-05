\name{ServiceFactory}
\alias{ServiceFactory}
\alias{ServiceFactory.createService}
\title{ServiceFactory}
\description{
  A factory for creating Service Instances.
}

\usage{
  ServiceFactory.createService(serviceName, typeR = TRUE)
  ServiceFactory.createService(serviceName, initData, typeR = TRUE)
  ServiceFactory.createService(serviceName, initData, options, description, typeR = TRUE)
  ServiceFactory.createService(serviceName, initData, options, description, condition, typeR = TRUE)
}
\details{
  A factory for creating Service Instances.\cr
  \code{Service} instances does not have to be explicitly deleted.

  \describe{
	\item{ServiceFactory.createService(serviceName, typeR = TRUE)}{
	  Creates a new \code{Service} instance.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{serviceName}              \tab The name of the registered \code{\link{Service}}\cr
			\code{typeR}                    \tab If the Service is implemented in R and is using the NATIVE_R interface
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Service.
		  }
		}
	  }
	}
  }

  \describe{
	\item{ServiceFactory.createService(serviceName, initData, typeR = TRUE)}{
	  Creates a new \code{Service} instance.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{serviceName}              \tab The name of the registered \code{\link{Service}}\cr
			\code{initData}                 \tab The state initialization data as list(), passed into the \code{Options.INIT_METHOD} if set. initData is not allowed when creating a Shared Service.\cr
			\code{typeR}                    \tab If the Service is implemented in R and is using the NATIVE_R interface
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Service.
		  }
		}
	  }
	}
  }

  \describe{
	\item{ServiceFactory.createService(serviceName, initData, options, description, typeR = TRUE)}{
	  Creates a new \code{Service} instance.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{serviceName}              \tab The name of the registered \code{\link{Service}}\cr
			\code{initData}                 \tab The state initialization data as list(), passed into the \code{Options.INIT_METHOD} if set. initData is not allowed when creating a Shared Service.\cr
			\code{options}                  \tab \code{\link{Service}} options. This will not override the binding. \code{options} is of type \code{\link{StringProperties}}\cr
			\code{description}              \tab \code{\link{Service}} description. This will not override the binding. \code{description} is of type \code{\link{StringProperties}} \cr
			\code{typeR}                    \tab If the Service is implemented in R and is using the NATIVE_R interface
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Service.
		  }
		}
	  }
	}
  }

  \describe{
	\item{ServiceFactory.createService(serviceName, initData, options, description, condition, typeR = TRUE)}{
	  Creates a new \code{Service} instance.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{serviceName}              \tab The name of the registered \code{Service}\cr
			\code{initData}                 \tab The state initialization data as list(), passed into the \code{Options.INIT_METHOD} if set. initData is not allowed when creating a Shared Service.\cr
			\code{options}                  \tab \code{\link{Service}} options. This will not override the binding. \code{options} is of type \code{\link{StringProperties}}\cr
			\code{description}              \tab \code{\link{Service}} description. This will not override the binding. \code{description} is of type \code{\link{StringProperties}} \cr
			\code{condition}                \tab A \code{\link{Condition}} for any Invocation to be executed.\cr
			\code{typeR}                    \tab If the Service is implemented in R and is using the NATIVE_R interface
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the Service.
		  }
		}
	  }
	}
  }
}
\examples{
	# "myservice" using NATIVE_R interface and without a initialization method
	cs <- ServiceFactory.createService("myservice")

	# "myservice" using non-NATIVE_R interface and without a initialization method
	cs <- ServiceFactory.createService("myservice", FALSE)

	# "myservice" using NATIVE_R interface with a initialization method init(f, g)
	cs <- ServiceFactory.createService("myservice", list(f, g))


	#
	options <- StringProperties()
	description <- StringProperties()
	condition <- SchedulingConditionFactory.createPropertyAffinity(EngineProperties.OS, ....)
	cs <- ServiceFactory.createService("myservice", list(f, g), options, description, condition)
}
\seealso{
\code{\link{StringProperties}}, \code{\link{SchedulingConditionFactory}}, \code{\link{EngineProperties}}, \code{\link{Options}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
