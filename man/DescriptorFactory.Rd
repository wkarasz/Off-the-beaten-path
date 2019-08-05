\name{DescriptorFactory}
\alias{DescriptorFactory}
\alias{DescriptorFactory.create}
\title{DescriptorFactory
}
\description{
  A factory for creating a Descriptor \code{Condition}, which is used to display descriptive information about an invocation request.
}
\usage{
  DescriptorFactory.create(information)
}
\details{
  A factory for creating a Descriptor \code{Condition}, which is used to display descriptive information about an invocation request.\cr
  This information is available on the Task Admin page and Admin API. This enables the user to set descriptive information from the Client, as opposed to using com.datasynapse.gridserver.engine.ServiceInvocationContext.updateInfo(), which sets it from the Engine. Also, unlike the Task Info, this information cannot be updated after it is set. By default, the information is not stored in the reporting database. However, you can choose to store it rather than the Task Information in the database on a Broker-wide level using the associated Database configuration setting.

  \describe{
	\item{DescriptorFactory.create(information)}{
	  Creates an information object.\cr
	  This object should be added to the context set for the request. This \code{Condition} may not be used in a set created by \code{ConditionSet.createORSet}, and can only be in the top level of a set created by \code{ConditionSet.createANDSet}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{information}   \tab The information string. It cannot include non-ASCII or XML reserved characters. Note that if this is longer than the column size of your database, the value will either be truncated or discarded. See your database documentation for more details.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The context object
		  }
		}
	  }
	}
  }
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
