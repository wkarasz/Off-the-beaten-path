\name{Service}
\alias{Service}
\alias{Service.cancelInvocation}
\alias{Service.destroy}
\alias{Service.destroyWhenInactive}
\alias{Service.execute}
\alias{Service.flush}
\alias{Service.getDescription}
\alias{Service.getId}
\alias{Service.getInvocationCount}
\alias{Service.getOption}
\alias{Service.isDestroyed}
\alias{Service.submit}
\alias{Service.updateState}
\alias{Service.waitUntilInactive}
\title{Service}
\description{
The base Service client class.
}
\usage{
Service.cancelInvocation(self, invocationId)
Service.destroy(self)
Service.destroyWhenInactive(self)
Service.execute(self, methodname, data, condition = NULL)
Service.flush(self)
Service.getDescription(self, name)
Service.getId(self)
Service.getInvocationCount(self)
Service.getOption(self, name)
Service.isDestroyed(self)
Service.submit(self, rs, methodname, data, condition = NULL)
Service.updateState(self, methodname, data, append)
Service.waitUntilInactive(self, millis = as.integer(60 * 1000))
}
\details{
  A \code{Service} enables dispatching and distribution of remote services across a Grid. These services can be implemented in a language independent manner; for example, a Java client can call a .NET \code{Service}. The only requirement is that the \code{Service} implementation methods follow allowed formats. Any \code{Service} can be implemented with methods that take String input and return String output, and byte[] input / byte[] output. Additionally, a Java or .NET service can perform serialized object I/O from a client implemented in its respective language.\cr
  A \code{Service} implemented in R and exposed as \code{NATIVE_R} can take native R types as input and return. Closure is \bold{not} supported.\cr
  All methods that require an instance of \code{Service} object can use the \code{$} method to invoke the desired method.\cr

  For example:  \cr
  \code{
	r <- cs$execute("myfunction", list())
  }

  \describe{
	\item{Service.cancelInvocation(self, invocationId)}{
	  Cancels an invocation specified by the given invocation ID.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{invocationId} \tab Cancels an invocation specified by the given invocation ID.
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.destroy(self)}{
	  Destroys this Service.\cr
	  This refers to deleting resources on the Engines and the Broker; it is the client's responsibility to delete the \code{Service}. This method has different implications for the various collection types.
	  \code{Options.Collection.IMMEDIATELY} (default): If there are no outstanding requests, the Session closes as 'Finished'. If there are, they will be cancelled, and the Session will close as 'Cancelled'\cr
	  \code{Options.Collection.LATER}: Must be called when the Instance is finished submitting, to indicate that results can now be collected from this Session.\cr
	  \code{Options.Collection.NEVER}: Indicates to the local Instance that submission is complete, so that it can free local resources associated with the session.\cr
	  ServiceFactory.getService(): Destroys the Session so that no more Instances can collect from from the Session.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.destroyWhenInactive(self)}{
	  Indicates that this Service will be destroyed as soon as all pending tasks have completed.\cr
	  This refers to deleting resources on the Engines and the Broker; it is the client's responsibility to delete the \code{Service}. This method returns immediately. Typically this is used to run a \code{Service} like a Job. That is, the \code{Service} is started, a set of tasks are submitted, then this call is made. addLifecycleHandler can be used for notification upon completion. Tasks may be added up until the session becomes inactive; typically this would only be done in the case of a failed task that the client decides to resubmit.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.execute(self, methodname, data, condition = NULL)}{
	  ynchronously invokes the method on the \code{Service}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{methodname}   \tab The method to call\cr
			\code{input}        \tab The input data\cr
			\code{condition}    \tab A condition for the invocation to be executed. May be NULL, meaning no condition
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The output data
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.flush(self)}{
	  Immediately sends any buffered task submissions that have not yet been sent.\cr
	  When \code{Options.INVOCATIONS_PER_MESSAGE} > 1, this should be called after submitting tasks if you do not intend to call \code{Service.waitUntilInactive()} or \code{Service.destroyWhenInactive()} immediately (both of which also flush submissions). Otherwise the submissions will not be sent until the next auto-flush poll.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.getDescription(self, name)}{
	  Gets a description property, used for discrimination and accounting. See \code{\link{Description}}
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{name}         \tab The name
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			A description
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.getId(self)}{
	  \describe{
		Get the ID of the service.
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			A unique identifier
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.getInvocationCount(self)}{
	  Gets the count of outstanding invocations.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the count
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.getOption(self, name)}{
	  Gets a service option.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{name}         \tab The option name 
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The option value
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.isDestroyed(self)}{
	  Whether the Service is destroyed.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			TRUE if destroyed
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.submit(self, rs, methodname, data, condition = NULL)}{
	  Asynchronously invokes the method on the Service.\cr
	  When using asynchronous submission, under most circumstances you should set \code{Options.INVOCATIONS_PER_MESSAGE} > 1 to improve performance. Also, you may want to consider \code{Options.Collection.AFTER_SUBMIT} depending on your architecture.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{rs}           \tab An instance of \code{ResultSet}\cr
			\code{methodnamee}  \tab The method to call\cr
			\code{data}         \tab The input data\cr
			\code{condition}    \tab A condition for the invocation to be executed. May be null, meaning no condition
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The handler ID
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.updateState(self, methodname, data, append)}{
	  Updates the state of the Service.\cr
	  This call is 'lazy', in that it is not called on any Engine until the Engine takes a Service invocation, at which time the update method invocation is made immediately prior to the service invocation. If append is set to false, the list of update calls for this service is reset, and an Engine will ignore any previous update calls if it has not yet made them. If true, it will execute all update calls it has not yet made in order, starting with the last non-appended call.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{methodnamee}  \tab The method to call\cr
			\code{data}         \tab The input data\cr
			\code{append}       \tab If calls should be appended
		  }
		}
	  }
	}
  }

  \describe{
	\item{Service.waitUntilInactive(self, millis = as.integer(60 * 1000))}{
	  Waits until there are no outstanding active Services.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{Service}\cr
			\code{millis}       \tab Timeout, in milliseconds. 0 means forever.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			\code{TRUE} if inactive, \code{FALSE} if not (always \code{TRUE} if \code{millis} == 0)
		  }
		}
	  }
	}
  }

}
\examples{
	cs <- ServiceFactory.createService("myservice")

	# call myfunction()
	cs$execute("myfunction", list())

	# creates a ResultSet object to collect the responses from asynchronous calls
	rs <- ResultSet(cs)

	# call myfunction() asynchronously
	# equivalent to Service.submit(cs, rs, "myfunction", list())
	cs$submit(rs, "myfunction", list())

	# wait at most for 10 secs.
	# equivalent to: LL <- ResultSet.waitForAll(rs, as.integer(10 * 1000))
	LL <- rs$waitForAll(as.integer(10 * 1000))

	# wait for all tasks to complete
	cs$waitUntilInactive(as.integer(60 * 1000))

	cs$destroy()

	# call to delete() is not required. C++ objects are free'd when R objects are free'd
	delete(rs)
	delete(cs)
}
\seealso{
\code{\link{ResultSet}},
\code{\link{SchedulingConditionFactory}},
\code{\link{ServiceFactory}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
