\name{ResultSet}
\alias{ResultSet}
\alias{ResultSet.getResult}
\alias{ResultSet.waitFor}
\alias{ResultSet.waitForAll}
\title{ResultSet
}
\description{
  ResultSet handles the asynchronous notification when tasks are completed.\cr
}

\usage{
  ResultSet(cs)
  ResultSet.getResult(self, id)
  ResultSet.waitFor(self, id, millis = as.integer(60 * 1000))
  ResultSet.waitForAll(self, millis = as.integer(300 * 1000))
  delete(rs)
}
\details{
  Results from asynchronous notifications are accumulated until they are retrieved by the code in R.  \cr
  All methods that require an instance of \code{ResultSet} object can use the \code{$} method to invoke the desired method.\cr
  For example:  \cr
  \code{
    rs$waitForAll()
  }

  \describe{
	\item{delete(rs)}{
	  Deletes the underlying C++ object. Calling \code{delete} is not required. Underlying C++ objects are freed when R GC cleans up the object.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{rs}         \tab An instance of \code{ResultSet}
		  }
		}
	  }
	}
  }

  \describe{
	\item{ResultSet(cs)}{
	  Creates a new \code{ResultSet} object.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{cs}         \tab An instance of \code{Service}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			An instance of \code{ResultSet}
		  }
		}
	  }
	}
  }

  \describe{
	\item{ResultSet.getResult(self, id)}{
	  Waits for the result of \code{id} to complete or to the max wait time to expire.\cr
	  Return value represents the current state of submission, which may not have completed yet.\cr
	  The returned value is a \code{list} with following attributes:\cr
	  \itemize{
		\item{\code{taskId}}
		\item{\code{status}}
		\item{\code{value}}
	  }
	  There are 3 possible conditions reported by the \code{status}:\cr
	  \itemize{   
		\item{task was submitted but there is no response yet.}
		\item{task was submitted and has failed.}
		\item{task was submitted and successfully.}
	  }
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}       \tab An instance of \code{ResultSet}\cr
			\code{id}         \tab \code{id} returned from \code{Service.submit}
		  }
		}
		\item{Returns}{
		  Current state of submission:
		  \tabular{ll}{
			\code{taskId} \tab {taskId returned by \code{\link{Service.submit}}}\cr
			\code{status} \tab {0: completed successfully, -1: failed, 1: submitted but no response yet.}\cr
			\code{value}  \tab {(status == 0): returned value from the Service, (status == -1): error message}
		  }
		}
	  }
	}
  }

  \describe{
	\item{ResultSet.waitFor(self, id, millis)}{
	  Waits for the result of \code{id} to complete or to the max wait time to expire.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}       \tab An instance of \code{ResultSet}\cr
			\code{id}         \tab \code{id} returned from \code{Service.submit}\cr
			\code{millis}     \tab Timeout, in milliseconds.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			Returns the result of submission as described in \code{ResultSet.getResult}.
		  }
		}
	  }
	}
  }

  \describe{
	\item{ResultSet.waitForAll(self, millis)}{
	  Waits for all submissions to complete or to the max wait time to expire.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}       \tab An instance of \code{ResultSet}\cr
			\code{millis}     \tab Timeout, in milliseconds.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			Returns a list of result of each submission as described in \code{ResultSet.getResult}.
		  }
		}
	  }
	}
  }

}
\examples{
	cs <- ServiceFactory.createService("myservice")

	# creates a ResultSet object to collect the responses from asynchronous calls
	rs <- ResultSet(cs)

	# call myfunction() asynchronously
	# equivalent to Service.submit(cs, rs, "myfunction", list())
	id <- cs$submit(rs, "myfunction", list())

	# wait at most for 10 secs.
	# equivalent to: r <- ResultSet.waitFor(rs, id, as.integer(10 * 1000))
	r <- rs$waitFor(id, as.integer(10 * 1000))

	# 
	# There are 3 attributes in r:
	# - r$taskId
	# - r$status
	# - r$value
	#
	# refer to getResult() for more details.
	#

	#
	# Service.submit() followed by ResultSet.waitFor()
	#
	#   id <- cs$submit(rs, "myfunction", list()
	#   r <- rs$waitFor(id, as.integer(10 * 1000))
	#
	# is almost the same as:
	#
	#   r <- cs$execute("myfunction", list())
	#
	# Service.execute() blocks until the test terminates.
	#

}


\seealso{
\code{\link{Service}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
