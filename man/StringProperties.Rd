\name{StringProperties}
\alias{StringProperties}
\alias{StringProperties.new}
\alias{StringProperties.delete}
\alias{StringProperties.get}
\alias{StringProperties.put}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{StringProperties
}
\description{
StringProperties is a collection of name-value pair where both name and value are strings.
}
\usage{
StringProperties()
StringProperties.get(self, key)
StringProperties.put(self, key, value)
StringProperties.delete(self)
}
\details{
  This class is used as the return for some of the Job interfaces.
\cr
  All methods that require an instance of \code{StringProperties} object can use the \code{$} method to invoke the desired method.\cr

  For example:\cr
  \code{
    obj$put("key", "value")
  }

  \describe{
	\item{StringProperties()}{
	  Identical to \code{StringProperties.new}
	}
  }

  \describe{
	\item{StringProperties.new()}{
	  Creates a new \code{StringProperties} object.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			An instance of \code{StringProperties}
		  }
		}
	  }
	}
  }

  \describe{
	\item{StringProperties.delete(self)}{
	  Deletes the underlying C++ object.\cr
	  Any subsequent use of "free'd" objects yields to no-op.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{StringProperties}
		  }
		}
	  }
	}
  }

  \describe{
	\item{StringProperties.get(self, key)}{
	  Returns the value for the specified keys.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{StringProperties}\cr
			\code{key}          \tab the key
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the value. The string will be empty if no matching key was found.
		  }
		}
	  }
	}
  }

  \describe{
	\item{StringProperties.put(self, key, value)}{
	  Insert the key with the value.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}         \tab An instance of \code{StringProperties}\cr
			\code{key}          \tab the key\cr
			\code{value}        \tab the value. value with \code{boolean} semantic must use \code{"true"} and \code{"false"}.
		  }
		}
	  }
	}
  }

}
\examples{
	sp <- StringProperties()

	# add key/value
	sp$put("key", "value")

        # boolean values
	sp$put("key_true", "true")
	sp$put("key_false", "false")

	# retrieve value
	value <- sp$get("key")

	# call to delete() is not required. C++ objects are free'd when R objects are free'd
	delete(sp)
}
\seealso{
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
