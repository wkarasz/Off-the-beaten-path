\name{BasicEngineDiscriminator}
\alias{BasicEngineDiscriminator}
\alias{BasicEngineDiscriminator.new}
\alias{BasicEngineDiscriminator.setCpus}
\alias{BasicEngineDiscriminator.setExcludeUsers}
\alias{BasicEngineDiscriminator.setFreeMemInKB}
\alias{BasicEngineDiscriminator.setIncludeUsers}
\alias{BasicEngineDiscriminator.setMegaFlops}
\alias{BasicEngineDiscriminator.setOS}
\alias{BasicEngineDiscriminator.setTotalMemInKB}
\title{BasicEngineDiscriminator
}

\description{
  Deprecated: Use SchedulingConditionFactory.\cr
  Discriminates based on the default Engine properties. This is a convience class for using \code{PropertyDiscriminator} with \code{EngineProperties}.
}
\usage{
BasicEngineDiscriminator()
}

\details{
  All methods that require an instance of \code{BasicEngineDiscriminator} object can use the \code{$} method to invoke the desired method.\cr

  For example:\cr
  \code{
    obj$setCpus("2,4")
  }


  \describe{
	\item{BasicEngineDiscriminator.setCpus(self, cpus)}{
	  Comma-delimited list of processor numbers.\cr
	  For example, '2,4' would allow machines with 2 or 4 processors. Convenience method for \code{EngineProperties.CPU_NO}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{cpus}    \tab Comma-delimited list of processor numbers.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setExcludeUsers(self, usernames)}{
	  Comma-delimited list of users allowed.\cr
	  Convenience method for \code{EngineProperties.USERNAME}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}      \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{usernames} \tab Comma-delimited list of users allowed.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setFreeMemInKB(self, freeMem)}{
	  The minimum total memory in kilobytes.
	  Convenience method for \code{EngineProperties.FREE_MEMORY_KB}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{freeMem} \tab The minimum total memory in kilobytes as integer.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setIncludeUsers(self, usernames)}{
	  Comma-delimited list of users allowed.\cr
	  Convenience method for \code{EngineProperties.USERNAME}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}      \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{usernames} \tab Comma-delimited list of users allowed.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setMegaFlops(self, mflops}{
	  The minimum megaflop rating.
	  Convenience method for \code{EngineProperties.TOTAL_CPU}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{mflops}  \tab The minimum megaflop rating as numeric.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setOS(self, os)}{
	  Comma-delimited list of OS platforms allowed.
	  Supported OS platforms can be viewed on the Engine Install page.
	  Convenience method for \code{EngineProperties.OS}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{os}      \tab Comma-delimited list of OS platforms allowed.
		  }
		}
	  }
	}
  }

  \describe{
	\item{BasicEngineDiscriminator.setTotalMemInKB(self, totalMem)}{
	  The minimum free memory in kilobytes.
	  Convenience method for \code{EngineProperties.TOTAL_MEMORY_KB}.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}     \tab An instance of \code{BasicEngineDiscriminator}\cr
			\code{totalMem} \tab The minimum free memory in kilobytes as integer
		  }
		}
	  }
	}
  }
}

\seealso{
\code{\link{SchedulingConditionFactory}}
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
