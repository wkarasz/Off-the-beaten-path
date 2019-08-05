\name{ServiceGroup}
\alias{ServiceGroup.addService}
\alias{ServiceGroup.create}
\alias{ServiceGroup.createService}
\alias{ServiceGroup.destroyAll}
\alias{ServiceGroup.getId}
\alias{ServiceGroup.getName}
\alias{ServiceGroup.getService}
\alias{ServiceGroup.getServiceCount}
\alias{ServiceGroup}
\alias{ServiceGroup.waitUntilInactive}
\title{ServiceGroup
}
\description{
A convenience class for grouping Service Instances together.
}
\usage{
ServiceGroup(groupName)
ServiceGroup.create(groupName)
ServiceGroup.addService(self, service)
ServiceGroup.createService(self, serviceName, typeR = TRUE)
ServiceGroup.createService(self, serviceName, initData, typeR = TRUE)
ServiceGroup.createService(self, serviceName, initData, option, description, typeR = TRUE)
ServiceGroup.destroyAll(self)
ServiceGroup.getId(self)
ServiceGroup.getName(self)
ServiceGroup.getService(self, i)
ServiceGroup.getServiceCount(self)
ServiceGroup.waitUntilInactive(self)
}
\details{
  The group autogenerates a random unique ID.
  All \bold{Service} Instances created within this group will automatically assign the \code{Description.SERVICE_GROUP_ID} value to this ID,
  and the \code{Description.SERVICE_GROUP_NAME} value to the given name.\cr
  All methods that require an instance of \code{ServiceGroup} object can use the \code{$} method to invoke the desired method.\cr

  For example:  \cr
  \code{
	obj$addService(ServiceFactory.createService("serviceName"))
  }

  \describe{
	\item{ServiceGroup.create(groupName)}{
	  Creates a new \code{Service} group, and generates the ID.\cr
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{groupName} \tab The group name. If empty, the name is set to the ID
		  }
		}
	  }
	}

	\item{ServiceGroup.addService(self, service)}{Adds an existing instance to this group.\cr
	  Useful for when the Service cannot be directly created from this group, such as when using a generated proxy class.
	  You must create this Service with the \code{Description.SERVICE_GROUP_ID} and \code{Description.SERVICE_GROUP_NAME} values set to this Service Group's values.\cr
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}\cr
			\code{service} \tab An instance of \code{Service}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			TheSservice instance. The instance is deleted when the \code{ServiceGroup} is deleted.
		  }
		}
	  }
	}

	\item{ServiceGroup.createService(self, serviceName, typeR = TRUE)}{Creates a new Service instance in this group.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}        \tab An instance of \code{ServiceGroup}\cr
			\code{serviceName} \tab The name of the registered Service.\cr
			\code{typeR}       \tab if this Service exposed as \code{NATIVE_R}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The service instance. The instance is deleted when the \code{ServiceGroup} is deleted.
		  }
		}
	  }
	}
	\item{ServiceGroup.createService(self, serviceName, initData, typeR = TRUE)}{Creates a new Service instance in this group.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}        \tab An instance of \code{ServiceGroup}\cr
			\code{serviceName} \tab The name of the registered service.\cr
			\code{initData}    \tab The initialization data\cr
			\code{typeR}       \tab if this Service is exposed as \code{NATIVE_R}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The Service instance. The instance is deleted when the \code{ServiceGroup} is deleted.
		  }
		}
	  }
	}
	\item{ServiceGroup.createService(self, serviceName, initData, option, description, typeR = TRUE)}{Creates a new Service instance in this group.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}        \tab An instance of \code{ServiceGroup}\cr
			\code{serviceName} \tab The name of the registered Service.\cr
			\code{initData}    \tab The initialization data\cr
			\code{option}      \tab Additional Service options. See \code{\link{Options}}. This will not override the registered options.\cr
			\code{description} \tab Additional Service description itmes. See \code{\link{Description}}. This will not override the registered descriptions.\cr
			\code{typeR}       \tab if this Service is exposed as \code{NATIVE_R}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The Service instance. The instance is deleted when the \code{ServiceGroup} is deleted.
		  }
		}
	  }
	}
	\item{ServiceGroup.destroyAll(self)}{Destroys all instances in this group.\cr
	  If instances were created with the same group ID in a separate ServiceGroup, they will not be destroyed.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}
		  }
		}
	  }
	}
	\item{ServiceGroup.getId(self)}{Gets the generated ID.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the group id
		  }
		}
	  }
	}
	\item{ServiceGroup.getName(self)}{Gets the group name.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the name
		  }
		}
	  }
	}
	\item{ServiceGroup.getService(self, i)}{Returns the i'th instance in this group.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}\cr
			\code{i}       \tab 	the index
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the instance
		  }
		}
	  }
	}
	\item{ServiceGroup.getServiceCount(self)}{Gets the number of Service instances in this group.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			the count
		  }
		}
	  }
	}
	\item{ServiceGroup.waitUntilInactive(self)}{Waits until all instances are inactive
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{self}    \tab An instance of \code{ServiceGroup}
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
