\name{DriverManager}
\alias{DriverManager}
\alias{DriverManager.setProperty}
\alias{DriverManager.getProperty}
\alias{DriverManager.getConfigProperties}
\alias{DriverManager.getUserProperties}
\alias{DriverManager.addProperties}
\alias{DriverManager.connect}
\alias{DriverManager.disconnect}
\alias{DriverManager.getCurrentBroker}
\alias{DriverManager.getCurrentBrokers}
\alias{DriverManager.isConnected}
\alias{DriverManager.BROKER_LIST}
\alias{DriverManager.BROKER_TIMEOUT}
\alias{DriverManager.BROKER_TIMEOUT_DEFAULT}
\alias{DriverManager.CONSOLE_LOG_LEVEL}
\alias{DriverManager.CURRENT_USER_PASSWORD}
\alias{DriverManager.DATA_REFERENCE_FILE_EXPIRATION}
\alias{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_DEFAULT}
\alias{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_FREQUENCY}
\alias{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_FREQUENCY_DEFAULT}
\alias{DriverManager.DDT_ENABLED}
\alias{DriverManager.DDT_ENABLED_DEFAULT}
\alias{DriverManager.DDT_FILE_EXPIRATION}
\alias{DriverManager.DDT_FILE_EXPIRATION_DEFAULT}
\alias{DriverManager.DDT_FILE_EXPIRATION_FREQUENCY}
\alias{DriverManager.DDT_FILE_EXPIRATION_FREQUENCY_DEFAULT}
\alias{DriverManager.DDT_TYPE}
\alias{DriverManager.DDT_TYPE_DEFAULT}
\alias{DriverManager.FAILOVER_TIMEOUT}
\alias{DriverManager.FAILOVER_TIMEOUT_DEFAULT}
\alias{DriverManager.FILE_LOG_LEVEL}
\alias{DriverManager.FILESERVER_BASE_DIR}
\alias{DriverManager.FILESERVER_BASE_DIR_DEFAULT}
\alias{DriverManager.FILESERVER_BASE_URL}
\alias{DriverManager.FILESERVER_BASE_URL_DEFAULT}
\alias{DriverManager.FILESERVER_KEEP_ALIVE_THREADS}
\alias{DriverManager.FILESERVER_KEEP_ALIVE_THREADS_DEFAULT}
\alias{DriverManager.FILESERVER_KEEP_ALIVE_TIME}
\alias{DriverManager.FILESERVER_KEEP_ALIVE_TIME_DEFAULT}
\alias{DriverManager.FILESERVER_PORT}
\alias{DriverManager.FILESERVER_PORT_DEFAULT}
\alias{DriverManager.FILESERVER_PORT_RANGE}
\alias{DriverManager.FILESERVER_PORT_RANGE_DEFAULT}
\alias{DriverManager.FILESERVER_TYPE}
\alias{DriverManager.FILESERVER_TYPE_DEFAULT}
\alias{DriverManager.GRIDCACHE_MAX_MEM_IN_MB}
\alias{DriverManager.GRIDCACHE_MAX_MEM_IN_MB_DEFAULT}
\alias{DriverManager.JOB_TIMEOUT}
\alias{DriverManager.LOCAL_IP_ADDRESS}
\alias{DriverManager.LOG_DIR}
\alias{DriverManager.LOG_DIR_DEFAULT}
\alias{DriverManager.LOG_EXPIRATION_DEFAULT}
\alias{DriverManager.LOG_EXPIRATION_HOURS}
\alias{DriverManager.LOG_LENGTH}
\alias{DriverManager.LOG_LENGTH_DEFAULT}
\alias{DriverManager.LOG_LEVEL}
\alias{DriverManager.LOG_LEVEL_DEFAULT}
\alias{DriverManager.LOG_STACK_TRACE}
\alias{DriverManager.LOG_STACK_TRACE_DEFAULT}
\alias{DriverManager.LOG_TIMEZONE}
\alias{DriverManager.LOG_TIMEZONE_DEFAULT}
\alias{DriverManager.LOG_TO_FILE}
\alias{DriverManager.LOG_TO_FILE_DEFAULT}
\alias{DriverManager.LOG_TO_SYSTEM}
\alias{DriverManager.MAX_CONCURRENT_COLLECT}
\alias{DriverManager.MAX_CONCURRENT_COLLECT_DEFAULT}
\alias{DriverManager.MAX_CONCURRENT_JOBS}
\alias{DriverManager.MAX_CONCURRENT_JOBS_DEFAULT}
\alias{DriverManager.NEGOTIATE_ENABLED}
\alias{DriverManager.PASSWORD}
\alias{DriverManager.PRIMARY_DIRECTOR}
\alias{DriverManager.SECONDARY_DIRECTOR}
\alias{DriverManager.SSL_TRUST_FILE}
\alias{DriverManager.SSL_VERIFY_HOSTNAME}
\alias{DriverManager.USERNAME}
\alias{Options.AFFINITY_DEPTH}
\alias{Options.AFFINITY_WAIT}
\alias{Options.ALLOWED_DRIVER_PROFILES}
\alias{Options.ALLOWED_ROLES}
\alias{Options.AUTO_PACK_MODE}
\alias{Options.AUTO_PACK_NUM}
\alias{Options.BLACKLIST_MACHINE}
\alias{Options.COLLECTION_TYPE}
\alias{Options.ENCRYPTION_ENABLED}
\alias{Options.FAILOVER_SHARED_SERVICE_MIN_TASK_ID}
\alias{Options.INVOCATIONS_PER_MESSAGE}
\alias{Options.RESULTS_PER_MESSAGE}
\alias{Options.SHARED_SERVICE_INACTIVITY_KEEP_ALIVE_TIME}
\alias{Options.SHARED_SERVICE_NAME}
\title{DriverManager}
\description{
Allows for runtime setting of Driver properties, both configuration and user-defined.
}
\usage{
DriverManager.setProperty(name, value)
value <- DriverManager.getProperty(name)
props <- DriverManager.getConfigProperties()
props <- DriverManager.getUserProperties()
DriverManager.addProperties(props)
DriverManager.connect()
DriverManager.disconnect(broker)
broker <- DriverManager.getCurrentBroker()
list <- DriverManager.getCurrentBrokers()
value <- DriverManager.isConnected()
}
\details{
Allows for runtime setting of Driver properties, both configuration and user-defined.\cr
The properties set here will override any properties set in the \code{driver.properties} file. If any of these values are not set in the driver.properties file and not set here, defaults will be used. All configuration property names and default values are specified by the static Strings.

  \describe{
	\item{DriverManager.addProperties(props)}{
	  Adds the properties, user-defined or configuration.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{props}         \tab The properties
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.connect(manager)}{
	  Specifies that the Driver should connect to the specified Manager (by Name or ID), or to clear the specification if blank.\cr
          Typically this is used to connect to a Broker by name. It can also be used to connect to a Standalone Director by Manager ID, for the sake of using the Admin API. If currently connected, it will log off and then attempt to connect to the specified Manager. This method can only be called if the Driver is idle; if there are any services running, it will not logoff, and execution stops.
	  \describe{
		\item{Arguments}{
		  \tabular{ll}{
			\code{manager}         \tab The Broker name, the Standalone Director ID, or blank. If blank, it clears the specified Broker so that it will log into any Broker that it is normally allowed.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.connect()}{
	  Initiates the connection of the Driver to the Manager.\cr
Since instantiating a Job will connect the Driver if the Driver has not yet initialized, it is only necessary to call this method if only intending to use the Managers and not executing Jobs.
	}
  }



  \describe{
	\item{DriverManager.disconnect()}{
	  Disconnects this Driver from the current Manager(s) and releases system resources.\cr
This will first cleanly log this Driver off of the current Manager(s). It will then shut down all threads, and release nearly all memory (except a small amount of static data that will be re-used upon reconnect). This call has no effect if the Driver is already disconnected. Any operation that reconnects the Driver, such as connect(), creating a new \code{\link{Service}}, etc., will bring up a new Driver session.
	}
  }

  \describe{
	\item{DriverManager.getConfigProperties()}{
	  Gets all configuration properties.\cr
Only valid after the Driver has initialized.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			The properties
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.getCurrentBroker()}{
	  The current logged in Broker.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			the Broker name if the Driver is connected, or blank otherwise.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.getCurrentBrokers()}{
	  The list of logged in Brokers.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			all Broker names if the Driver is connected to.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.getProperty(name)}{
	  \describe{
		Gets a property, configuration or user-defined.
		\item{Arguments}{
		  \tabular{ll}{
			\code{name}         \tab The property, can be a defined configuration property, or a user-defined property for Broker discrimination.
		  }
		}
		\item{Returns}{
		  \tabular{l}{
			The property value, or null if not defined.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.getUserProperties()}{
	  Gets all user-defined properties.\cr
Only valid after the Driver has initialized.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			The properties.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.isConnected()}{
	  The connection status.
	  \describe{
		\item{Returns}{
		  \tabular{l}{
			\code{TRUE} if the Driver is connected to a Broker, \code{FALSE} otherwise.
		  }
		}
	  }
	}
  }

  \describe{
	\item{DriverManager.setProperty(name, name)}{
	  \describe{
		Sets a property, configuration or user-defined.\cr
This can only be called when the Driver is disconnected.
		\item{Arguments}{
		  \tabular{ll}{
			\code{name}         \tab The property, can be a defined configuration property, or a user-defined property for Broker discrimination.
			\code{value}        \tab The property value
		  }
		}
	  }
	}
  }

}
\value{
\item{DriverManager.BROKER_LIST}{Semicolon-delimited list of Broker names.\cr If set, Broker Spanning will be used for this Driver.}
\item{DriverManager.BROKER_TIMEOUT}{The timeout in seconds used by the Driver to determine that the current Broker is unavailable.\cr After the Driver has been logged off for this amount of time, it will reconnect to any Broker and resubmit any outstanding work.}
\item{DriverManager.BROKER_TIMEOUT_DEFAULT}{ \code{"300"}}
\item{DriverManager.CONSOLE_LOG_LEVEL}{Level of console messages.\cr
(deprecated, use LOG_LEVEL instead) \cr
1: Initialization only \cr
2: Driver events \cr
3: 2 + Engine events \cr
4: 3 + Job events \cr
5: 4 + high-level debug messages \cr
6: 5 + low-level debug messages}
\item{DriverManager.CURRENT_USER_PASSWORD}{The password of the current user running the Driver.\cr
This is only used for the Run-as feature; see the Administration Guide for more details on using Run-as. This value MUST be set if running a \code{Service} as this user. It can only be set prior to Driver initialization. For security, null is always returned from \code{DriverManager.getProperty(DriverManager.CURRENT_USER_PASSWORD)}.}
\item{DriverManager.DATA_REFERENCE_FILE_EXPIRATION}{The amount of time Data Reference files will persist before being deleted.}
\item{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_DEFAULT}{\code{"120"}}
\item{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_FREQUENCY}{The number of times per day that the Driver will check for and clean expired Data Reference files.\cr
Set to \code{"0"} to disable cleaning of these files.}
\item{DriverManager.DATA_REFERENCE_FILE_EXPIRATION_FREQUENCY_DEFAULT}{\code{"24"}}
\item{DriverManager.DDT_ENABLED}{Whether data is transferred directly to the Engines from this Driver.\cr
\code{"true"}: Engines will download data directly from the Driver via the local fileserver. \cr
\code{"false"}: Drivers will send data to the Broker, from which Engines will download it.
}
\item{DriverManager.DDT_ENABLED_DEFAULT}{\code{"true"}}
\item{DriverManager.DDT_FILE_EXPIRATION}{The amount of time DDT files will persist before being deleted.\cr
Under normal conditions, files and directories are deleted when they are no longer needed. Files used by a currently executing \code{Service} will not be deleted regardless of this value or how old they are. This cleaner will delete files left over by a \code{Service} that did not complete normally.}
\item{DriverManager.DDT_FILE_EXPIRATION_DEFAULT}{\code{"120"}}
\item{DriverManager.DDT_FILE_EXPIRATION_FREQUENCY}{The number of times per day that the Driver will check for and clean expired DDT files.
Set to \code{"0"} to disable cleaning of these files.}
\item{DriverManager.DDT_FILE_EXPIRATION_FREQUENCY_DEFAULT}{}
\item{DriverManager.DDT_TYPE}{  Deprecated:\cr
No longer used, use \code{Options.SHARED_UNIX_DIR/SHARED_WIN_DIR}
}
\item{DriverManager.DDT_TYPE_DEFAULT}{\code{"http"}}
\item{DriverManager.FAILOVER_TIMEOUT}{The timeout in seconds used by the Driver to determine that no Brokers are available and all services should fail.\cr
If the Driver cannot log in to any Broker for this amount of time after the BROKER_TIMEOUT has expired, it assumes that no Brokers will become available, and will cancel any outstanding work and throw an exception if any new work is started. A negative number means never time out.}
\item{DriverManager.FAILOVER_TIMEOUT_DEFAULT}{\code{"-1"}\cr Any negative number means never timeout.}
\item{DriverManager.FILE_LOG_LEVEL}{Level of log messages.\cr
(deprecated, use LOG_LEVEL instead) \cr
1: Initialization only \cr
2: Driver events \cr
3: 2 + Engine events \cr
4: 3 + Job events \cr
5: 4 + high-level debug messages \cr
6: 5 + low-level debug messages\cr
}
\item{DriverManager.FILESERVER_BASE_DIR}{The directory in which fileserver data is placed.\cr
This property can only be changed before Driver initialization.
}
\item{DriverManager.FILESERVER_BASE_DIR_DEFAULT}{\code{"./data"}}
\item{DriverManager.FILESERVER_BASE_URL}{If using "external", the URL which maps to the \code{FILESERVER_BASE_DIR}.\cr
Ignored otherwise\cr
This property can only be changed before Driver initialization.}
\item{DriverManager.FILESERVER_BASE_URL_DEFAULT}{\code{"none"}}
\item{DriverManager.FILESERVER_KEEP_ALIVE_THREADS}{Number of threads to wait for keep alive connections.\cr
default: \code{"200"}
}
\item{DriverManager.FILESERVER_KEEP_ALIVE_THREADS_DEFAULT}{\code{"200"}}
\item{DriverManager.FILESERVER_KEEP_ALIVE_TIME}{Number of seconds an idle file server keep-alive connection should listen to the socket.\cr
Set to \code{"0"} to disable keep-alive connection. default: \code{"60"}}
\item{DriverManager.FILESERVER_KEEP_ALIVE_TIME_DEFAULT}{\code{"60"}}
\item{DriverManager.FILESERVER_PORT}{The port of the fileserver.\cr
If \code{"0"}, any available port will be used on the internal fileserver. 
This property can only be changed before Driver initialization.
}
\item{DriverManager.FILESERVER_PORT_DEFAULT}{\code{"0"}}
\item{DriverManager.FILESERVER_PORT_RANGE}{The port range of the internal fileserver, in the format [startPort]-[endPort].\cr
For example, \code{3456}-\code{3466}. If specified, property \code{FILESERVER_PORT} will be ignored. The first available port in range is used. This property can only be changed before Driver initialization.}
\item{DriverManager.FILESERVER_PORT_RANGE_DEFAULT}{\code{"0"}}
\item{DriverManager.FILESERVER_TYPE}{Type of fileserver to use for DDT.\cr
\code{"none"}: no fileserver, You cannot use DDT or StreamJobs in this case, so \code{DDT_TYPE} must be \code{"false"} \cr
\code{"internal"}: uses the included internal lightweight fileserver. \cr
\code{"external"}: connects to another webserver that is already running on this machine. If your application already uses a fileserver such as Apache or Resin, you may want to link to it rather than using the internal webserver.\cr
This property can only be changed before driver initialization.}
\item{DriverManager.FILESERVER_TYPE_DEFAULT}{\code{"internal"}}
\item{DriverManager.GRIDCACHE_MAX_MEM_IN_MB}{\code{"none"}}
\item{DriverManager.GRIDCACHE_MAX_MEM_IN_MB_DEFAULT}{\code{"none"}}
\item{DriverManager.JOB_TIMEOUT}{       
Deprecated:\cr
Use \code{BROKER_TIMEOUT}
}
\item{DriverManager.LOCAL_IP_ADDRESS}{Set this to override the IP address used when other clients connect to the Driver for DDT.\cr
A typical case is when there is more than one network adapter and you do not want to use the default adapter.
}
\item{DriverManager.LOG_DIR}{   
The log directory.\cr
This should be set to a directory that is not shared by any other application. For Windows, forward slashes are translated to back slashes, and the directory root is the current Windows drive. For example, the directory \code{/sdk/log} should be used instead of \code{c:}.}
\item{DriverManager.LOG_DIR_DEFAULT}{\code{./logs}}
\item{DriverManager.LOG_EXPIRATION_DEFAULT}{\code{120}}
\item{DriverManager.LOG_EXPIRATION_HOURS}{Log files older than this are automatically deleted. }
\item{DriverManager.LOG_LENGTH}{The length, in characters, after which the log file will rollover to a new file.}
\item{DriverManager.LOG_LENGTH_DEFAULT}{\code{1000000}}
\item{DriverManager.LOG_LEVEL}{
Level of log messages.\cr
This value superscedes both FILE and CONSOLE LOG_LEVEL \cr
1: Initialization only \cr
2: Driver events \cr
3: 2 + Engine events \cr
4: 3 + Job events \cr
5: 4 + high-level debug messages \cr
6: 5 + low-level debug messages}
\item{DriverManager.LOG_LEVEL_DEFAULT}{\code{4}}
\item{DriverManager.LOG_STACK_TRACE}{Whether error messages should print out exception stack traces, if any. }
\item{DriverManager.LOG_STACK_TRACE_DEFAULT}{\code{true}}
\item{DriverManager.LOG_TIMEZONE}{The log timezone.\cr
Set this to \code{"GMT"} to use Greenwich Mean Time in the timestamps on log messages. Setting this property to any other value will cause the logger to use the current timezone when generating timestamps.
default: \code{"default"}}
\item{DriverManager.LOG_TIMEZONE_DEFAULT}{\code{"default"}}
\item{DriverManager.LOG_TO_FILE}{Whether to log to file.\cr
\code{"true} or \code{"false}.}
\item{DriverManager.LOG_TO_FILE_DEFAULT}{\code{true}}
\item{DriverManager.LOG_TO_SYSTEM}{     
Whether to log to console.\cr
\code{"true} or \code{"false}.}
\item{DriverManager.MAX_CONCURRENT_COLLECT}{    
The maximum number of threads that can be downloading output data concurrently.\cr
If output data is large, you may need to decerease this value to prevent the Driver from running out of memory.
}
\item{DriverManager.MAX_CONCURRENT_COLLECT_DEFAULT}{\code{25}}
\item{DriverManager.MAX_CONCURRENT_JOBS}{The maximum number of concurrent active Jobs/Services that can be in progress at the same time.\cr
Additional Jobs started will be queued. Note that once submission for a deferred-collection service has been closed, it is not considered active on this Driver.}
\item{DriverManager.MAX_CONCURRENT_JOBS_DEFAULT}{\code{1000}}
\item{DriverManager.NEGOTIATE_ENABLED}{Set to \code{"true"} to enable Negotiate Single sign-on.\cr
See the documentation for configuration instructions. If \code{"false"}, username and password are used. default: \code{"false"}
}
\item{DriverManager.PASSWORD}{The Driver password, if using authentication.\cr
This property can only be changed before Driver initialization.
}
\item{DriverManager.PRIMARY_DIRECTOR}{The primary Director, in the form [address]:[port].\cr
If set after the Driver has initialized, it will log off from the current Manager and relogin to this Director.
}
\item{DriverManager.SECONDARY_DIRECTOR}{        
The secondary Director, in the form [address]:[port].\cr
If set after the Driver has initialized, it will log off from the current Manager and relogin to this Director.
}
\item{DriverManager.SSL_TRUST_FILE}{The location of the keystore which contains additional trusted certs.\cr
The keyword \code{auto\%} indicates that the Driver should automatically search for the keystore, first in the directory in which the \code{driver.properties} file is located. If not found there, the working directory will be checked. Default: \code{auto\%ssl.pem} }
\item{DriverManager.SSL_VERIFY_HOSTNAME}{Whether the server hostname is verified to match the CN in its certificate, when initiating an HTTPS connection.\cr
If \code{false}, the default, HostnameVerifier on \code{javax.net.ssl.HttpsURLConnection} is set to a verifier that always allows mismatches. If set to \code{true}, mismatch is not allowed and Driver will not be able to login.\cr
default: \code{false} }
\item{DriverManager.USERNAME}{The username associated with this Driver.\cr
If using client authentication, it must be a user associated with a Driver Profile.\cr
If set after the Driver has initialized, the Driver will disconnect from the Broker, and reconnect using this new username. In this way, setting this value at runtime allows the Driver to dynamically change Brokers. The password must be set prior to setting the username, if it is different and authentication is enabled..}
}
\examples{
	DriverManager.connect()
	DriverManager.disconnect()
}
\seealso{
\code{\link{Service}},
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
