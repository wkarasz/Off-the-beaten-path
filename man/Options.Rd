\name{Options}
\alias{Options}
\alias{Options.AFFINITY_DEPTH}
\alias{Options.AFFINITY_WAIT}
\alias{Options.ALLOWED_DRIVER_PROFILES}
\alias{Options.ALLOWED_ROLES}
\alias{Options.AUTO_CANCEL_ALLOWED_FAILURES}
\alias{Options.AUTO_CANCEL}
\alias{Options.AUTO_PACK_MODE}
\alias{Options.AUTO_PACK_NUM}
\alias{Options.AutoCancel.ALWAYS}
\alias{Options.AutoCancel.LIBRARY_LOAD}
\alias{Options.AutoCancel.NEVER}
\alias{Options.AutoCancel}
\alias{Options.AutoPack.PARALLEL}
\alias{Options.AutoPack.SERIAL}
\alias{Options.AutoPack}
\alias{Options.BLACKLIST_MACHINE}
\alias{Options.CHECKPOINT}
\alias{Options.Collection.AFTER_SUBMIT}
\alias{Options.Collection.EXISTING}
\alias{Options.Collection.IMMEDIATELY}
\alias{Options.Collection.LATER}
\alias{Options.Collection.NEVER}
\alias{Options.COLLECTION_TYPE}
\alias{Options.Collection}
\alias{Options.COMPRESS_DATA}
\alias{Options.EMAIL}
\alias{Options.ENCRYPTION_ENABLED}
\alias{Options.ENGINE_BLACKLISTING}
\alias{Options.ENGINE_GREYLISTING}
\alias{Options.FAILOVER_SHARED_SERVICE_MIN_TASK_ID}
\alias{Options.FAILURES_BEFORE_BLACKLIST}
\alias{Options.GRID_LIBRARY_STRICT_VERSIONING}
\alias{Options.GRID_LIBRARY_VERSION}
\alias{Options.GRID_LIBRARY}
\alias{Options.INVOCATION_MAX_TIME}
\alias{Options.INVOCATIONS_PER_MESSAGE}
\alias{Options.KILL_ENGINE_ON_CANCEL}
\alias{Options.MAX_ENGINES}
\alias{Options.MAX_TASK_RESCHEDULES}
\alias{Options.MAX_TASK_RETRIES}
\alias{Options.PENDING_TASK_TIMEOUT}
\alias{Options.Priority.HIGH}
\alias{Options.Priority.LOW}
\alias{Options.Priority.MEDIUM}
\alias{Options.Priority.URGENT}
\alias{Options.Priority}
\alias{Options.PRIORITY}
\alias{Options.PURGE_INVOCATION_DATA}
\alias{Options.PurgeInvocationData.INVOCATION_COMPLETED_SUCCESSFULLY}
\alias{Options.PurgeInvocationData.INVOCATION_COMPLETED}
\alias{Options.PurgeInvocationData.MANUAL}
\alias{Options.PurgeInvocationData.SERVICE_COMPLETED}
\alias{Options.PurgeInvocationData.SERVICE_REMOVED}
\alias{Options.PurgeInvocationData}
\alias{Options.REDUNDANT_RESCHEDULING_ENABLED}
\alias{Options.RESCHEDULE_ON_TIMEOUT}
\alias{Options.RESUBMIT_ON_DDT_FAILURE}
\alias{Options.RESULTS_PER_MESSAGE}
\alias{Options.SERVICE_FAIL_RESCHEDULE}
\alias{Options.SERVICE_FAIL_RESTART}
\alias{Options.SERVICE_FAIL_RETRY}
\alias{Options.SHARED_SERVICE_INACTIVITY_KEEP_ALIVE_TIME}
\alias{Options.SHARED_SERVICE_NAME}
\alias{Options.SHARED_UNIX_DIR}
\alias{Options.SHARED_WIN_DIR}
\alias{Options.STATE_AFFINITY}
\alias{Options.STATUS_EXPIRES}
\alias{Options.TASK_EXPIRATION_START_OFFSET}
\alias{Options.UNLOAD_NATIVE_LIBRARY}
\docType{data}
\title{
Service Options.
}
\description{
Service Options.\cr
If an option is set in the registry, it cannot be overridden on a Service Instance by the Client.\cr
If it is left as [not set], and it is not set on an instance, the default value is used.\cr
Use \code{"true"} and \code{"false"} for value with \code{boolean} semantic.
}
\value{

\item{Options.AutoCancel}{Values used with Options.AUTO_CANCEL}
\item{Options.AutoCancel.ALWAYS }{Auto-cancel the Service if any invocation fails.}
\item{Options.AutoCancel.LIBRARY_LOAD }{Auto-cancel the Service if an Engine fails loading a class or library. This is the default behavior.
Typically occurs when libraries are not deployed or up to date.
This also cancels the Service if a Task is retried and the last failure is a library load failure.}
\item{Options.AutoCancel.NEVER }{Never auto-cancel on an invocation failure.}

\item{options.Autopack}{Values used with Options.AUTO_PACK_MODE}
\item{Options.AutoPack.PARALLEL }{Invocations are executed in parallel}
\item{Options.AutoPack.SERIAL }{Invocations are executed serially}

\item{Options.Collection}{Values used with \code{Options.COLLECTION_TYPE}}
\item{Options.Collection.AFTER_SUBMIT }{
Collection of the results does not begin until either \code{waitUntilInactive(long)} or
\code{destroyWhenInactive()} is called.
This mode would typically be used when all of the work is submitted upon Service creation,
so that results collection does not slow down the speed of submission. This is typically used
in conjunction with the setting \code{Options.INVOCATIONS_PER_MESSAGE > 1}, which also speeds
submission time.
\cr 
Calling \code{Service#execute(String methodname, Object data)} will not work with this collection type,
}
\item{Options.Collection.EXISTING }{ }
\item{Options.Collection.IMMEDIATELY }{Normal operation. Results of requests are collected as soon
as they are ready.}
\item{Options.Collection.LATER }{Results are not collected in this Service; another Service
retrieved using \code{ServiceFactory.getService(...)}
will collect them later.\cr When this Service Session
type is created, it must call \code{destroy} when
finished submitting, to indicate that all requests have been
submitted, and collection may begin. If it is not called,
results cannot be collected, as it is assumed that there was
a problem during submission.\cr If the application
running the Driver is shut down, and \code{destroy}
was not called, the Session will be automatically destroyed
when the Driver connection is closed after the connection
timeout.  This is because the lack of a call indicates a
submission failure.  Otherwise, it is left open for later
collection. \code{execute} is not allowed on
this Service, since it is not collecting results.  The
\code{ServiceInvocationHandler} must be
\code{null} on all \code{submit} calls.\cr If
the Driver is being used in a mode in which is submits its
requests and then exits, either DDT should not be used, or
an external fileserver should be used by the Driver so that
Engines can retrieve submissions after the Driver is shut
down.
\cr
The \code{LATER} mode cannot be used with a Service
proxy.  Also, this mode should only be used with
short-running Tasks. If used for Services with indefinite
duration, there will be no way to clean up the inputs.}
\item{Options.Collection.NEVER }{}

\item{Options.PurgeInvocationData}{Values used with \code{Options.PURGE_INVOCATION_DATA}.
Used to determine when the data in the Task Admin for this Service is deleted from the Broker.}
\item{Options.PurgeInvocationData.INVOCATION_COMPLETED }{Purge as soon as the invocation completes.\cr In this case, the data is only available when the request is queued or being processed.}
\item{Options.PurgeInvocationData.INVOCATION_COMPLETED_SUCCESSFULLY }{Purge as soon as the invocation completes successfully. This allows you to view data for an invocation that has fails, at the
cost of memory overhead on the Broker. The data is removed when the instance is removed either automatically according to
the Service Cleanup removed either automatically according to the Service Cleanup settings on the Broker, or when it is
manually removed.}
\item{Options.PurgeInvocationData.MANUAL }{}
\item{Options.PurgeInvocationData.SERVICE_COMPLETED }{Purge all when the Service instance completes.\cr
As long as the instance is running, the data is available.}
\item{Options.PurgeInvocationData.SERVICE_REMOVED }{Purge when the completed Service instance is removed from the Administration Tool.
Allows you to view data after an instance completes, at the cost of memory overhead on the Broker.
The instance is removed either automatically according to the Service Cleanup
settings on the Broker, or when it is manually removed.}

\item{Options.Priority}{Values used with \code{Options.PRIORITY}}
\item{Options.Priority.HIGH }{Indicates that the Service should be run at the highest priority, preempting other Services if necessary}
\item{Options.Priority.LOW }{0, which means the Service is suspended and no Tasks will be taken until the priority is increased.}
\item{Options.Priority.MEDIUM }{Indicates that the Service should be run at the priority midway between HIGH and LOW.}
\item{Options.Priority.URGENT }{}

\item{Options.AFFINITY_DEPTH }{Used for Invocation-level affinity to determine how deep into the request queue the affinity score should be calculated between all available engines.\cr
Must be > 0 if adding affinity to Tasks. Larger values can result in longer scheduling episodes, so this number should be chosen wisely.\cr
Default: \code{0}
}
\item{Options.AFFINITY_WAIT }{The amount of time, in msecs, for the scheduler to wait to assign a Task to an Engine that has already worked on this \code{Service}, before assigning it to any Engine.\cr
The count starts on each task when the scheduler first attempts to assign the task to an Engine that otherwise would have taken it.\cr
Default: \code{0}
}

\item{Options.ALLOWED_DRIVER_PROFILES }{Comma-delimited list of allowed Driver Profiles that may use this \code{Service}.\cr
All Drivers are allowed if empty.}
\item{Options.ALLOWED_ROLES }{Comma-delimited list of allowed Roles that may use this \code{Service}.\cr
All Drivers are allowed if empty.}

\item{Options.AUTO_PACK_MODE }{ 
Deprecated:\cr
\code{AutoPack.PARALLEL} is no longer supported. Use Multiplexed Engines to run tasks in parallel in the same process.
}
\item{Options.AUTO_PACK_NUM }{If set to greater than 0, the number of requests processed per invocation.\cr In this mode, a \code{Service} will process multiple requests in one routine.\cr
This is accomplished by packing input data into a single request, and the request invokes the requests on all inputs. This mode should be used if there will be many more inputs than Engines and/or invocations are of short duration, to maximize efficient use of memory and Engine processing power.\cr

The Service instance will push out the data when this many inputs have been submitted, or if it has been one second since a submission, if less than this amount have been submitted.\cr

Packs will not span multiple tasks separated by updateState calls. For example, if AUTO_PACK_NUM is 5 and the client code submits 2 tasks, then 2 updateState calls, and then 3 more tasks, there will be two packs sent to the server: one with the 2 tasks submitted before the updateState, and one with the other 3. In other words, all the tasks in any one pack are all associated with a single "state".\cr

The invocations IDs that Administration Tool will report will be the IDs of the packages, so they will not directly correspond to the invocation ID from the Driver and Engine point of view.\cr

If any request fails on an Engine, the entire pack will be retried. If a request is submitted with a Scheduling or Dependency condition, it will not be packed.\cr

Default: \code{0}
}

\item{Options.BLACKLIST_MACHINE }{If \code{"true"}, when an Engine is blacklisted, all other Engines running on that machine are also blacklisted.\cr
Only meaningful if \code{ENGINE_BLACKLISTING} is enabled. Note that this does not make the number of failures cumulative across all Engines. For example, if \code{FAILURES_BEFORE_BLACKLIST} is set to \code{4}, an Engine Daemon is not blacklisted if \code{4} Engines have a single task failure, but would be blacklisted if one Engine had a task fail \code{4} times.\cr
Default: \code{false}
}

\item{Options.AUTO_CANCEL }{If the Service is automatically cancelled on a single Invocation failure.\cr See \code{Options.AutoCancel}}
\item{Options.AUTO_CANCEL_ALLOWED_FAILURES }{The amount of allowed Task failures before auto-cancellation\cr If \code{AUTO_CANCEL} is \code{Options.AutoCancel.LIBRARY_LOAD} only library load failures are counted.\cr If \code{AUTO_CANCEL} is \code{Options.AutoCancel.ALWAYS} any failure is counted.\cr If \code{AUTO_CANCEL} is \code{Options.AutoCancel.NEVER} this value has no meaning.}

\item{Options.CHECKPOINT }{If an Service invocation can checkpoint data, so that it can pick up where it left off on a reschedule.\cr Checkpointing can be useful for long-running invocations that may be interrupted.\cr Note that filenames used in checkpointing cannot contain commas.\cr See the documentation for more information on how to use checkpointing.}
\item{Options.COLLECTION_TYPE }{The mode of result collection.\cr
See \code{\link{Options.Collection}} for an explanation of the different types.\cr

Default: \code{Collection.IMMEDIATELY}
}

\item{Options.COMPRESS_DATA }{If the initial, input, and output data will be compressed. For data sizes > 10K per input or output, compression time is minimal and is recommended.}
\item{Options.EMAIL }{The optional email address, for notification when the Service is created and closes.}
\item{Options.ENGINE_BLACKLISTING }{If \code{"true"}, an Engine that fails in a Service invocation will be prevented from taking subsequent requests from the Service instance.}
\item{Options.ENCRYPTION_ENABLED }{If \code{"true"}, the tasklet, input, updates, and output will be encrypted.\cr
A new AES key is created for each \code{Service}, and is used for all encryption. The key is encrypted by the communication channels from Driver to broker to Engine. The highest keysize available will be used, \code{256} or \code{192}.\cr

Default: \code{false}
}

\item{Options.ENGINE_GREYLISTING }{If \code{"true"}, when an Engine fails a Task it will have affinity to the Service subtracted from it, to make it more likely that other Engines work on that Service.\cr Every time an Engine fails a Task, the value of the Greylist Affinity as configured on the Broker will be subtracted from the overall Service affinity.\cr This can also be used alongside blacklisting when \code{FAILURES_BEFORE_BLACKLIST} is nonzero, so that an Engine becomes less likely to work on this Service as is has failures, and is not allowed at all after too many.}
\item{Options.FAILOVER_SHARED_SERVICE_MIN_TASK_ID }{The min task ID of shared services on failover Broker.\cr
When submitting shared services tasks to failover Broker, if the primary Broker becomes available, the driver will fail back to the primary. To avoid task ID conflicts, this lower bound will be used on failover Broker. If not set, we will use \code{1610612735} as default.}
\item{Options.FAILURES_BEFORE_BLACKLIST }{The number of total Task failures allowed on this session before an Engine is blacklisted.\cr Only applied if \code{ENGINE_BLACKLISTING} is enabled.}
\item{Options.GRID_LIBRARY }{The Grid Library for this Service.}
\item{Options.GRID_LIBRARY_STRICT_VERSIONING }{When a Service is set to use the latest version of a Grid Library, this will guarantee that all Engines work on the version that was the latest at Service creation time.\cr Otherwise, if \code{"false"}, and a newer version of the library is deployed while the Service is running, Engines that have not yet worked on the Service will use the newer version, while Engines that worked on it prior to deployment will continue to use the older version.}
\item{Options.GRID_LIBRARY_VERSION }{The Grid Library version for this Service.}
\item{Options.INVOCATION_MAX_TIME }{If a running Service Invocation exceeds this amount of time (seconds), it will be rescheduled or retried, depending on the \code{RESCHEDULE_ON_TIMEOUT} setting.\cr Note that the rescheduler regularly polls for timed out requests, so it may take as long as this poll period to be rescheduled or retried.}
\item{Options.INVOCATIONS_PER_MESSAGE }{The number of invocations that will be sent per submission message.\cr
In most cases this value should be set in the range of \code{25}-\code{50} for any service that performs asynchronous task submission (as opposed to synchronous execution) for optimal throughput. When set, the Driver will buffer this number of invocation requests before sending to the Manager. To flush any remaining invocations after submitting, call \code{\link{Service.flush}}. They will also be flushed if \code{\link{Service.waitUntilInactive}} or \code{\link{Service.destroyWhenInactive}} is called. Otherwise they are flushed every second if the buffer is not full.}
\item{Options.KILL_ENGINE_ON_CANCEL }{If an Engine will be killed and restarted if a Service Invocation is cancelled.\cr Invocations are cancelled when cancelled from the Administration Tool, when a Service is cancelled, and when another Engine completes the call due to redundant rescheduling.\cr If this value is \code{"false"}, the \code{cancelled} method will be called rather than killing the Engine, to provide user-defined interruption of the Task and any necessary cleanup.}
\item{Options.MAX_ENGINES }{The maximum number of Engines that can be concurrently executing Service invocations.}
\item{Options.MAX_TASK_RESCHEDULES }{The maximum number of redundant reschedules allowed for any Task, if any of the rescheduler strategies are in effect.}
\item{Options.MAX_TASK_RETRIES }{The maximum number of retries allowed of any Task that fails.\cr A retry occurs if the Task failed and \code{SERVICE_FAIL_RETRY} is \code{"true"}, or if the Task exceeds the \code{INVOCATION_MAX_TIME} and \code{KILL_ON_TIMEOUT} is \code{"true"}.}
\item{Options.PENDING_TASK_TIMEOUT }{The number of minutes a Task will wait for an Engine before timing out.\cr This must be greater than 0 for a Task to time out.  If a Task times out, it will be scheduled on another Broker with a priority based on the Task ID.\cr The lower the Task ID, the higher the priority.}
\item{Options.PRIORITY }{The priority level of this Service instance.\cr This is a positive integer value less than the number of priority levels available on the Broker (by default 10), or one of the \code{Options.Priority} constants.}
\item{Options.PURGE_INVOCATION_DATA }{Determines when Service invocation status data on the Task Admin page is purged from the Administration Tool.
This data is stored in memory on the Broker, so this setting may need to be modified depending on how
the Service will be used.
\cr
\code{INVOCATION_COMPLETED} means that the data for each request is purged as soon as it has completed.
In this case, the data is only available when the request is queued or being processed. This should be used
when the Service may be a long-running open Service, as each entry requires memory.
\cr
\code{SERVICE_COMPLETED} means that the data is purged when the Service instance has completed.
As long as the instance is running, the data is available. This is the default.
\cr
\code{SERVICE_REMOVED} means that the data is removed when the completed Service instance is removed from the Administration Tool.
This allows you to view data after an instance completes, at the cost of memory overhead on the Broker.
The instance is removed either automatically according to the Service Cleanup
settings on the Broker, or when it is manually removed.}

\item{Options.REDUNDANT_RESCHEDULING_ENABLED }{If Redundant Rescheduling is enabled for this Service.
If \code{"false"}, any redundant strategies enabled on the \code{Services > Service Rescheduler}
configuration will be ignored.\cr
Such strategies are 'Average Limit', 'Standard Deviation Limit', and
'Remaining Task Percent'. Note the by default these are all disabled, so setting
this value to \code{"true"} will not enable redundant rescheduling unless a strategy is enabled.}
\item{Options.RESCHEDULE_ON_TIMEOUT }{How a Task is dealt with if it exceeds the
\code{INVOCATION_MAX_TIME}.  If \code{"true"}, the request is
rescheduled, and the current one continues.  If \code{"false"}, the
Engine running the Task is killed, and the Task is retried.}
\item{Options.RESUBMIT_ON_DDT_FAILURE }{The action the Driver will take if there is a failure downloading an output.
The most likely failure is a problem with the Engine Daemon fileserver,
so it usually does not indicate an error in the Invocation.
Thus, the default value is \code{"true"}, meaning that the invocation request will be resubmitted
with the same original data and ID.
If \code{"false"}, it will fail, with \code{ServiceInvocationException.getMessage()}
set to the error message that occured during download.}
\item{Options.RESULTS_PER_MESSAGE }{The maximum number of results that will be collected per message.\cr
Default is \code{100}}

\item{Options.SERVICE_FAIL_RESCHEDULE }{deprecated Use \code{SERVICE_FAIL_RETRY}.}
\item{Options.SERVICE_FAIL_RESTART }{If an Engine will restart on a Service failure.}
\item{Options.SERVICE_FAIL_RETRY }{If an Service request will be retried on a failure.
If \code{"true"}, it will only be retried up to the maximum numbers of times,
as set by \code{MAX_TASK_RETRIES}.
Note that if there is a problem downloading the input, or a problem
loading libraries, the request is always retried.}
\item{Options.SHARED_SERVICE_INACTIVITY_KEEP_ALIVE_TIME }{
The time in seconds that a Shared \code{Service} will be kept alive after the service become inactive.\cr

A Shared \code{Service} becomes inactive when all attached clients have destroyed their services, whether they are directly destroyed or destroyed by Driver timeout.\cr

Default: \code{0}
}
\item{Options.SHARED_SERVICE_NAME }{    
Sets this service to be a Shared \code{Service} with the given name.\cr

If set, this service is a Shared \code{Service}, meaning that more than one client can be attached to the same service instance. If an instance with this name does not yet exist on the Broker, it will be created by this \code{Service}. If it does, this service will attach to the instance. Results of requests are always returned to the client that made the request. If the instance is destroyed on the client, only that client's outstanding requests are cancelled, and the service instance on the Broker will remain active if there are other clients attached.
}
\item{Options.SHARED_UNIX_DIR }{If this value is set, Drivers and Engines will transfer data
using a shared directory rather than using their respective
fileservers.  This directory must be a mounted directory to
which all Unix Engines working on this Service have read/write
access.  It is optimally a directory local to the Driver, for
minimum network bandwidth.\cr If set and using Windows Engines,
the Windows shared directory must also be set to the equivalent of
this directory.}
\item{Options.SHARED_WIN_DIR }{If this value is set, Drivers and Engines will transfer data using a shared directory rather than using their respective fileservers.
This directory must be a mounted directory to which all Windows Engines working on this Service have read/write access.
It is optimally a directory local to the Driver, for minimum network bandwidth.\cr
Typically, the share will be Windows UNC format, such as \code{//server/data}"\cr
If set and using Unix Engines, the Unix shared directory must also be set to the equivalent of this directory.}
\item{Options.STATE_AFFINITY }{Overrides the default scheduler state affinity for this
Service. If higher than the default, this Service is more
likely to have Engines with state affinity for this Service
assigned to it than other Service.}
\item{Options.STATUS_EXPIRES }{If the status of the Service instance in the
Administration Tool expires.  If \code{"false"}, the status must be
manually removed. If \code{"true"}, it is removed according to the
Service Cleanup settings on the Broker.}
\item{Options.TASK_EXPIRATION_START_OFFSET }{The number of minutes a Task will wait before Task expiration is started.
The pending Task timeout counter does not start until the Task expiration offset has completed
Must be greater than or equal to 0, default is 0.
This applies only to spanned Services.}
\item{Options.UNLOAD_NATIVE_LIBRARY }{If the native library should be unloaded once the Service is finished.
Set the value to \code{"false"} for sharing global objects in the library.}
}
\seealso{
}
\keyword{datasets}
