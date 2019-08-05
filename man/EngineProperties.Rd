\name{EngineProperties}
\alias{EngineProperties}
\alias{EngineProperties.AMI_AVAIL_ZONE}
\alias{EngineProperties.AMI_AVAIL_ZONE}
\alias{EngineProperties.AMI_INSTANCE_ID}
\alias{EngineProperties.AMI_PRODUCT_CODE}
\alias{EngineProperties.AVAILABLE_DISK_MB}
\alias{EngineProperties.CONFIGURATION_NAME}
\alias{EngineProperties.CPU_CORE_COUNT}
\alias{EngineProperties.CPU_ID_STRING}
\alias{EngineProperties.CPU_NO}
\alias{EngineProperties.CPU_SOCKET_COUNT}
\alias{EngineProperties.CPU_THREAD_COUNT}
\alias{EngineProperties.CRT_VERSIONS}
\alias{EngineProperties.CUDA_DEVICES}
\alias{EngineProperties.CUDA_FIRST_GPU_NAME}
\alias{EngineProperties.CUDA_GLOBAL_MEMORY}
\alias{EngineProperties.CUDA_PROCESSORS}
\alias{EngineProperties.CUDA_VERSION}
\alias{EngineProperties.CURRENT_CPU}
\alias{EngineProperties.DATA_DIR}
\alias{EngineProperties.DATA_URL}
\alias{EngineProperties.FREE_DISK_MB}
\alias{EngineProperties.FREE_MEMORY_KB}
\alias{EngineProperties.GRID_LIBRARIES}
\alias{EngineProperties.GRID_LIBRARIES_UPDATE_TIME}
\alias{EngineProperties.GRID_LIBRARY}
\alias{EngineProperties.GRID_LIBRARY_DIR}
\alias{EngineProperties.GRID_LIBRARY_VERSION}
\alias{EngineProperties.GUID}
\alias{EngineProperties.HIGHEST_NUMA_NODE}
\alias{EngineProperties.HOME_BROKERS}
\alias{EngineProperties.HOME_DIR}
\alias{EngineProperties.ID}
\alias{EngineProperties.INSTANCE}
\alias{EngineProperties.IP}
\alias{EngineProperties.MIC_CORES}
\alias{EngineProperties.MIC_DEVICES}
\alias{EngineProperties.MIC_MEMORY}
\alias{EngineProperties.MIC_MODEL}
\alias{EngineProperties.MULTICORE_ENGINE_COUNT}
\alias{EngineProperties.MULTICORE_MASTER_ID}
\alias{EngineProperties.NET_FRAMEWORK}
\alias{EngineProperties.NET_FRAMEWORK_VERSIONS}
\alias{EngineProperties.OS}
\alias{EngineProperties.OS_NAME}
\alias{EngineProperties.OS_USERNAME}
\alias{EngineProperties.OS_VERSION}
\alias{EngineProperties.PID}
\alias{EngineProperties.SCHEDULING_POLICY}
\alias{EngineProperties.SHARED_DIR_ID}
\alias{EngineProperties.TOTAL_CPU}
\alias{EngineProperties.TOTAL_MEMORY_KB}
\alias{EngineProperties.USERNAME}
\alias{EngineProperties.VALID_LIBRARIES}
\alias{EngineProperties.WORK_URL}
\docType{data}
\title{
EngineProperties used with Discriminators
}
\description{
These are the \code{key} used with Discriminators.setComparator\cr
}
\value{
\item{EngineProperties.AMI_AVAIL_ZONE }{The AMI Availability Zone.}
\item{EngineProperties.AMI_INSTANCE_ID }{The AMI Instance ID.}
\item{EngineProperties.AMI_PRODUCT_CODE }{The AMI Product Code.}
\item{EngineProperties.AVAILABLE_DISK_MB }{The amount of available disk space on the Engine's host computer, in megabytes, updated when any Engine instance is launched.\cr  The value does not include disk space reserved for superusers in UNIX.  }
\item{EngineProperties.CONFIGURATION_NAME }{The current Engine configuration name.\cr The value is [os]:[name] like on the Engine Configuration page.}
\item{EngineProperties.CPU_CORE_COUNT }{The total number of CPU cores on the Engine's host computer, if available.}
\item{EngineProperties.CPU_ID_STRING }{The CPUID string of the first physical CPU on the Engine's host computer, if available.\cr A value of -1 indicates that a value could not be determined.}
\item{EngineProperties.CPU_NO }{The number of CPUs on the Engine's host computer.}
\item{EngineProperties.CPU_SOCKET_COUNT }{The total number of physical CPUs on the Engine's host computer, if available.\cr  A value of -1 indicates that a value could not be determined.}
\item{EngineProperties.CPU_THREAD_COUNT }{The total number of hardware threads on the Engine's host computer, if available.\cr A value of -1 indicates that a value could not be determined.}
\item{EngineProperties.CRT_VERSIONS }{Installed Windows C Runtime Versions (VC71, VC8, VC9, VC10)}
\item{EngineProperties.CUDA_DEVICES }{The number of CUDA devices detected on this Engine.}
\item{EngineProperties.CUDA_FIRST_GPU_NAME }{The name of the first CUDA device.}
\item{EngineProperties.CUDA_GLOBAL_MEMORY }{The amount of CUDA global memory supported.\cr If there is more than one device, this is the minimum amount supported.}
\item{EngineProperties.CUDA_PROCESSORS }{This is the number of CUDA processors supported.\cr If there is more than one device, this is the minimum number of processors.}
\item{EngineProperties.CUDA_VERSION }{The is the CUDA capability version.\cr If there is more than one device, this is the minimum version.}
\item{EngineProperties.CURRENT_CPU }{}
\item{EngineProperties.DATA_DIR }{Data directory, the directory which is home to the DDT (direct data transfer) data.}
\item{EngineProperties.DATA_URL }{Data URL, the URL corresponding to the DDT directory.}
\item{EngineProperties.FREE_DISK_MB }{The amount of free disk space on the Engine's host computer, in megabytes, updated when any Engine instance is launched.\cr In some UNIX Systems, the value also includes disk space reserved for superusers.}
\item{EngineProperties.FREE_MEMORY_KB }{The amount of free physical memory on the Engine's host computer, in kilobytes, updated when any Engine instance is launched.}
\item{EngineProperties.GRID_LIBRARIES }{All downloaded Grid libraries, in the format "gridlibName1[,version1.1,version1.2];gridlibName2[,version2.1]...".\cr If a version does not exist, the version part may be empty.}
\item{EngineProperties.GRID_LIBRARIES_UPDATE_TIME }{The local time on the Engine when the grid libraries were last updated, in milliseconds from Unix Epoch, that is the value of Java System.currentTimeMillis().}
\item{EngineProperties.GRID_LIBRARY }{The current Grid Library.}
\item{EngineProperties.GRID_LIBRARY_DIR }{The full path to the directory of the root Grid Library for the current Session.\cr Null if using default resources.}
\item{EngineProperties.GRID_LIBRARY_VERSION }{The current Grid Library version.}
\item{EngineProperties.GUID }{The globally unique identifier for the Engine (network card address).}
\item{EngineProperties.HIGHEST_NUMA_NODE }{Highest NUMA node available in the system}
\item{EngineProperties.HOME_BROKERS }{The comma-delimited list of home Brokers.\cr This is set in the Engine Configuration.}
\item{EngineProperties.HOME_DIR }{Home directory, the installation directory on the Engine's host computer.}
\item{EngineProperties.ID }{The numeric ID unique to each Engine, generated upon installation.}
\item{EngineProperties.INSTANCE }{The instance of the Engine,  starting with zero.\cr Multi-CPU hosts may have multiple instances of Engines running concurrently. Multiplexed Engines will return 0 for all Engine instances.}
\item{EngineProperties.IP }{The IP address of the Engine.}
\item{EngineProperties.MIC_CORES }{The number of processors supported.\cr If there is more than one device, this is the number of logical processors of the oldest co-processor.}
\item{EngineProperties.MIC_DEVICES }{The number of MIC co-processors detected on this Engine.}
\item{EngineProperties.MIC_MEMORY }{The amount of memory supported.\cr If there is more than one device, this is the amount supported by the oldest co-processor.}
\item{EngineProperties.MIC_MODEL }{The identification of the oldest MIC device.}
\item{EngineProperties.MULTICORE_ENGINE_COUNT }{The number of Engines instances in a multiplexed Engine.}
\item{EngineProperties.MULTICORE_MASTER_ID }{The Engine session ID of the 0th instance of the multiplexed Engine.  (This is used internally only.)}
\item{EngineProperties.NET_FRAMEWORK }{If the .NET Framework is installed on this Engine.}
\item{EngineProperties.NET_FRAMEWORK_VERSIONS }{A comma-separated list of the supported versions of the .NET Framework available on this Engine.}
\item{EngineProperties.OS }{The operating system.\cr  Supported OS platforms can be viewed on the Engine Install page.}
\item{EngineProperties.OS_NAME }{The specific name of the operating system of the engine's host.}
\item{EngineProperties.OS_USERNAME }{The owner of the process running the Engine.}
\item{EngineProperties.OS_VERSION }{The version of the operation system.}
\item{EngineProperties.PID }{The process ID of the Engine process.}
\item{EngineProperties.SCHEDULING_POLICY }{Engine scheduling policy}
\item{EngineProperties.SHARED_DIR_ID }{A number which identifies a group of Engines running from the same home directory.}
\item{EngineProperties.TOTAL_CPU }{The amount of processing power on the Engine's host computer, measured in millions of floating-point arithmetic operations per second.\cr Calculated on the Engine Daemon using a Linpack performance calculation and reported when the Engine Daemon logs into the Manager and updated on a frequency that is set in the Engine Configuration.}
\item{EngineProperties.TOTAL_MEMORY_KB }{The amount of total physical memory on the Engine's host computer, in kilobytes.}
\item{EngineProperties.USERNAME }{The username of the Engine, which will be the hostname of the Engine unless installed with tracking properties.}
\item{EngineProperties.VALID_LIBRARIES }{}
\item{EngineProperties.WORK_URL }{Work URL, the URL of work directory on the Engine's fileserver that stores log and temporary files.}
%% ...
}
\seealso{
\code{\link{Discriminator.EComparator}}, \code{\link{Discriminator.setComparator}}
}
\examples{
}
\keyword{datasets}
