#' start, access, stop or restart a cluster for parallel computation with rtop
#' 
#' Convenience function for using parallel computation with rtop. The function
#' is usually not called by the user.
#' 
#' It is usually not necessary for the user to call this function for starting
#' or accessing a cluster. This is done automatically by the different
#' rtop-functions when needed if the parameter nclus is larger than one (see
#' \code{\link{getRtopParams}}). If the user actually starts the cluster by a
#' call to this function, it will also be necessary to set the nclus parameter
#' to a value larger than one for the cluster to be used by different
#' functions.
#' 
#' Restarting the cluster might be necessary if the cluster has a problem (e.g.
#' does not return memory) or if the user wants to change to a different
#' cluster type.
#' 
#' Stopping the cluster is useful when the user does not want to continue with
#' parallel computation and wants to close down the workers.
#' 
#' @param nclus The number of workers in the cluster
#' @param ... Arguments for \code{\link[parallel]{clusterEvalQ}}; commands to
#' be evaluated for each worker, such as library-calls
#' @param action Defines the action of the function. There are three options:
#' \describe{ \item{"start"}{Starts a new cluster if necessary, reuses an
#' existing if it has already been started} \item{"restart"}{Stops the cluster
#' and starts it again. To be used in case there are difficulties with the
#' cluster, or if the user wants to change the type of the cluster} }
#' @param type The type of cluster; see \code{\link[parallel]{makeCluster}} for
#' more details.  The default of makeCluster is used if type is missing or NA.
#' @param outfile File to direct the output,
#' \code{\link[parallel]{makeCluster}} for more details.
#' @return If the function is called with action equal to "start" or "restart",
#' the result is a cluster with nclus workers. The cluster is also added to the
#' global options with the name rtopCluster \cr
#' (\code{getOption("rtopCluster")}).
#' 
#' If the function is called with action equal to "stop", the function stops
#' the cluster, sets the rtopCluster of options to NULL and returns NULL to the
#' user.
#' @author Jon Olav Skoien
#' @keywords spatial
#' @export
rtopCluster = function(nclus, ..., action = "start", type, outfile = NULL ) {
  cl = getOption("rtopCluster")
  if (length(cl) > 0 && (action == "stop" | action == "restart")) {
    parallel::stopCluster(cl)
    options(rtopCluster = NULL)
  } 
  if (length(cl) > 0 && action == "start") {
    if (length(list(...)) > 0) parallel::clusterEvalQ(cl, ...)
  } else if (action == "start" | action == "restart") {
    if (!requireNamespace("parallel")) stop("Not able to start cluster, parallel not available")
    if (missing(type) || is.null(type)) {
      cl <- parallel::makeCluster(nclus, outfile = outfile) 
    } else {
      cl <- parallel::makeCluster(nclus, type, outfile = outfile)
    }
#    doParallel::registerDoParallel(cl, nclus)
    if (length(list(...)) > 0) parallel::clusterEvalQ(cl, ...)
    options(rtopCluster = cl)    
  }
  getOption("rtopCluster")
}

