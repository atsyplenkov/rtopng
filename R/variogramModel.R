#' create or update variogram model
#'
#' This gives an easier interface to the parameters of the variogram model
#'
#'
#' @aliases variogramModel rtopVariogramModel updateRtopVariogram updateRtopVariogram.rtop
#' updateRtopVariogram.rtopVariogramModel
#' @param model variogram model, currently "Ex1" is the only implemented, see
#' Skoien et al (2006)
#' @param sill sill of variogram
#' @param range range of variogram
#' @param exp the exponent of the fractal part of the variogram, see Skoien et
#' al (2006)
#' @param exp0 gives the angle of the first part of the variogram in a log-log
#' plot (weibull type), should be between 0 and 2. See Skoien et al (2006)
#' @param nugget nugget of point variogram
#' @param formulaString formula that defines the dependent variable as a linear
#' model of independent variables, see e.g. \code{\link{createRtopObject}} for
#' more details.
#' @param object either: object of class \code{rtop} (see
#' \link{rtopng-package}), or an rtopVariogramModel.
#' @param action character variable defining whether the new parameters should
#' be \code{add}(-ed), \code{mult}(-iplied) or \code{replace} the former
#' parameters.  Leaving the parameters equal to NULL will cause no change.
#' @param checkVario logical, will issue a call to\code{\link{checkVario}} if
#' TRUE
#' @param sampleVariogram a sample variogram of the data
#' @param observations a set of observations
#' @param ... parameters to lower level functions
#' @return The function helps creating and updating the parameters of the
#' variogram, by using common names and simple update methods. This is mainly
#' for manual fitting of the variogram. The automatic call to checkVario makes
#' it easier to visualize the effect of the changes to the variogram
#' @author Jon Olav Skoien
#' @seealso \code{\link{rtopng-package}}
#' @keywords spatial
#' @examples
#'
#' \dontrun{
#' library(sf)
#' rpath <- system.file("extdata",package="rtopng")
#' observations <- st_read(rpath,"observations")
#' # Create a column with the specific runoff:
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#' predictionLocations <- st_read(rpath,"predictionLocations")
#' rtopObj <- createRtopObject(observations,predictionLocations)
#'  # Fit a variogram (function also creates it)
#' rtopObj <- rtopFitVariogram(rtopObj)
#' rtopObj <- updateRtopVariogram(rtopObj, exp = 1.5, action = "mult",
#'               checkVario = TRUE)
#' }
#' @export
rtopVariogramModel <- function(
  model = "Ex1",
  sill = NULL,
  range = NULL,
  exp = NULL,
  nugget = NULL,
  exp0 = NULL,
  observations = NULL,
  formulaString = obs ~ 1
) {
  if (tolower(model) == "ex1") {
    model <- "Ex1"
    if (!is.null(observations)) {
      parInit <- findParInit(formulaString, observations, model)$par0
      if (is.null(sill)) {
        sill <- parInit[1]
      }
      if (is.null(range)) {
        range <- parInit[2]
      }
      if (is.null(nugget)) {
        nugget <- parInit[3]
      }
      if (is.null(exp)) {
        exp <- parInit[4]
      }
      if (is.null(exp0)) exp0 <- parInit[5]
    } else {
      if (is.null(sill)) {
        sill <- 1
      }
      if (is.null(range)) {
        range <- 1
      }
      if (is.null(exp)) {
        exp <- 0
      }
      if (is.null(nugget)) {
        nugget <- 0
      }
      if (is.null(exp0)) exp0 <- 1
    }

    variogramModel <- list(
      model = model,
      params = c(sill, range, nugget, exp, exp0)
    )
    class(variogramModel) <- "rtopVariogramModel"
  }
  variogramModel
}

#' @rdname rtopVariogramModel
#' @export
updateRtopVariogram.rtop <- function(object, ...) {
  object$variogramModel <- updateRtopVariogram(
    object$variogramModel,
    sampleVariogram = object$variogram,
    observations = object$observations,
    ...
  )
  object
}

#' @rdname rtopVariogramModel
#' @export
updateRtopVariogram.rtopVariogramModel <- function(
  object,
  action = "mult",
  ...,
  checkVario = FALSE,
  sampleVariogram = NULL,
  observations = NULL
) {
  variogramModel <- object
  dots <- list(...)

  if (variogramModel$model == "Ex1") {
    if ("sill" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[1] <- variogramModel$params[1] * dots$sill
      } else if (action == "add") {
        variogramModel$params[1] <- variogramModel$params[1] * dots$sill
      } else if (action == "replace") {
        variogramModel$params[1] <- dots$sill
      }
    }
    if ("range" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[2] <- variogramModel$params[2] * dots$range
      } else if (action == "add") {
        variogramModel$params[2] <- variogramModel$params[2] * dots$range
      } else if (action == "replace") {
        variogramModel$params[2] <- dots$range
      }
    }
    if ("nugget" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[3] <- variogramModel$params[3] * dots$nugget
      } else if (action == "add") {
        variogramModel$params[3] <- variogramModel$params[3] * dots$nugget
      } else if (action == "replace") {
        variogramModel$params[3] <- dots$nugget
      }
    }
    if ("exp" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[4] <- variogramModel$params[4] * dots$exp
      } else if (action == "add") {
        variogramModel$params[4] <- variogramModel$params[4] * dots$exp
      } else if (action == "replace") {
        variogramModel$params[4] <- dots$exp
      }
    }
    if ("exp0" %in% names(dots)) {
      if (action == "mult") {
        variogramModel$params[5] <- variogramModel$params[5] * dots$exp0
      } else if (action == "add") {
        variogramModel$params[5] <- variogramModel$params[5] * dots$exp0
      } else if (action == "replace") {
        variogramModel$params[5] <- dots$exp0
      }
    }
    if (checkVario) {
      checkVario(
        variogramModel,
        sampleVariogram = sampleVariogram,
        observations = observations
      )
    }
  }
  variogramModel
}


#' Plot and Identify Data Pairs on Sample Variogram Cloud
#'
#' Plot a sample variogram cloud, possibly with identification of individual
#' point pairs
#'
#'
#' @param x object of class \code{variogramCloud}
#' @param ...  parameters that are passed through to
#' \code{\link[gstat]{plot.variogramCloud}} The most important are: \describe{
#' \item{identify}{ logical; if TRUE, the plot allows identification of a
#' series of individual point pairs that correspond to individual variogram
#' cloud points (use left mouse button to select; right mouse button ends) }
#' \item{digitize}{ logical; if TRUE, select point pairs by digitizing a region
#' with the mouse (left mouse button adds a point, right mouse button ends) }
#' \item{xlim}{ limits of x-axis } \item{ylim}{ limits of y-axis } \item{xlab}{
#' x axis label } \item{ylab}{ y axis label } \item{keep}{ logical; if TRUE and
#' \code{identify} is TRUE, the labels identified and their position are kept
#' and glued to object x, which is returned. Subsequent calls to plot this
#' object will now have the labels shown, e.g. to plot to hardcopy } }
#' @note This function is mainly a wrapper around
#' \code{\link[gstat]{plot.variogramCloud}}, necessary because of different
#' column names and different class names. The description of arguments and
#' value can therefore be found in the help page of
#' \code{\link[gstat]{plot.variogramCloud}}.
#' @author Jon Olav Skoien
#' @seealso \link[gstat]{plot.gstatVariogram}
#' @references \url{http://www.gstat.org/}
#' @keywords dplot
#' @examples
#'
#' \donttest{
#' rpath <- system.file("extdata",package="rtopng")
#' library(sf)
#' observations <- st_read(rpath, "observations")
#'
#' observations$obs <- observations$QSUMMER_OB/observations$AREASQKM
#'
#' # Create the sample variogram
#' rtopVario <- rtopVariogram(observations, params = list(cloud = TRUE))
#' plot(rtopVario)
#'
#' }
#'
#' @method plot rtopVariogramCloud
#' @exportS3Method plot rtopVariogramCloud
plot.rtopVariogramCloud <- function(x, ...) {
  x$np <- x$ord
  class(x) <- "variogramCloud"
  plot(x, ...)
}
