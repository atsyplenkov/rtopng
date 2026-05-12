# Remarks
# Is it reasonable to have ainfo <<- read.area.info(finfo,...) to make sure that ainfo
# is also available at the top level after being delivered to function read.areas?

#' create SpatialPointsDataFrame with observations of data with a spatial
#' support
#'
#' readAreaInfo will read a text file with observations and descriptions of
#' data with a spatial support.
#'
#' The function is of particular use when data are not available as
#' shape-files, or when the observations are not part of the shape-files. This
#' function is mainly for compatibility with the former FORTRAN-version. The
#' simplest way to read the data in that case is through
#' \code{\link[sf]{st_read}}. See also \code{\link{rtopng-package}}.
#'
#' @param fname name of file with areal information
#' @param id name of column with observation id
#' @param iobs name of column with number of observations
#' @param obs name of column with observations
#' @param unc name of column with possible uncertainty of observation
#' @param filenames name of column with filenames of areas if different names
#' than id should be used.
#' @param sep separator in csv-file
#' @param debug.level used for giving additional output
#' @param moreCols name of other column names the user wants included in ainfo
#' @return \code{\link[sp:SpatialPoints]{SpatialPointDataFrame}} with
#' information about observations and/or predictionLocations.
#' @author Jon Olav Skoien
#' @keywords spatial
#' @export
readAreaInfo <- function(
  fname = "ainfo.txt",
  id = "id",
  iobs = "iobs",
  obs = "obs",
  unc = "unc",
  filenames = "filenames",
  sep = "\t",
  debug.level = 1,
  moreCols = list(NULL)
) {
  # Separate function to read in information about the areas, with possibility to define column names
  # fname = name of file with information
  # inum = area number
  # id = internal identity - the way areas are stored on hard drive
  # iobs = number of observations for each station
  # obs = the actual observation
  # unc = The uncertainty of the observation, standard deviation
  #       This variable is optional
  # filenames = filenames for areas
  # MoreCols = other variables the user wants to pass on to ainfo
  #  cat(paste(fname))
  if (debug.level > 1) {
    print(paste(fname, id, iobs, obs, unc, filenames, debug.level, moreCols))
  }
  ainfot <- read.csv(fname, header = TRUE, sep = sep)
  if (debug.level > 1) {
    print(summary(ainfot))
  }
  ainfo <- data.frame(
    id = ainfot[, names(ainfot) == id],
    iobs = ainfot[, names(ainfot) == iobs],
    obs = ainfot[, names(ainfot) == obs]
  )
  if (unc %in% names(ainfot)) {
    ainfo <- data.frame(ainfo, unc = ainfot[, names(ainfot) == unc])
  }
  if (filenames %in% names(ainfot)) {
    ainfo <- data.frame(ainfo, filenames = ainfot[, names(ainfot) == filenames])
  }
  # Including remaining arguments, if user wants to include
  ncols <- length(moreCols)
  if (!is.null(moreCols[[1]])) {
    for (icol in 1:ncols) {
      col <- moreCols[[icol]]
      if (debug.level > 1) {
        cat(paste("moreCols", icol, col, "\n"))
      }
      ainfo <- data.frame(ainfo, col1 = ainfot[, names(ainfot) == col])
      names(ainfo) <- c(names(ainfo[1:(length(names(ainfo)) - 1)]), col)
      if (debug.level > 1) {
        cat(paste("moreCols2", icol, col, "\n"))
      }
      if (debug.level > 1) cat(paste("names(ainfo)", names(ainfo), "\n"))
    }
  }
  return(ainfo)
}


#' help file for creating SpatialPolygonsDataFrame with observations and/or
#' predictionLocations of data with a spatial support
#'
#' readAreas will read area-files, add observations and convert the result to
#' \cr \code{\link[sp]{SpatialPolygonsDataFrame}}
#'
#' If \code{object} is a file name, \code{\link{readAreaInfo}} will be called.
#' If it is a \cr \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}} with
#' observations and/or predictionLocations, the function will read areal data
#' from files according to the ID associated with each
#' observation/predictionLocation.
#'
#' The function is of particular use when data are not available as
#' shape-files, or when the observations are not part of the shape-files. This
#' function is mainly for compatibility with the former FORTRAN-version. The
#' simplest way to read the data in that case is through
#' \code{\link[sf]{st_read}}. See also \code{\link{rtopng-package}}.
#'
#' @param object either name of file with areal information or
#' \code{\link[sp]{SpatialPointsDataFrame}} with observations
#' @param adir directory where the files with areal information are to be found
#' @param ftype type of file, the only type supported currently is "xy",
#' referring to x- and y-coordinates of boundaries
#' @param projection add projection to the object if input is boundary-files
#' @param ... further parameters to be passed to \code{\link{readAreaInfo}}
#' @return The function creates a
#' \code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}} of observations
#' and/or predictionLocations, depending on the information given in
#' \code{object}.
#' @author Jon Olav Skoien
#' @keywords spatial
#' @export
readAreas <- function(object, adir = ".", ftype = "xy", projection = NA, ...) {
  # ainfo is 1 - ainfo e.g. read by readAreaInfo
  #          2 - name of the file to pass to readAreaInfo. ainfo is in that case delivered as a top level data.frame
  # pdif gives directory to areal information
  # need option to use other separators as well
  # Output of this function is a list consisting of
  #      1 Updated version of ainfo
  #      2 A list of Spatial polygons defining the borders of the areas or a set of Spatial grids
  if (is.character(object)) {
    cat(paste("calling readAreaInfo with filename ", object, "\n"))
    ainfo <- readAreaInfo(object, ...)
  } else {
    ainfo <- object
  }
  cat(paste(names(ainfo)))
  cat(paste("\n"))
  if (sum(names(ainfo) == "filenames") == 1) {
    fnames <- paste(adir, "/", ainfo$filenames, sep = "")
  } else {
    fnames <- ainfo$id
  }
  areas <- list()
  row.names(ainfo) <- c(1:dim(ainfo)[1])
  if (ftype == "xy") {
    fnames <- paste(adir, "/", fnames, ".xy", sep = "")
    Srl <- list()
    for (i in seq_along(fnames)) {
      cat(paste("reading first polygon", i, length(fnames), "\n"))
      boun <- read.table(fnames[i], header = FALSE)
      names(boun) <- c("x", "y")
      coordinates(boun) <- ~ x + y
      boun <- Polygon(boun)
      cat(paste("adding data to ainfo", i, boun@area, boun@labpt[1], "\n"))
      cat(paste(" Finished polygon\n"))
      Srl[[i]] <- Polygons(list(boun), ID = as.character(i))
    }
    Sr <- SpatialPolygons(Srl, proj4string = CRS(as.character(projection)))
    #  } else if (ftype == "shp") {
    # This part is when each file has a single shape
    # NOT PROPERLY IMPLEMENTED - need testing with real shapes
    # The files will probably be read as lists of polygons, necessary
    # to extract the actual polygon
    #    require(maptools)
    #    fnames = paste(adir,"/",fnames,".shp",sep="")
    #    Srl = list()
    #    for (i in 1:dim(fnames)) {
    #      boun = readShapePoly(fnames[i])
    #      Srl[[i]] = Polygons(boun,ID = as.character(i))
    #    }
    #    Sr = SpatialPolygons(Srl, proj4string=CRS(as.character(projection)))
    #  } else if (ftype == "shps") {
    #This clause is when one shapefile includes all the shapes
    # NOT properly tested yet
    # Necessary to split all shapes into single Polygons
    #    require(maptools)
    #    Sr = readShapePoly(adir)
    #    SPDF = SpatialPolygonsDataFrame(Sr,data = ainfo)
  } else {
    stop(paste("Filetype", ftype, "not recognized"))
  }
  ainfo$area <- unlist(lapply(Sr@polygons, FUN = function(poly) poly@area))
  ainfo$labx <- unlist(lapply(Sr@polygons, FUN = function(poly) poly@labpt[1]))
  ainfo$laby <- unlist(lapply(Sr@polygons, FUN = function(poly) poly@labpt[2]))
  ainfo$bdim <- unlist(lapply(Sr@polygons, FUN = function(poly) dim(poly)[1]))
  SPDF <- SpatialPolygonsDataFrame(Sr, data = ainfo, match.ID = TRUE)
  SPDF
}
