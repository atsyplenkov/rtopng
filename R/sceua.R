#' Optimisation with the Shuffle Complex Evolution method
#'
#' Calibration function which searches a parameter set which is minimizing the
#' value of an objective function
#'
#' sceua is an R-implementation of the Shuffle Complex Evolution - University
#' of Arizona (Duan et al., 1992), a global optimization method which "combines
#' the strengths of the simplex procedure of Nelder and Mead (1965) with the
#' concepts of controlled random search (Price, 1987), competetive evolusion
#' (Holland, 1975)" with the concept of complex shuffling, developed by Duan et
#' al. (1992).
#'
#' This implementation follows the Fortran implementation relatively close, but
#' adds the possibility of searching in log-space for one or more of the
#' parameters, and it uses the capability of R to pass functions as arguments,
#' making it possible to pass implicit conditions to the parameter selection.
#'
#' The objective function \code{OFUN} is a function which should give an error
#' value for each parameter set. It should never return non-numeric values such
#' as NA, NULL, or Inf. If some parameter combinations can give such values,
#' the return value should rather be a large number.
#'
#' The function works with fixed upper and lower boundaries for the parameters.
#' If the possible range of a parameter might span several orders of magnitude,
#' it might be better to search in log-space for the optimal parameter, to
#' reduce the risk of being trapped in local optima. This can be set with the
#' argument \code{plog}, which is either a single value (FALSE/TRUE) or a
#' vector for all parameters.  \code{plog = c(TRUE, FALSE, FALSE, TRUE, TRUE)}
#' means that the search for parameters 1,4 and 5 should be in log10-space,
#' whereas the search for parameters 2 and 3 are in normal space.
#'
#' Implicit boundaries can be evoked by passing a function \code{implicit} to
#' \code{sceua}. This function should give 0 when parameters are acceptable and
#' 1 if not. If, for example, the condition is that the following sum of
#' parameters four and five should be limited:
#'
#' `sum(pars[4]+pars[5]) <= 1`
#'
#' then the function will be implicit = function(pars) `(2*pars[4] + pars[5]) > 1`
#'
#' @param OFUN A function to be minimized, with first argument the vector of
#' parameters over which minimization is to take place. It should return a
#' scalar result as an indicator of the error for a certain parameter set
#' @param pars a vector with the initial guess the parameters
#' @param lower the lower boundary for the parameters
#' @param upper the upper boundary for the parameters
#' @param maxn the maximum number of function evaluations
#' @param kstop number of shuffling loops in which the criterion value must
#' change by the given percentage before optimization is terminated
#' @param pcento percentage by which the criterion value must change in given
#' number (kstop) of shuffling loops to continue optimization
#' @param ngs number of complexes in the initial population
#' @param npg number of points in each complex
#' @param nps number of points in a sub-complex
#' @param nspl number of evolution steps allowed for each complex before
#' complex shuffling
#' @param mings minimum number of complexes required, if the number of
#' complexes is allowed to reduce as the optimization proceeds
#' @param iniflg flag on whether to include the initial point in population.
#' iniflg <- 0, not included. iniflg= 1, included
#' @param iprint flag for controlling print-out after each shuffling loop.
#' iprint < 0: no output. iprint = 1: print information on the best point of
#' the population. iprint > 0: print information on every point of the
#' population
#' @param iround number of significant digits in print-out
#' @param peps convergence level for parameter set (lower number means smaller
#' difference between parameters of the population required for stop)
#' @param plog whether optimization should be done in log10-domain. Either a
#' single TRUE value for all parameters, or a vector with TRUE/FALSE for the
#' different parameters
#' @param implicit function for implicit boundaries for the parameters (e.g.
#' `sum(pars[4]+pars[5]) < 1`). See below for details
#' @param timeout if different from NULL: maximum time in seconds for execution
#' before the optimization returns with the parameters so far.
#' @param ... arguments for the objective function, must be named
#' @return The function returns a list with the following elements \describe{
#' \item{par}{ - a vector of the best parameters combination } \item{value}{ -
#' the value of the objective function for this parameter set}
#' \item{convergence}{ - a list of two values \describe{ \item{funConvergence}{-
#' the function convergence relative to pcento} \item{parConvergence}{ - the
#' parameter convergence relative to peps} }} \item{counts}{ - the number of
#' function evaluations} \item{iterations}{ - the number of shuffling loops}
#' \item{timeout}{ - logical; TRUE if the optimization was aborted because the
#' timeout time was reached, FALSE otherwise} } There are also two elements
#' returned as attributes: \describe{ \item{parset}{ - the entire set of
#' parameters from the last evolution step } \item{xf}{ - the values of the
#' objective function from the last evolution step } } The last two can be
#' accessed as \code{attr(sceuares, "parset")} and \code{attr(sceuares, "xf")},
#' if the result is stored as \code{sceuares}.
#' @author Jon Olav Skoien
#' @references Duan, Q., Sorooshian, S., and Gupta, V.K., 1992. Effective and
#' efficient global optimization for conceptual rainfall-runoff models. Water
#' Resour. Res. 28 (4), 1015?1031.
#'
#' Holland, H.H., 1975. Adaptation in natural and artificial systems,
#' University of Michigan Press, Ann Arbor.
#'
#' Nelder, J.A. and Mead, R., 1965. A simplex method for function minimization,
#' Comput. J., 7(4), 308-313.
#'
#' Price, W.L., 1987. Global optimization algorithms for a CAD workstation, J.
#' Optim. Theory Appl., 55(1), 133-146.
#'
#' Skoien, J. O., Bloschl, G., Laaha, G., Pebesma, E., Parajka, J., Viglione,
#' A., 2014. Rtop: An R package for interpolation of data with a variable
#' spatial support, with an example from river networks. Computers &
#' Geosciences, 67.
#' @examples
#'
#' set.seed(1)
#' # generate example data from a function with three parameters
#' # with some random noise
#' fun <- function(x, pars) pars[2]*sin(x*pars[1])+pars[3]
#' x <- rnorm(50, sd = 3)
#' y <- fun(x, pars = c(5, 2, 3)) +  rnorm(length(x), sd = 0.3)
#' plot(x,y)
#'
#' # Objective function, summing up squared differences
#' OFUN <- function(pars, x, yobs) {
#'   yvals <- fun(x, pars)
#'   sum((yvals-yobs)^2)
#' }
#'
#' sceuares <- sceua(OFUN, pars = c(0.1,0.1,0.1), lower = c(-10,0,-10),
#'                  upper = c(10,10,10), x = x, yobs = y)
#' sceuares
#' xx <- seq(min(x), max(x), 0.1)
#' lines(xx, fun(xx, pars = sceuares$par))
#'
#'
#' @export
sceua <- function(
  OFUN,
  pars,
  lower,
  upper,
  maxn = 10000,
  kstop = 5,
  pcento = 0.01,
  ngs = 5,
  npg = 2 * length(pars) + 1,
  nps = length(pars) + 1,
  nspl = 2 * length(pars) + 1,
  mings = ngs,
  iniflg = 1,
  iprint = 0,
  iround = 3,
  peps = 0.0001,
  plog = rep(FALSE, length(pars)),
  implicit = NULL,
  timeout = NULL,
  ...
) {
  # OFUN - objective function
  # pars - starting values
  # lower - lower bounds
  # upper - upper bounds
  # maxn - maximum number of iterations
  # kstop - number of shuffling loops in which the criterion value must change
  #         by the given percentage before optimization is terminated
  # pcento - percentage by which the criterion value must change in given number of shuffling loops
  # ngs - number of complexes in the initial population
  # npg = number of points in each complex
  # npt = total number of points in initial population (npt=ngs*npg)
  # nps = number of points in a sub-complex
  # nspl = number of evolution steps allowed for each complex before
  #     complex shuffling
  # mings = minimum number of complexes required, if the number of
  #     complexes is allowed to reduce as the optimization proceeds
  # iniflg = flag on whether to include the initial point in population
  #     = 0, not included
  #     = 1, included
  # iprint = flag for controlling print-out after each shuffling loop
  #     = 0, print information on the best point of the population
  #      = 1, print information on every point of the population
  # implicit = function for implicit boundaries (e.g. sum(par[4]+par[5]) < 1)

  oofun <- function(pars) OFUN(pars, ...)

  if (!is.null(timeout)) {
    tstart <- Sys.time()
  }
  npars <- length(pars)
  if (length(plog) == 1) {
    plog <- rep(plog, npars)
  }
  if (
    length(upper) != npars || length(lower) != npars || length(plog) != npars
  ) {
    stop(
      "pars, upper, lower and plog must be of same length, plog can alternatively be of length 1"
    )
  }
  pars <- ifelse(plog, log10(pars), pars)
  upper <- ifelse(plog, log10(upper), upper)
  lower <- ifelse(plog, log10(lower), lower)

  nloop <- 0
  npt <- ngs * npg
  loop <- 0
  bound <- upper - lower
  criter <- rep(1e10, 20)
  parset <- matrix(nrow = npt, ncol = npars)
  xf <- rep(1e10, npt)
  icall <- 1

  lpars <- ifelse(plog, 10^pars, pars)
  fa <- oofun(lpars)
  if (iprint > 0 && icall %% iprint == 0) {
    cat(icall, signif(fa, iround), "\n")
  }
  parset[1, ] <- pars
  xf[1] <- fa
  stdinit <- rep(1, npars)
  for (ii in ifelse(iniflg == 1, 2, 1):npt) {
    parset[ii, ] <- getpnt(idist = 1, lower, upper, stdinit, lower, implicit)
    lpars <- ifelse(
      plog,
      10^parset[ii, , drop = FALSE],
      parset[ii, , drop = FALSE]
    )
    xf[ii] <- oofun(lpars)
    icall <- icall + 1
    if (iprint > 0 && icall %% iprint == 0) {
      cat(icall, round(xf[ii], iround), "\n")
    }
  }

  parset <- parset[order(xf), , drop = FALSE]
  xf <- sort(xf)
  bestpar <- parset[1, , drop = FALSE]
  worstpar <- parset[npt, , drop = FALSE]
  bestf <- xf[1]
  worstf <- xf[npt]

  parsttout <- parstt(npt, npars, parset, bound, peps)
  ipcnvg <- parsttout$ipcnvg
  gnrng <- parsttout$gnrng
  parstd <- parsttout$parstd

  repeat {
    nloop <- nloop + 1
    for (igs in 1:ngs) {
      karr <- (c(1:npg) - 1) * ngs + igs
      cx <- parset[karr, , drop = FALSE]
      cf <- xf[karr]
      for (loop in 1:nspl) {
        kpos <- 1
        lcs <- 1
        if (nps == npg) {
          lcs <- c(1:nps)
        } else {
          repeat {
            lpos <- 1 +
              floor(
                npg + 0.5 - sqrt((npg + 0.5)^2 - npg * (npg + 1) * runif(1))
              )
            lcs[kpos] <- lpos
            if (sum(duplicated(lcs)) == 0) {
              kpos <- kpos + 1
            }
            if (kpos > nps) break
          }
        }
        lcs <- sort(lcs)

        soc <- cx[lcs, , drop = FALSE]
        sf <- cf[lcs]
        cceout <- cce(
          oofun,
          npars,
          nps = nps,
          soc = soc,
          sf = sf,
          lower = lower,
          upper = upper,
          parstd = parstd,
          icall = icall,
          maxn = maxn,
          iprint = iprint,
          iround = iround,
          bestf = bestf,
          plog = plog,
          implicit = implicit
        )
        soc <- cceout$soc
        sf <- cceout$sf
        icall <- cceout$icall
        cx[lcs, ] <- soc
        cf[lcs] <- sf
        cx <- cx[order(cf), , drop = FALSE]
        cf <- sort(cf)
        if (!is.null(timeout)) {
          if (difftime(Sys.time() - tstart, "secs") > timeout) {
            return(list(
              par = ifelse(plog, 10^bestpar, bestpar),
              value = xf[1],
              convergence = list(
                funConvergence = signif(concrit, iround) / pcento,
                parConvergence = gnrng / peps
              ),
              counts = icall,
              iterations = nloop,
              timeout = TRUE
            ))
          }
        }
      }
      parset[karr, ] <- cx
      xf[karr] <- cf
    }
    parset <- parset[order(xf), , drop = FALSE]
    xf <- sort(xf)
    bestpar <- parset[1, , drop = FALSE]
    worstpar <- parset[npt, , drop = FALSE]
    bestf <- xf[1]
    worstf <- xf[npt]
    parsttout <- parstt(npt, npars, parset, bound, peps)
    ipcnvg <- parsttout$ipcnvg
    gnrng <- parsttout$gnrng
    parstd <- parsttout$parstd
    fbestf <- criter[kstop]
    concrit <- 2 * (fbestf - bestf) / (fbestf + bestf)
    criter[2:length(criter)] <- criter[1:(length(criter) - 1)]
    criter[1] <- bestf
    if (iprint >= 0) {
      cat(
        icall,
        "best",
        signif(bestf, iround),
        "function convergence",
        signif(concrit, iround) / pcento,
        "parameter convergence",
        gnrng / peps,
        "\n"
      )
    }

    if (concrit < pcento && ipcnvg == 1) {
      break
    }
    if (icall > maxn) {
      break
    }
    if (ngs > mings) {
      compout <- comp(npars, npt, ngs, npg, parset, xf)
      ngs <- ngs - 1
      parset <- compout$parset
      xf <- compout$xf
    }
  }
  bestpar <- ifelse(plog, 10^bestpar, bestpar)
  retList <- list(
    par = bestpar,
    value = xf[1],
    convergence = list(
      funConvergence = signif(concrit, iround) / pcento,
      parConvergence = gnrng / peps
    ),
    counts = icall,
    iterations = nloop,
    timeout = FALSE
  )
  parset <- ifelse(plog, 10^parset, parset)
  attr(retList, "parset") <- parset
  attr(retList, "xf") <- xf
  retList
}

#' @noRd
comp <- function(npars, npt, ngs, npg, parset, xf) {
  xn <- parset
  xfn <- xf
  for (igs in 1:ngs) {
    karr1 <- (c(1:npg) - 1) * ngs + igs
    karr2 <- (c(1:npg) - 1) * (ngs - 1) + igs
    xn[karr2, ] <- parset[karr1, , drop = FALSE]
    xfn[karr2] <- xf[karr1]
  }
  return(list(parset = xn, xf = xfn))
}

#' @noRd
cce <- function(
  oofun,
  npars,
  nps,
  soc,
  sf,
  lower,
  upper,
  parstd,
  icall,
  maxn,
  iprint,
  iround,
  bestf,
  plog,
  implicit
) {
  alpha <- 1.
  beta <- 0.5
  n <- dim(soc)[1]
  sb <- soc[1, , drop = FALSE]
  sw <- soc[n, , drop = FALSE]
  ce <- colMeans(soc)
  fw <- sf[n]
  snew <- ce + alpha * (ce - sw)
  #  print(icall)
  if (chkcst(snew, lower, upper, implicit) > 0) {
    snew <- getpnt(2, lower, upper, parstd, sb, implicit)
  }
  #  print(snew)
  lpars <- ifelse(plog, 10^snew, snew)
  fnew <- oofun(lpars)
  icall <- icall + 1
  if (iprint > 0 && icall %% iprint == 0) {
    cat(icall, lpars, signif(fnew, iround), signif(bestf, iround), "\n")
  }
  if (fnew > fw) {
    snew <- ce - beta * (ce - sw)
    lpars <- ifelse(plog, 10^snew, snew)
    fnew <- oofun(lpars)
    icall <- icall + 1
    if (iprint > 0 && icall %% iprint == 0) {
      cat(icall, signif(fnew, iround), signif(bestf, iround), "\n")
    }
    if (fnew > fw) {
      snew <- getpnt(2, lower, upper, parstd, sb, implicit)
      lpars <- ifelse(plog, 10^snew, snew)
      fnew <- oofun(lpars)
      icall <- icall + 1
      if (iprint > 0 && icall %% iprint == 0) {
        cat(icall, signif(fnew, iround), signif(bestf, iround), "\n")
      }
    }
  }
  soc[n, ] <- snew
  sf[n] <- fnew
  return(list(soc = soc, sf = sf, icall = icall))
}

#' @noRd
chkcst <- function(parlocal, lower, upper, implicit) {
  ibound <- if (
    sum(mapply(
      FUN = function(x, y, z) {
        max(y - x, x - z, 0)
      },
      parlocal,
      lower,
      upper
    )) >
      0
  ) {
    1
  } else {
    0
  }
  if (ibound == 0 && length(parlocal) > 1 && !is.null(implicit)) {
    # Possibility to include implicit constraints
    if (!is.function(implicit)) {
      stop("implicit has to be a function")
    }
    ibound <- implicit(parlocal)
  }
  ibound
}


#' @noRd
getpnt <- function(idist, lower, upper, std, pari, implicit) {
  #  rand = (ifelse(rep(idist,npars) == 1,runif(npars),rnorm(npars)))
  #  print(xi)
  #  print(rand)
  ic <- 0
  repeat {
    parj <- mapply(
      FUN = get1p,
      pari,
      std = std,
      lower = lower,
      upper = upper,
      MoreArgs = list(idist = idist, implicit = implicit)
    )
    if (chkcst(parj, lower, upper, implicit) == 0) {
      break
    }
    ic <- ic + 1
    if (ic > 100) {
      stop(
        "Cannot find a parameter set respecting the fixed or implicit boundaries after 100 iterations"
      )
    }
  }
  return(parj)
}

#' @noRd
get1p <- function(pari, std, lower, upper, idist, implicit) {
  #  print(paste(xi,std,rand,lower,upper))
  ic <- 0
  repeat {
    rand <- ifelse(idist == 1, runif(1), rnorm(1))
    parj <- pari + std * rand * (upper - lower)
    #    print(x)
    #    print(chkcst(x,lower,upper))
    if (chkcst(parj, lower, upper, implicit) == 0) {
      break
    }
    ic <- ic + 1
    if (ic > 100) {
      stop(
        "Not possible to find a parameter that respect the fixed or implicit boundaries after 100 iterations"
      )
    }
    #    print(acdf)
  }
  return(parj)
}

#' @noRd
parstt <- function(npt, npars, parset, bound, peps) {
  parstd <- apply(parset, MARGIN = 2, FUN = function(x) sd(x)) / bound
  parmin <- apply(parset, MARGIN = 2, FUN = function(x) min(x))
  parmax <- apply(parset, MARGIN = 2, FUN = function(x) max(x))
  gsum <- sum(log((parmax - parmin) / bound))
  gnrng <- exp(gsum / npars)
  ipcnvg <- ifelse(gnrng <= peps, 1, 0)
  return(list(ipcnvg = ipcnvg, gnrng = gnrng, parstd = parstd))
}

#FUN = function(pars,target) (pars[1]*pars[2]*pars[3]-target)^2

#p0 = c(1,1,1)
#upper = c(20,20,20)
#lower = c(-20,-20,-20)
#best = sceua(FUN,p0,bl,bu,plog=FALSE,target = 3.75)
# plog is a logical to define if parameters should be logarithmized or not
