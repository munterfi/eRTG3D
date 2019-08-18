.setUpPSOCKcluster <- function(nNodes, packages, export, envir) {
  cl <- parallel::makePSOCKcluster(nNodes)
  if (!is.null(packages)) {
    parallel::clusterExport(cl, c("packages"), envir = environment())
    invisible(utils::capture.output(
      parallel::clusterEvalQ(cl, sapply(packages, require, character.only = TRUE))
      )
    )
  }
  if (!is.null(export)) {
    parallel::clusterExport(cl, export, envir = envir)
  }
  return(cl)
}

.nNodes <- function(nNodes) {
  if (is.logical(nNodes)) {
    if (nNodes) {return(parallel::detectCores() - 1)} else {return(1)}
  }
  if (is.numeric(nNodes)) {
    return(min(parallel::detectCores(), round(nNodes), na.rm = TRUE))
  } else {
    stop(sprintf("'nNodes' must be of class numeric or logical, not %s", class(nNodes)))
  }
}

.clpbapply <- function(applyFun, X, FUN, packages, export, MARGIN, nNodes, envir)
{
  start.time <- Sys.time()
  os <- .Platform$OS.type
  pbapply::pboptions(type = "txt", style = 3, char = "=", txt.width = getOption("width") - 10)
  if (nNodes <= 1) {
    if (!is.null(MARGIN)) {
      res <- suppressMessages(applyFun(X = X, FUN = FUN, MARGIN = MARGIN))
    } else {
      res <- suppressMessages(applyFun(X = X, FUN = FUN))
    }
  } else {
    if (os == "unix") {
      message(sprintf("  |ForkCluster (%s): Running %s processes in parallel", os, nNodes))
      cl <- parallel::makeForkCluster(nNodes)
      if (!is.null(MARGIN)) {
        res <- applyFun(cl = cl, X = X, FUN = FUN, MARGIN = MARGIN)
      } else {
        res <- applyFun(cl = cl, X = X, FUN = FUN)
      }
    } else {
      message(sprintf("  |PSOCKcluster (%s): Running %s processes in parallel", os, nNodes))
      cl <- .setUpPSOCKcluster(nNodes = nNodes, packages = packages, export = export, envir = envir)
      if (!is.null(MARGIN)) {
        res <- applyFun(cl = cl, X = X, FUN = FUN, MARGIN = MARGIN)
      } else {
        res <- applyFun(cl = cl, X = X, FUN = FUN)
      }
    }
    parallel::stopCluster(cl)
  }
  message(sprintf("  |Elapsed time: %ss", round(as.numeric(Sys.time()) - as.numeric(start.time), 1)))
  return(res)
}

#' Parallel lapply with progressbar
#'
#' Function detects the operating system and chooses the approximate kind of process for parallelizing the task:
#' Windows: PSOCKCluster, Unix: Forking.
#'
#' @param X a vector (atomic or list) or an expression object. Other objects (including classed objects) will be coerced by base::as.list
#' @param FUN function, the function to be applied to each element of X
#' @param packages character vector, Only relevant for Windows: the packages needed in the function provided, eg. c("MASS", "data.table")
#' @param export character vector, Only relevant for Windows: the varibales needed in the function provided, eg. c("df", "vec")
#' @param envir environment, Only relevant for Windows: Environment from which the variables should be exported from
#' @param nNodes numeric, Number of processes to start (unix: best to fit with the available Cores)
#'
#' @return
#' A list with the results.
#' @export
#'
#' @examples
#' square <- function(x){x*x}
#' l <- parpblapply(X = 1:1000, FUN = square, export = c("square"), nNodes = 2)
parpblapply <- function(X, FUN, packages = NULL, export = NULL, envir = environment(),
                        nNodes = parallel::detectCores() - 1)
{
  .clpbapply(applyFun = pbapply::pblapply,
            X = X, FUN = FUN, packages = packages, export = export,
            MARGIN = NULL, nNodes = nNodes, envir = envir)
}

#' Parallel sapply with progressbar
#'
#' Function detects the operating system and chooses the approximate kind of process for parallelizing the task:
#' Windows: PSOCKCluster, Unix: Forking.
#'
#' @param X a vector (atomic or list) or an expression object. Other objects (including classed objects) will be coerced by base::as.list.
#' @param FUN function, the function to be applied to each element of X
#' @param packages character vector, Only relevant for Windows: the packages needed in the function provided, eg. c("MASS", "data.table")
#' @param export character vector, Only relevant for Windows: the varibales needed in the function provided, eg. c("df", "vec")
#' @param envir environment, Only relevant for Windows: Environment from which the variables should be exported from
#' @param nNodes numeric, Number of processes to start (unix: best to fit with the available Cores)
#'
#' @return
#' A vector with the results.
#' @export
#'
#' @examples
#' square <- function(x){x*x}
#' s <- parpbsapply(X = 1:1000, FUN = square, export = c("square"), nNodes = 2)
parpbsapply <- function(X, FUN, packages = NULL, export = NULL, envir = environment(),
                        nNodes = parallel::detectCores() - 1)
{
  .clpbapply(applyFun = pbapply::pbsapply,
            X = X, FUN = FUN, packages = packages, export = export,
            MARGIN = NULL, nNodes = nNodes, envir = envir)
}

#' Parallel apply with progressbar
#'
#' Function detects the operating system and chooses the approximate kind of process for parallelizing the task:
#' Windows: PSOCKCluster, Unix: Forking.
#'
#' @param X an array, including a matrix.
#' @param FUN function, the function to be applied to each element of X
#' @param MARGIN a vector giving the subscripts which the function will be applied over. E.g., for a matrix 1 indicates rows, 2 indicates columns, c(1, 2) indicates rows and columns. Where X has named dimnames, it can be a character vector selecting dimension names.
#' @param packages character vector, Only relevant for Windows: the packages needed in the function provided, eg. c("MASS", "data.table")
#' @param export character vector, Only relevant for Windows: the varibales needed in the function provided, eg. c("df", "vec")
#' @param envir environment, Only relevant for Windows: Environment from which the variables should be exported from
#' @param nNodes numeric, Number of processes to start (unix: best to fit with the available Cores)
#'
#' @return
#' Returns a vector or array or list of values obtained by applying a function to margins of an array or matrix.
#'
#' @export
#'
#' @examples
#' n <- 1000
#' df <- data.frame(
#' x = seq(1, n, 1),
#' y = -seq(1, n, 1)
#' )
#' a <- parpbapply(X = df, FUN = sum, MARGIN = 1, nNodes = 2)
parpbapply <- function(X, FUN, MARGIN, packages = NULL, export = NULL, envir = environment(),
                       nNodes = parallel::detectCores() - 1)
{
  .clpbapply(applyFun = pbapply::pbapply,
            X = X, FUN = FUN, packages = packages, export = export,
            MARGIN = MARGIN, nNodes = nNodes, envir = envir)
}
