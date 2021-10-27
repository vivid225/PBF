#'
#' Generate descriptive summary for objects returned by other functions
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#'
#'
#' @details \code{print()} prints the objects returned by other functions.
#'
#' @return \code{print()} prints the objects returned by other functions.
#'
#' @author
#'
#' @examples
#'
#' ###### single-agent trial ######
#'
#' ## sprint the object returned by get.boundary()
#' bound <- get.boundary(target=0.3, ncohort=10, cohortsize=3)
#' print(bound)
#'
#'
#' ## print the object returned by get.oc()
#' oc.single <- get.oc(target=0.3, p.true=c(0.05, 0.15, 0.3, 0.45, 0.6), ncohort=10,
#' cohortsize=3, ntrial=1000)
#' print(oc.single)
#'
#'
#' ## print the object returned by select.mtd()
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
#' print(selmtd)
#'
#'
#' @export

print.pbf<-function(x,...){
  print.default(x)
}
