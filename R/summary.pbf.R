#'
#' Generate descriptive summary for objects returned by other functions
#'
#' Generate descriptive summary for objects returned by other functions.
#'
#' @param object the object returned by other functions.
#' @param ... ignored arguments
#'
#' @import knitr
#' @details \code{summary()} prints the objects returned by other functions.
#'
#' @return \code{summary()} prints the objects returned by other functions.
#'
#' @author Suyu Liu, Liangcai Zhang and Ying Yuan
#'
#' @examples
#'
#' ###### single-agent trial ######
#'
#' ## summarize the object returned by get.boundary()
#' bound <- get.boundary(target=0.3, ncohort=10, cohortsize=3)
#' summary(bound)
#'
#'
#' ## summarize the object returned by get.oc()
#' oc.single <- get.oc(target=0.3, p.true=c(0.05, 0.15, 0.3, 0.45, 0.6), ncohort=10,
#' cohortsize=3, ntrial=1000)
#' summary(oc.single)
#'
#'
#' ## summarize the object returned by select.mtd()
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
#' summary(selmtd)
#'
#'
#' ###### drug-combination trial######
#'
#' ###### drug-combiation trial to find a single MTD ######
#'
#' ## summarize the object returned by next.comb()
#' n <- matrix(c(3, 0, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol=4, byrow=TRUE)
#' nxt.comb <- next.comb(target=0.25, npts=n, ntox=y, dose.curr=c(1, 1))
#' summary(nxt.comb)
#'
#'
#' ## summarize the object returned by next.comb()
#' n <- matrix(c(3, 3, 0, 0, 0, 0, 0, 0,  0, 0, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), ncol=4, byrow=TRUE)
#' nxt.comb <- next.comb(target=0.25, npts=n, ntox=y, dose.curr=c(1, 2))
#' summary(nxt.comb)
#'
#'
#' ## summarize the object returned by get.oc.comb() when mtd.contour=FALSE
#' p.true <- matrix(c(0.02,0.04,0.08,0.14,
#'                0.08,0.25,0.42,0.48,
#'                0.25,0.45,0.50,0.60), byrow=TRUE, ncol=4)
#'
#' oc.comb <- get.oc.comb(target=0.25, p.true=p.true, ncohort=16, cohortsize=3,
#'            ntrial=100)
#' summary(oc.comb)
#'
#'
#' ## summarize the object returned by select.mtd.comb()
#' n <- matrix(c(6, 3, 0, 0, 6, 24, 9, 0, 0, 0, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 0, 0, 0, 1, 5, 4, 0, 0, 0, 0, 0), ncol=4, byrow=TRUE)
#' sel.comb <- select.mtd.comb(target=0.25, npts=n, ntox=y)
#' summary(sel.comb)
#'
#'
#'
#' ###### drug-combiation trial to find the MTD contour ######
#'
#' ## summarize the object returned by next.subtrial()
#' n <- matrix(c(6, 0,  0, 0,
#'            6, 0, 0, 0,
#'            9, 12, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 0, 0, 0,
#'            1, 0, 0, 0,
#'            2, 3, 0, 0), ncol=4, byrow=TRUE)
#' nxt.trial <- next.subtrial(target=0.3, npts=n, ntox=y)
#' summary(nxt.trial)
#'
#'
#' ## summarize the object returned by get.oc.comb() when mtd.contour=TRUE.
#' p.true <- matrix(c(0.01,0.03,0.10,0.20,0.30,
#'                0.03,0.05,0.15,0.30,0.60,
#'                0.08,0.10,0.30,0.60,0.75), byrow=TRUE, ncol=5)
#'
#' oc.comb <- get.oc.comb(target=0.3, p.true, ncohort=c(10,5,5), cohortsize=3,
#'    n.earlystop=12, startdose=c(1,1),ntrial=100, mtd.contour=TRUE)
#' summary(oc.comb)
#'
#'
#' ## summarize the object returned by select.mtd.comb()
#' n <- matrix(c(6, 9, 24, 0,
#'            6,  24, 9, 0,
#'            12, 18, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 1, 5, 0,
#'            1, 5, 4, 0,
#'            1, 5, 0, 0), ncol=4, byrow=TRUE)
#' sel.comb2 <- select.mtd.comb(target=0.3, npts=n, ntox=y, mtd.contour=TRUE)
#' summary(sel.comb2)
#'
#'
#'
#' @export

summary.pbf<- function (object, ...) {

  ## get.boundary -----
  if (!is.null(object$out.boundary)) {
    cat("The decision boundaries are\n")
    print(object$out.boundary)
    cat("\n")
    cat("A more completed version of the decision boundaries is given by\n")
    print(object$out.full.boundary)
  }

  ## get.oc -----
  if (!is.null(object$sel.pct)) {
    cat("selection percentage at each dose level (%):\n")
    cat(formatC(object$sel.pct*100, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of patients treated at each dose level:\n")
    cat(formatC(object$num.p, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of toxicity observed at each dose level: ")
    cat(formatC(object$num.tox, digits = 1, format = "f"),
        sep = "  ", "\n")
    cat("average number of toxicities:", formatC(sum(object$num.tox),
                                                 digits = 1, format = "f"), "\n")
    cat("average number of patients:", formatC(sum(object$num.p),
                                               digits = 1, format = "f"), "\n")
    cat("percentage of early stopping due to toxicity:",
        formatC(object$early*100, digits = 1, format = "f"),
        "% \n")
    cat("risk of underdosing (>80% of patients treated below the MTD):",
        formatC(object$risk.under*100, digits = 1, format = "f"),
        "% \n")
    cat("risk of overdosing (>80% of patients treated above the MTD):",
        formatC(object$risk.over*100, digits = 1, format = "f"),
        "% \n")
  }


  ## select.mtd -----
  if (!is.null(object$MTD)) {
    tbl <- data.frame(dose = 1:length(object$p_est),
                      dlt = round(object$p_est,2))
    colnames(tbl) <- c("Dose level", "Posterior DLT")
    cat("The MTD is dose level ", object$MTD, "\n")
    kable(tbl, "simple")
  }



}





