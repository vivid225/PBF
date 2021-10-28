#'
#' Plot the flowchart and simulation results for BOIN designs
#'
#' Plot the objects returned by other functions, including (1) flowchart of BOIN design;
#' (2) operating characteristics of the design, including selesction percentage and the
#' number of patients treated at each dose;
#' (3) the estimate of toxicity probability for each dose and corresponding 95\% credible interval
#'
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User doesn't need to input this parameter.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#'
#' @author
#'
#' @examples
#'
#' ###### single-agent trial ######
#'
#' ## get dose escalation and deescalation boundaries for conducting the trial
#' bound <- get.boundary(target=0.3, ncohort=10, cohortsize=3)
#' plot(bound)
#'
#'
#' ## get the operating characteristics for BOIN single agent trial
#' oc <- get.oc(target=0.3, p.true=c(0.05,0.15,0.3,0.45,0.6),
#'    ncohort=10, cohortsize=3, ntrial=1000)
#' summary(oc)
#' plot(oc)
#'
#'
#' ## select the MTD based on the trial data
#' n <- c(3, 3, 15, 9, 0)
#' y <- c(0, 0, 4, 4, 0)
#' selmtd <- select.mtd(target=0.3, npts=n, ntox=y)
#' summary(selmtd)
#' plot(selmtd)
#'
#'
#' ###### drug-combination trial ######
#'
#' ##### combination trial to find a single MTD  ######
#'
#' ## get the operating characteristics for BOIN combination trial
#' p.true <- matrix(c(0.01,0.03,0.10,0.20,0.30,
#'                 0.03,0.05,0.15,0.30,0.60,
#'                 0.08,0.10,0.30,0.60,0.75), byrow=TRUE, ncol=5)
#'
#' oc.comb <- get.oc.comb(target=0.3, p.true, ncohort=20, cohortsize=3, n.earlystop=12,
#'      startdose=c(1,1),ntrial=100)
#' summary(oc.comb)
#' plot(oc.comb)
#'
#'
#' ## select a MTD based on the trial data
#' n <- matrix(c(3, 5, 0, 0, 0, 7, 6, 15, 0, 0, 0, 0, 4, 0, 0), ncol=5, byrow=TRUE)
#' y <- matrix(c(0, 1, 0, 0, 0, 1, 1, 4, 0, 0, 0, 0, 2, 0, 0), ncol=5, byrow=TRUE)
#' sel.comb <- select.mtd.comb(target=0.3, npts=n, ntox=y)
#' summary(sel.comb)
#' plot(sel.comb)
#'
#'
#' ##### combination trial to find a MTD contour (e.g., multiple MTDs)  #####
#'
#' ## get the operating characteristics for BOIN waterfall design
#' p.true <- matrix(c(0.01, 0.10, 0.20, 0.30,
#'                 0.03, 0.15, 0.30, 0.60,
#'                 0.08, 0.30, 0.60, 0.75), byrow=TRUE, ncol=4)
#'
#' oc.comb2 <- get.oc.comb(target=0.3, p.true, ncohort=c(8,6,6), cohortsize=3, n.earlystop=12,
#'        startdose=c(1,1), ntrial=100, mtd.contour=TRUE)
#' summary(oc.comb2)
#' plot(oc.comb2)
#'
#'
#' ## select the MTD contour based on the trial data
#' n <- matrix(c(6, 9, 24, 0,  6, 24, 9, 0,  12, 18, 0, 0), ncol=4, byrow=TRUE)
#' y <- matrix(c(0, 1,  5, 0,  1,  5, 4, 0,  1, 5, 0, 0), ncol=4, byrow=TRUE)
#' sel.comb2 <- select.mtd.comb(target=0.3, npts=n, ntox=y, mtd.contour=TRUE)
#' summary(sel.comb2)
#' plot(sel.comb2)
#'
#' @import knitr
#' @export

plot.pbf<- function (object){
  ## get.boundary -----
  if (!is.null(object$out.boundary)) {
    # pp <- readPNG("/Users/rrrrrita/Documents/GitHub/PBF/inst/Flowchart/PBF_Flowchart.png")
    # plot.new()
    # pp
    # rasterImage(pp,0,-0.5,1,0.5)
    # include_graphics("/Users/rrrrrita/Documents/GitHub/PBF/inst/Flowchart/PBF_Flowchart.png")
    # plot(pp)
    # print()
    kable(object$out.boundary,"simple")

    }

  ## get.oc -----
  if (!is.null(object$sel.pct)) {
    dt.plot <- data.frame(dose = 1:length(object$num.p),
                          pct = object$sel.pct)
    p <- barplot(dt.plot$pct,names.arg=c("1", "2", "3","4"),
                 xlab="Dose level", ylab="Selection percentage (%)")
    p
  }


  ## select.mtd -----
  if (!is.null(object$MTD)) {
    plot(object$p_est, ylab="DLT rate", xlab="Dose level",pch=20)
    abline(h=object$target, col="Red", lty=3)
  }
}


