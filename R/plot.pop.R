#'
#' Plot the flowchart and simulation results for BOIN designs
#'
#' Plot the objects returned by other functions, including (1) flowchart of BOIN design;
#' (2) operating characteristics of the design, including selection percentage and the
#' number of patients treated at each dose;
#' (3) the estimate of toxicity probability for each dose and corresponding 95\% credible interval
#'
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User doesn't need to input this parameter.
#'
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered
#'
#'
#' @importFrom graphics abline barplot axis
#' @export

plot.pop<- function (x, ...){
  ## get.boundary -----
  if (!is.null(x$out.boundary)) {
    library(magick)

    pp <- readPNG("/Users/rrrrrita/Documents/GitHub/PoP/inst/Flowchart/PoP_Flowchart.png")
    # pp <- readPNG(img)
    # plot(pp)
    # plot.new()
    # pp
    # rasterImage(pp,0,-0.5,1,0.5)
    # include_graphics("/Users/rrrrrita/Documents/GitHub/pop/inst/Flowchart/pop_Flowchart.png")
    # plot(pp)
    # print()
    # kable(object$out.boundary,"simple")
    # library("imager")
    # original = imager::load.image("/Users/rrrrrita/Documents/GitHub/PoP/inst/Flowchart/PoP_Flowchart_scaled.png")
    # d = dim(original)[1:2]
    # d
    #
    # scale = max(d / 400)
    # img = imager::resize(original,
    #                      imager::width(original) / scale,
    #                      imager::height(original) / scale,
    #                      interpolation_type = 6)
    # img
    # plot(original)

    # square_img = imager::pad(img,
    #                          nPix = 400 - height(img),
    #                          axes = "y",
    #                          val = "white")
    # imager::save.image(img, file = "/Users/rrrrrita/Documents/GitHub/PoP/inst/Flowchart/PoP_Flowchart_scaled.png")
    # square_img
    # out.width="200px", fig.retina = 1
    # knitr::include_graphics("/Users/rrrrrita/Documents/GitHub/PoP/inst/Flowchart/PoP_Flowchart.png")

    }

  ## get.oc.pdf -----
  if (!is.null(x$sel.pct)) {
    dt.plot <- data.frame(dose = 1:length(x$num.p),
                          pct = x$sel.pct)
    p <- barplot(dt.plot$pct,names.arg=c("1", "2", "3","4"),
                 xlab="Dose level", ylab="Selection percentage (%)")
    p
  }


  ## select.mtd.pdf -----
  if (!is.null(x$MTD)) {
    plot(x$p_est, ylab="DLT rate", xlab="Dose level",pch=20,xaxt = "n")
    axis(1, at=1:length(x$p_est))
    abline(h=x$target, col="Red", lty=3)
  }
}


