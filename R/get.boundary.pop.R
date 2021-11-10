#'
#' Generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' Use this function to generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' @usage get.boundary.pop(target, n.cohort, cohortsize)
#'
#' @param target the target DLT rate
#' @param n.cohort the total number of cohorts
#' @param cohortsize the cohort size
#'
#' @details TBD
#'
#' @return  \code{get.boundary.pop()} returns a list object, including the corresponding decision tables
#'          \code{$out.boundary} and \code{$out.full.boundary}.
#'
#'
#' @note TBD
#'
#' @references TBD
#'
#'
#' @examples
#'
#' ## get the dose escalation and deescalation boundaries for BOIN design with
#' ## the target DLT rate of 0.3, maximum sample size of 30, and cohort size of 3
#' bound <- get.boundary.pop(n.cohort = 10, cohortsize = 3, target=0.5)
#' summary(bound) # get the descriptive summary of the boundary
#' plot(bound)    # plot the flowchart of the design with boundaries
#'
#' @import stats
#' @export

get.boundary.pop = function(target,n.cohort,cohortsize){

  if (target < 0.05) {
    stop("the target is too low! ")
  }
  if (target > 0.6) {
    stop("the target is too high!")
  }

  out.boundary <- rbind(cohortsize * (1:n.cohort),
                        bound(target=target,n.cohort=n.cohort,cohortsize = cohortsize))
  total = n.cohort*cohortsize
  out.full.boundary <- rbind(1:total,
                             bound(target=target,n.cohort=total,cohortsize = 1))

  rownames(out.boundary) <- c("Number of patients treated",
                              "Escalation if # of DLT (C1) <=","Deescalation if # of DLT (C2) >=",
                              "Subtherapeutic exclusion if # of DLT (E1) <=",
                              "Overly toxic exclusion if # of DLT (E2) >=")

  rownames(out.full.boundary) <- c("Number of patients treated",
                                   "Escalation if # of DLT (C1) <=","Deescalation if # of DLT (C2) >=",
                              "Subtherapeutic exclusion if # of DLT (E1) <=",
                              "Overly toxic exclusion if # of DLT (E2) >=")

  out = list(out.boundary = out.boundary,
             out.full.boundary = out.full.boundary)
  class(out)<-"pop"
  return(out)
}
#
# bd <-  get.boundary.pop(n.cohort = 10, cohortsize = 3, target=0.5)
# bd



