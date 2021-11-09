#'
#' Generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' Use this function to generate the optimal dose escalation and deescalation boundaries for conducting the trial.
#'
#' @usage get.boundary.pbf(target, n.cohort, cohortsize)
#'
#' @param target the target DLT rate
#' @param n.cohort the total number of cohorts
#' @param cohortsize the cohort size
#'
#' @details TBD
#'
#' @return  \code{get.boundary()} returns a list object, including the corresponding decision tables
#'          \code{$out.boundary} and \code{$out.full.boundary}.
#'
#'
#' @note TBD
#'
#' @references TBD
#'
#' @seealso
#'
#' @author
#'
#' @examples
#'
#' ## get the dose escalation and deescalation boundaries for BOIN design with
#' ## the target DLT rate of 0.3, maximum sample size of 30, and cohort size of 3
#' bound <- get.boundary.pbf(n.cohort = 10, cohortsize = 3, target=0.5)
#' summary(bound) # get the descriptive summary of the boundary
#' plot(bound)    # plot the flowchart of the design with boundaries
#'
#' @import stats
#' @export

get.boundary.pbf = function(target,n.cohort,cohortsize){

  prbf01 <- function(n,y,target){
    p <- beta(y+2,n-y+1)/beta(1+y,n-y+1)
    dbinom(y,n,prob=target)*exp(1)/dbinom(y,n,prob=p)
  }

  bound <- function(target,n.cohort,cohortsize){
    lower <- upper <- 0
    lower.ex <- upper.ex <- 0
    for(i in 1:n.cohort){
      t <- i*cohortsize
      x <- 0:(i*cohortsize)
      y <- lapply(x,prbf01,n=(i*cohortsize),target=target)
      a <- 0
      for(j in 1:length(x)){
        if(y[[j]]<exp(1)){ #e,exp(1)*(1/t*log(1+t))^(1/(2*t)),exp(1)*(1/t*log(1+t))^(1/(t))
          if(x[j]/(i*cohortsize)<target){
            a[j] <- 1
          }else{
            a[j] <- -1
          }
        }else{
          a[j] <- 0
        }
      }

      if (any(a==1)){
        lower[i] <- max(which(a==1))-1
      } else {
        lower[i] <- NA
      }
      if (any(a==-1)){
        upper[i] <- min(which(a==-1))-1
      } else {
        upper[i] <- NA
      }

      ex <- 0
      for(j in 1:length(x)){
        if(y[[j]]<0.3){
          if(x[j]/(i*cohortsize)<target){
            ex[j] <- 1 # exclude for being subtherapeutic
          }else{
            ex[j] <- -1 # exclude for being too toxic
          }
        }else{
          ex[j] <- 0
        }
      }
      if (any(ex==1)){
        lower.ex[i] <- max(which(ex==1))-1
      } else {
        lower.ex[i] <- NA
      }
      if (any(ex==-1)){
        upper.ex[i] <- min(which(ex==-1))-1
      } else {
        upper.ex[i] <- NA
      }
    }
    b <- rbind(lower,upper,lower.ex,upper.ex)
    return(b)
  }

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
                              "Escalation if # of DLT <=","Deescalation if # of DLT >=",
                              "Subtherapeutic exclusion if # of DLT <=",
                              "Overly toxic exclusion if # of DLT >=")

  rownames(out.full.boundary) <- c("Number of patients treated",
                                   "Escalation if # of DLT <=","Deescalation if # of DLT >=",
                              "Subtherapeutic exclusion if # of DLT <=",
                              "Overly toxic exclusion if # of DLT >=")

  out = list(out.boundary = out.boundary,
             out.full.boundary = out.full.boundary)
  class(out)<-"pbf"
  return(out)
}
#
# bd <-  get.boundary.pbf(n.cohort = 10, cohortsize = 3, target=0.5)
# bd



