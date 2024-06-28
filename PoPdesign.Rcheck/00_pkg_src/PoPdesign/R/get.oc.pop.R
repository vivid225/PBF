#' Operating characteristics for single-agent trials
#'
#' Generate the operating characteristics of the PoP design by simulating trials.
#'
#' @usage get.oc.pop(target,n.cohort,cohortsize,titration,skeleton,n.trial,cutoff,cutoff_e,
#'                      risk.cutoff,earlyterm,start,seed)
#'
#' @param target the target DLT rate
#' @param n.cohort the total number of cohorts
#' @param cohortsize the cohort size
#' @param titration default is TRUE. Set \code{titration=TRUE} to perform dose
#'                  escalation with cohort size = 1 to accelerate dose escalation
#'                  at the beginning of the trial.
#' @param skeleton a vector containing the true toxicity probabilities of the
#'                 investigational dose levels.
#' @param n.trial the total number of trials to be simulated
#' @param cutoff the cutoff for the predictive Bayes Factor (PrBF). Users can specify either a value or a function
#'               for cutoff. If PrBF < cutoff, we assign the next cohort of patients to an adjacent dose based on observed DLT.
#'               Otherwise, the evidence is in favor of \eqn{H_{0j}} and we need to retain the current dose.
#' @param cutoff_e the cutoff for the dose exclusion rule. If \eqn{PrBF_{0,1}<E(n_j)}, the evidence is in favor of \eqn{H_{1j}}. If \eqn{\hat{\pi}_j < \phi},
#'                 the current dose is deemed as subtherapeutic and we exclude the current dose and lower doses; If \eqn{\hat{\pi}_j > \phi}, the current dose
#'                 is overly toxic and we exclude the current dose and higher doses.
#' @param risk.cutoff the cutoff to eliminate an over/under toxic dose.
#'                  We recommend the default value of (\code{risk.cutoff=0.8}) for general use.
#' @param earlyterm the early termination parameter.
#' @param start specify the starting dose level. Default value is 1.
#' @param seed the seed for random number generation. Default is 123.
#'
#' @import Iso
#'
#' @details TBD
#'
#' @return \code{get.oc.pop()} returns the operating characteristics of the PoP design as a list,
#'        including:
#'
#'        (1) selection percentage at each dose level (\code{$sel.pct}),
#'
#'        (2) the number of patients treated at each dose level (\code{$num.p}),
#'
#'        (3) the number of toxicities observed at each dose level (\code{$num.tox}),
#'
#'        (4) the average number of toxicities,
#'
#'        (5) the average number of patients,
#'
#'        (6) the percentage of early stopping without selecting the MTD (\code{$early}),
#'
#'        (7) risk of underdosing 80\% or more of patients (\code{$risk.under}),
#'
#'        (8) risk of overdosing 80\% or more of patients (\code{$risk.over})
#'
#' @references Brunk, H., Barlow, R. E., Bartholomew, D. J. & Bremner, J. M (1972, ISBN-13: 978-0471049708).
#'
#' @examples
#'
#' ## get the operating characteristics for single-agent trials
#' oc <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,titration=TRUE,
#'                  cutoff=2.5,cutoff_e=5/24,
#'                  skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
#'                      risk.cutoff=0.8,earlyterm=TRUE,start=1, seed=123)
#'
#' summary(oc) # summarize design operating characteristics
#' plot(oc)
#'
#' @export
#'

get.oc.pop = function(target,n.cohort,cohortsize,titration=TRUE,
                      skeleton,n.trial=1000,
                      cutoff=2.5,cutoff_e=5/24,
                      risk.cutoff=0.8,earlyterm=TRUE,start=1, seed=123){

  set.seed(seed)

  iso.pop <- function(p1,p0){
    l <- which(p0>0)
    p <- (p1[l]+0.05)/(p0[l]+0.1)
    p.var = (p1[l] + 0.05) * (p0[l] - p1[l] + 0.05)/((p0[l] +
                                                        0.1)^2 * (p0[l] + 0.1 + 1))
    p.iso <- pava(p, wt = 1/p.var)
    p.iso = p.iso + (1:length(p.iso)) * 1e-10

    ## eliminate dose based on posterior probability
    elimi = rep(0, K)
    for (i in 1:K) {
      if (1 - pbeta(phi, p1[i] + 1, p0[i] - p1[i] + 1) > 0.95) {
        elimi[i:K] = 1
        break
      }
    }
    m <- which(elimi!=1)

    l <- l[m]
    p.iso <- p.iso[m]
    if(length(l)==0) {return(99)}
    l[sort(abs(p.iso - phi), index.return = T)$ix[1]]
  }

  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }

  cont <- function(x,n.level){
    ret <- rep(0,n.level+1)
    for(i in c(1:n.level)){
      ret[i] <- mean(x==i,na.rm = T)
    }
    ret[n.level+1] <- mean(x==99,na.rm = T)
    ret
  }

  scene <- function(target,K){
    MTD <- sample(K,1)
    M <- rbeta(1,max(K-MTD,0.5),1)
    B <- target+(1-target)*M
    skeleton <- rep(0,K)
    if(MTD == 1){
      skeleton[2:K] <- sort(runif(K-1,target,B))
      skeleton[1] <- runif(1,max(0,2*target-skeleton[2]),skeleton[2])
    }else if(MTD == K){
      skeleton[1:(K-1)] <- sort(runif(K-1,0,target))
      skeleton[K] <- runif(1,skeleton[K-1],min(2*target-skeleton[K-1],B))
    }else{
      skeleton[1:(MTD-1)] <- sort(runif(MTD-1,0,target))
      skeleton[(MTD+1):K] <- sort(runif(K-MTD,target,B))
      d <- min(target-skeleton[MTD-1],skeleton[MTD+1]-target)
      skeleton[MTD] <- runif(1,target-d,target+d)
    }
    round(skeleton,digits = 2)
  }


  ## Titration
  trial_titration <- function(target,lower,upper,elim.lower,elim.upper,skeleton,start=1,
                              earlyterm,n.trial,n.cohort,cohortsize,risk.cutoff=0.8)
  {
    n <- n.cohort*cohortsize
    true.mtd <- which.min(abs(skeleton-target))
    K <- length(skeleton)
    mtd <- rep(NA,n.trial)
    num.p <- num.tox <- matrix(nrow = n.trial,ncol = K)
    risk.over <- rep(NA,n.trial)
    risk.under <- rep(NA,n.trial)
    early <- rep(0,n.trial)

    for(count in 1:n.trial){
      toxic <- matrix(nrow = n,ncol = K)
      for(i in 1:n){
        toxic[i,1] <- rbinom(1,1,prob = skeleton[1])
        for(j in 2:K){
          if(toxic[i,j-1]==1){
            toxic[i,j] <- 1
          }else{
            target1 <- skeleton[j-1]
            target2 <- skeleton[j]
            toxic[i,j] <- rbinom(1,1,prob = (target2-target1)/(1-target1))
          }
        }
      }

      # }
      # Starting dose
      start.dose <- start

      dose.treated <- rep(0,K)
      dose.dlt <- rep(0,K)
      dose.next <- start.dose
      dose.elim <- rep(1,K)

      s <- n

      while(s>0){
        dose.treated[dose.next] <- dose.treated[dose.next]+1
        dlt <- toxic[s,dose.next]
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-1

        if(dlt==1){
          if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){ # if observed dlt <= boundary, escalate
            dose.next <- min(K,dose.next+1)
          }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
            dose.next <- max(dose.next-1,1)
          }
          break
        }else{
          if(dose.next==K){
            break
          }else{
            dose.next <- min(K,dose.next+1)
          }
        }
      }

      while(s>0){
        s.tr <- min(s,cohortsize)
        dose.treated[dose.next] <- dose.treated[dose.next]+s.tr
        dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-s.tr

        ## Exclusion decision
        if (earlyterm){
          if(dose.dlt[dose.next]<=elim.lower[dose.treated[dose.next]]){
            dose.elim[1:dose.next] <- -1 ## Exclude for being subtherapeutic
            if(sum(dose.elim==1)==0){
              early[count] <- 1
              # mtd[count] <- dose.next
              break
            }
          }
          if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
            dose.elim[dose.next:K] <- 0  ## Exclude for being overly toxic
            if(sum(dose.elim==1)==0){
              early[count] <- 1
              # mtd[count] <- dose.next
              break
            }
          }
        }


        ## Transition decision
        if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){
          if(dose.next < K){
            if(dose.elim[dose.next+1]==1){
              dose.next <- dose.next+1
            }
          }
        }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
          if(dose.next > 1){
            if(dose.elim[dose.next-1]==1){
              dose.next <- dose.next-1
            }
          }
        }
      }


      if(is.na(mtd[count])){
        mtd[count] <- iso.pop(dose.dlt,dose.treated)
      }
      if(is.na(mtd[count])){
        next
      }else{
        num.p[count,] <- dose.treated
        num.tox[count,] <- dose.dlt
        risk.over[count] <- 0
        risk.under[count] <- 0
        if(true.mtd==1){
          if(sum(dose.treated[2:K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }else if(true.mtd==K){
          if(sum(dose.treated[1:(K-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }
        }else{
          if(sum(dose.treated[1:(true.mtd-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }else if(sum(dose.treated[(true.mtd+1):K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }
      }
    }
    return(list(num.p=num.p,num.mtd=mtd,early=early,num.tox=num.tox,
                risk.over=risk.over,risk.under=risk.under))
  }

  ## Simulate trials, with early termination rules
  trial_early <- function(target,lower,upper,elim.lower,elim.upper,skeleton,start=1,
                          n.trial,n.cohort,cohortsize,risk.cutoff=0.8)
  {
    n <- n.cohort*cohortsize
    true.mtd <- which.min(abs(skeleton-target))
    K <- length(skeleton)
    mtd <- rep(NA,n.trial)
    num.p <- num.tox <- matrix(nrow = n.trial,ncol = K)
    risk.over <- rep(NA,n.trial)
    risk.under <- rep(NA,n.trial)
    early <- rep(0,n.trial)

    for(count in 1:n.trial){
      toxic <- matrix(nrow = n,ncol = K)
      for(i in 1:n){
        toxic[i,1] <- rbinom(1,1,prob = skeleton[1])
        for(j in 2:K){
          if(toxic[i,j-1]==1){
            toxic[i,j] <- 1
          }else{
            target1 <- skeleton[j-1]
            target2 <- skeleton[j]
            toxic[i,j] <- rbinom(1,1,prob = (target2-target1)/(1-target1))
          }
        }
      }

      # }
      # Starting dose
      start.dose <- start

      dose.treated <- rep(0,K)
      dose.dlt <- rep(0,K)
      dose.next <- start.dose
      dose.elim <- rep(1,K)

      s <- n

      while(s>0){
        dose.treated[dose.next] <- dose.treated[dose.next]+1
        dlt <- toxic[s,dose.next]
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-1

        if(dlt==1){
          if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){ # if observed dlt <= boundary, escalate
            dose.next <- min(K,dose.next+1)
          }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
            dose.next <- max(dose.next-1,1)
          }
          break
        }else{
          if(dose.next==K){
            break
          }else{
            dose.next <- min(K,dose.next+1)
          }
        }
      }

      while(s>0){
        s.tr <- min(s,cohortsize)
        dose.treated[dose.next] <- dose.treated[dose.next]+s.tr
        dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-s.tr

        ## Exclusion decision
        if(dose.dlt[dose.next]<=elim.lower[dose.treated[dose.next]]){
          dose.elim[1:dose.next] <- -1 ## Exclude for being subtherapeutic
          if(sum(dose.elim==1)==0){
            early[count] <- 1
            #mtd[count] <- dose.next
            break
          }
        }
        if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
          dose.elim[dose.next:K] <- 0  ## Exclude for being overly toxic
          if(sum(dose.elim==1)==0){
            early[count] <- 1
            #mtd[count] <- dose.next
            break
          }
        }

        ## Transition decision
        if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){
          if(dose.next < K){
            if(dose.elim[dose.next+1]==1){
              dose.next <- dose.next+1
            }
          }
        }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
          if(dose.next > 1){
            if(dose.elim[dose.next-1]==1){
              dose.next <- dose.next-1
            }
          }
        }
      }


      if(is.na(mtd[count])){
        mtd[count] <- iso.pop(dose.dlt,dose.treated)
      }
      if(is.na(mtd[count])){
        next
      }else{
        num.p[count,] <- dose.treated
        num.tox[count,] <- dose.dlt
        risk.over[count] <- 0
        risk.under[count] <- 0
        if(true.mtd==1){
          if(sum(dose.treated[2:K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }else if(true.mtd==K){
          if(sum(dose.treated[1:(K-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }
        }else{
          if(sum(dose.treated[1:(true.mtd-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }else if(sum(dose.treated[(true.mtd+1):K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }
      }
    }
    return(list(num.p=num.p,num.mtd=mtd,early=early,num.tox=num.tox,
                risk.over=risk.over,risk.under=risk.under))
  }


  ## Simulate trials, without early termination rules
  trial <- function(target,lower,upper,skeleton,risk.cutoff=0.8,
                    n.trial,n.cohort,cohortsize)
  {
    n <- n.cohort*cohortsize
    true.mtd <- which.min(abs(skeleton-target))
    K <- length(skeleton)
    mtd <- rep(NA,n.trial)
    num.p <- num.tox <- matrix(nrow = n.trial,ncol = K)
    risk.over <- rep(NA,n.trial)
    risk.under <- rep(NA,n.trial)
    early <- rep(0,n.trial)

    for(count in 1:n.trial){
      toxic <- matrix(nrow = n,ncol = K)
      for(i in 1:n){
        toxic[i,1] <- rbinom(1,1,prob = skeleton[1])
        for(j in 2:K){
          if(toxic[i,j-1]==1){
            toxic[i,j] <- 1
          }else{
            target1 <- skeleton[j-1]
            target2 <- skeleton[j]
            toxic[i,j] <- rbinom(1,1,prob = (target2-target1)/(1-target1))
          }
        }
      }

      # }
      # Starting dose
      start.dose <- start

      dose.treated <- rep(0,K)
      dose.dlt <- rep(0,K)
      dose.next <- start.dose
      dose.elim <- rep(1,K)

      s <- n

      while(s>0){
        dose.treated[dose.next] <- dose.treated[dose.next]+1
        dlt <- toxic[s,dose.next]
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-1

        if(dlt==1){
          if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){ # if observed dlt <= boundary, escalate
            dose.next <- min(K,dose.next+1)
          }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
            dose.next <- max(dose.next-1,1)
          }
          break
        }else{
          if(dose.next==K){
            break
          }else{
            dose.next <- min(K,dose.next+1)
          }
        }
      }

      while(s>0){
        s.tr <- min(s,cohortsize)
        dose.treated[dose.next] <- dose.treated[dose.next]+s.tr
        dlt <- sum(toxic[s+1-(1:s.tr),dose.next])
        dose.dlt[dose.next] <- dose.dlt[dose.next]+dlt
        s <- s-s.tr

        # ## Exclusion decision
        # if(dose.dlt[dose.next]<=elim.lower[dose.treated[dose.next]]){
        #   dose.elim[1:dose.next] <- -1 ## Exclude for being subtherapeutic
        #   if(sum(dose.elim==1)==0){
        #     early[count] <- 1
        #     #mtd[count] <- dose.next
        #     break
        #   }
        # }
        # if(dose.dlt[dose.next]>=elim.upper[dose.treated[dose.next]]){
        #   dose.elim[dose.next:K] <- 0  ## Exclude for being overly toxic
        #   if(sum(dose.elim==1)==0){
        #     early[count] <- 1
        #     #mtd[count] <- dose.next
        #     break
        #   }
        # }

        ## Transition decision
        if(dose.dlt[dose.next]<=lower[dose.treated[dose.next]]){
          if(dose.next < K){
            if(dose.elim[dose.next+1]==1){
              dose.next <- dose.next+1
            }
          }
        }else if(dose.dlt[dose.next]>=upper[dose.treated[dose.next]]){
          if(dose.next > 1){
            if(dose.elim[dose.next-1]==1){
              dose.next <- dose.next-1
            }
          }
        }
      }


      if(is.na(mtd[count])){
        mtd[count] <- iso.pop(dose.dlt,dose.treated)
      }
      if(is.na(mtd[count])){
        next
      }else{
        num.p[count,] <- dose.treated
        num.tox[count,] <- dose.dlt
        risk.over[count] <- 0
        risk.under[count] <- 0
        if(true.mtd==1){
          if(sum(dose.treated[2:K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }else if(true.mtd==K){
          if(sum(dose.treated[1:(K-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }
        }else{
          if(sum(dose.treated[1:(true.mtd-1)])>risk.cutoff*n){
            risk.under[count] <- 1
          }else if(sum(dose.treated[(true.mtd+1):K])>risk.cutoff*n){
            risk.over[count] <- 1
          }
        }
      }
    }
    return(list(num.p=num.p,num.mtd=mtd,num.tox=num.tox,
                risk.over=risk.over,risk.under=risk.under))

  }



  ## get.oc.pop function starts -----
  phi <- target
  K <- length(skeleton)
  if (titration){
    res = get.boundary.pop(target=target,n.cohort=n.cohort,cohortsize = cohortsize,
                           cutoff=cutoff,K=K,cutoff_e=cutoff_e)$out.full.boundary
  } else {
    if (cohortsize > 1) {
      res = get.boundary.pop(target=target,n.cohort=n.cohort,cohortsize = cohortsize,
                             cutoff=cutoff,K=K,cutoff_e=cutoff_e)$out.boundary
    }
    else {
      res = get.boundary.pop(target=target,n.cohort=n.cohort,cohortsize = cohortsize,
                             cutoff=cutoff,K=K,cutoff_e=cutoff_e)$out.full.boundary
    }
  }

  res[2,which(is.na(res[c(2),]))] <- -Inf
  res[c(4),which(is.na(res[c(4),]))] <- -Inf
  res[3,which(is.na(res[c(3),]))] <- Inf
  res[c(5),which(is.na(res[c(5),]))] <- Inf

  if (earlyterm){
    if (titration){
      out <- trial_titration(target=target,lower=res[2,],upper=res[3,],
                             elim.lower=res[4,],elim.upper=res[5,],
                             skeleton=skeleton,start=start,earlyterm=TRUE,
                             n.trial=n.trial,
                             n.cohort = n.cohort,cohortsize=cohortsize)
    } else {
      out <- trial_early(target=target,lower=res[2,],upper=res[3,],
                         elim.lower=res[4,],elim.upper=res[5,],
                         skeleton=skeleton,start=start,n.trial=n.trial,
                         n.cohort = n.cohort,cohortsize=cohortsize)
    }
    out$sel.pct <- cont(x = out$num.mtd, n.level = length(skeleton))
    out$num.p <- colSums(out$num.p)/n.trial
    out$num.tox <- colSums(out$num.tox)/n.trial
    out$early <- mean(out$early)
    out$risk.over <- mean(out$risk.over)
    out$risk.under <- mean(out$risk.under)

  } else {

    if (titration){
      out <- trial_titration(target=target,lower=res[2,],upper=res[3,],
                             elim.lower=res[4,],elim.upper=res[5,],
                             skeleton=skeleton,start=start,earlyterm=TRUE,
                             n.trial=n.trial,
                             n.cohort = n.cohort,cohortsize=cohortsize)
    } else {
      out <- trial(target=target,lower=res[2,],upper=res[3,],skeleton=skeleton,start=start,
                   n.trial = n.trial,n.cohort = n.cohort,cohortsize=cohortsize)
    }

    out$sel.pct <- cont(x = out$num.mtd, n.level = length(skeleton))
    out$num.p <- out$num.p/n.trial
    out$num.tox <- out$num.tox/n.trial
    out$risk.over <- mean(out$risk.over)
    out$risk.under <- mean(out$risk.under)
  }
  class(out)<-"pop"
  return(out)

}
# octest <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
#                      risk.cutoff=0.8,earlyterm=TRUE,start=1)
# summary(octest)

