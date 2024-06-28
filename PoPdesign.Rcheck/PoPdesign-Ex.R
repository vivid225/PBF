pkgname <- "PoPdesign"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('PoPdesign')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("get.boundary.pop")
### * get.boundary.pop

flush(stderr()); flush(stdout())

### Name: get.boundary.pop
### Title: Generate the dose escalation and de-escalation boundaries for
###   single-agent trials.
### Aliases: get.boundary.pop

### ** Examples


## get the dose escalation and deescalation boundaries for PoP design with
## the target DLT rate of 0.3, maximum sample size of 30, and cohort size of 3
bound <- get.boundary.pop(target=0.5, n.cohort = 10, cohortsize = 3,
                          cutoff=2.5,K=4,cutoff_e=5/24)
summary(bound) # get the descriptive summary of the boundary
plot(bound)    # plot the flowchart of the design along with decision boundaries




cleanEx()
nameEx("get.oc.pop")
### * get.oc.pop

flush(stderr()); flush(stdout())

### Name: get.oc.pop
### Title: Operating characteristics for single-agent trials
### Aliases: get.oc.pop

### ** Examples


## get the operating characteristics for single-agent trials
oc <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,titration=TRUE,
                 cutoff=2.5,cutoff_e=5/24,
                 skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
                     risk.cutoff=0.8,earlyterm=TRUE,start=1, seed=123)

summary(oc) # summarize design operating characteristics
plot(oc)




cleanEx()
nameEx("select.mtd.pop")
### * select.mtd.pop

flush(stderr()); flush(stdout())

### Name: select.mtd.pop
### Title: Maximum tolerated dose (MTD) selection for single-agent trials
### Aliases: select.mtd.pop

### ** Examples


### select the MTD for PoP trial
n <- c(4, 4, 16, 8, 0)
y <- c(0, 0, 5, 5, 0)
selmtd <- select.mtd.pop(target=0.3,n.pts=n, n.tox=y)
summary(selmtd)
plot(selmtd)




cleanEx()
nameEx("summary.pop")
### * summary.pop

flush(stderr()); flush(stdout())

### Name: summary.pop
### Title: Generate descriptive summary for objects returned by other
###   functions in PoPdesign
### Aliases: summary.pop

### ** Examples

## summarize the results returned by get.boundary.pop()
bound <- get.boundary.pop(n.cohort = 10, cohortsize = 3, target=0.3,
                          cutoff=exp(1), K=3,cutoff_e=exp(-1))
summary(bound)

## summarize the results returned by get.oc.pop()
oc <- get.oc.pop(target=0.3,n.cohort=10,cohortsize=3,titration=TRUE,
                 cutoff=TRUE,cutoff_e=exp(-1),skeleton=c(0.3,0.4,0.5,0.6),n.trial=1000,
                 risk.cutoff=0.8,earlyterm=TRUE,start=1)
summary(oc)

### summarize the results returned by select.mtd.pop()
n <- c(3, 3, 15, 9, 0)
y <- c(0, 0, 4, 4, 0)
selmtd <- select.mtd.pop(target=0.3,n.pts=n, n.tox=y)
summary(selmtd)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
