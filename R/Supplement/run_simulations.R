library(sf)
library(spdep)
library(INLA)
library(inlabru)

hfix1 <- list(
    initial = 0, fixed = TRUE)
pprior <- list(
    prior = "pc.prec",
    param = c(1, 0.05))

ctrc <- list(po = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
ctri <- list(int.strategy = "ccd")
    
fCollectResults <- function(r) {
    nsl <- function(x) {
        x[x<sqrt(.Machine$double.eps)] <- sqrt(.Machine$double.eps)
        return(-sum(log(x)))
    }
    return(c(dic = r$dic$dic,
             waic = r$waic$waic,
             po = -sum(log(r$po$po)),
             cpo = nsl(r$cpo$cpo),
             sapply(r$gcvlist, function(cv)
                 mean(sapply(cv$groups, function(g) length(g$idx)))),
             sapply(r$gcvlist, function(cv) nsl(cv$cv))))
}

SimulationCase <- 1
inla.setOption(num.threads = "2:1")

makeDataFile <- paste0("SimulationCase", SimulationCase, "data.R")

nSim <- 100
system.time(source(makeDataFile))

Outfiles <- 
    paste0("out_files/SimulationCase",
           SimulationCase, Scenarios$label, 
           ".txt")

mlevels1 <- c("m1" = 1, "m5" = 5, "m10" = 10, "m30" = 30)
mlevels <- setdiff(mlevels1, 1)

gcvlnams <- c(
    if(any(mlevels1==1)) "LOGCV1" else NULL, 
    paste0("LOGCV",
           rep(mlevels, each = length(mm) + 1),
           c(names(mm), "u")))
gcvlnams

for(ofl in Outfiles) {
    if (!file.exists(ofl)) {
### Write head of the output file 
        cat("s", "m", "DIC", "WAIC",
            "PO", "LOOCV",
            gsub("LOGCV", "m", gcvlnams),
            gcvlnams, "\n", file = ofl)
    }
}

for(sample.id in 1:nSim) {

    for(sc.id in Scenarios$id) {
        
        Outfile <- Outfiles[sc.id]
        
        Scenario <- lapply(
            Scenarios, function(x) x[sc.id])
    
        source(makeDataFile)
        
### fit models
        fits <- lapply(
            mm, bru, lhood, 
            options=list(
                control.compute = ctrc,
                control.inla = ctri
            )
        )

### self automatic group CV
        for(i in 1:nM) {
            fits[[i]]$gcvlist <- vector("list", length(gcvlnams))
            names(fits[[i]]$gcvlist) <- gcvlnams
            for(m in mlevels1) {
                if(m==1) {
                    fits[[i]]$gcvlist[[1]] <- inla.group.cv(
                        result = fits[[i]],
                        num.level.sets = -1,
                        size.max = 50)
                } else {
                    igr <- grep(m, gcvlnams)[1] + i -1
                    fits[[i]]$gcvlist[[igr]] <- inla.group.cv(
                        result = fits[[i]],
                        num.level.sets = m,
                        size.max = 50)
                }
            }
        }

### use each other model automatic group
        for(m in mlevels) {
            for(i in 1:nM) {
                for(j in 1:nM) {
                    if(!(i==j)) {
                        igr <- grep(m, gcvlnams)[j]
                        fits[[i]]$gcvlist[[igr]] <- inla.group.cv(
                            result = fits[[i]],
                            num.level.sets = m,
                            group.cv = fits[[j]]$gcvlist[[igr]],
                            size.max = 50)
                    }
                }
            }
        }
        
### union of the groups
        uGroups <- lapply(mlevels, function(m) {
            igr <- intersect(grep(m, gcvlnams),
                             sapply(names(mm), grep, gcvlnams))
            lapply(1:nd, function(i) {
                u <- union(fits[[1]]$gcvlist[[igr[1]]]$groups[[i]]$idx,
                           fits[[1]]$gcvlist[[igr[2]]]$groups[[i]]$idx)
                if(length(igr)>2)
                    for(k in 3:length(igr))
                        u <- union(u,
                                   fits[[1]]$gcvlist[[igr[k]]]$groups[[i]]$idx)
                return(u)
            })
        })

### use the union groups    
        for(i in 1:nM) {
            for(g in 1:length(mlevels)) {
                igr <- tail(grep(mlevels[g], gcvlnams), 1)
                fits[[i]]$gcvlist[[igr]] <- inla.group.cv(
                    result = fits[[i]],
                    num.level.sets = mlevels[g],
                    groups = uGroups[[g]])
            }
        }
        
        if(FALSE) {
            
            sapply(fits, function(r) 
                sapply(r$gcvlist, function(cv) {
                    if(is.null(cv$groups)) NA
                    else mean(sapply(cv$groups, function(g) length(g$idx)))
                }))

        }

        r <- sapply(fits, fCollectResults)
        
        for(i in 1:nM) 
            cat(sample.id, i, r[, i], "\n",
                file = Outfile, 
                append = TRUE)
        
    }
    
}
