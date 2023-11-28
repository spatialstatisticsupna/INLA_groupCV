level_sets <- c(m1 = 1, m3 = 3, m5 = 5, m10 = 10)
tresol <- 1
sresol <- 30
source("prepare_fit.R")

th.ini <- c(4, 4, 10, 10, -2)
print(th.ini)


nm <- length(models)
results <- vector("list", nm)
names(results) <- paste0("u", models)
names(results)

rfls <- paste0(
    rdir,
    "fg_u", models, 
    "_n", ndata,
    "_ns", sspde$n.spde,
    "_nt", tmesh$n,
    ".rds")
rfls

for(i.u in 1:4) {

    cat("Results file:", rfls[i.u], "\n")

### file to save the result
    if(file.exists(rfls[i.u])) {
        
        t0 <- Sys.time()
        cat("Loading the previously saved result ... ")
        rtmp <- readRDS(rfls[i.u])
        print(Sys.time() - t0)
        
    } else {

        t0 <- Sys.time()
        cat("Defining model",
            models[i.u], "... ")
        stmodel <- stModel.define(
            smesh, tmesh, models[i.u],
            control.priors = list(
                prs = c(S0, 0.05),
                prt = c(ifelse(models[i.u]=='121', R0t*2, R0t), 0.05),
                psigma = c(S0, 0.05)),
            constr = TRUE)
        cat(stmodel$f$cgeneric$data$matrices$xx[1:2], " ")
        print(Sys.time() - t0)
                
        cat("Fitting model", models[i.u], "\n")

        rtmp <- bru(
            comps,
            lhood,
            options = list(
                control.mode = list(
                    theta = th.ini, restart = TRUE),
                verbose = TRUE,
                inla.call = "remote",
                control.compute = ctrc)
        )
        
        cat(i.u, "")
        cat(c(th = unname(rtmp$mode$theta)), "\n")
        print(Sys.time() - t0)
        rtmp$cpu <- c(rtmp$cpu, nfn = max(rtmp$misc$nfunc))
        cat(i.u, "cpu, n:", rtmp$cpu, "\n")
        
        t0 <- Sys.time()
        cat(i.u, " computing GCPO ... ")        
        rtmp$.args$verbose <- FALSE
        if(i.u==1) {
            rtmp$gcvlist <- lapply(level_sets, function(m) 
                inla.group.cv(
                    result = rtmp, 
                    num.level.sets = m, 
                    size.max = 50,
                    strategy = "posterior")
                )
            cat("Define the reference automatic groups for \n",
                i.u, "", "model", models[i.u], "\n")            
            reference.groups <- lapply(rtmp$gcvlist, function(rgcv) {
                return(lapply(rgcv$groups, function(x) x$idx))
            })
        } else {
            if(!any(ls() == "reference.groups")) {
                cat("Collect the reference automatic groups from \n",
                    "model", models[1], "\n")
                refgcv <- readRDS(rfls[1])$gcvlist
                reference.groups <- lapply(refgcv, function(gcv) {
                    return(lapply(gcv$groups, function(x) x$idx))
                })
                rm(refgcv)
            } 
            rtmp$gcvlist <- lapply(1:length(level_sets), function(im) 
                inla.group.cv(
                    result = rtmp, 
                    num.level.sets = level_sets[im],
                    groups = reference.groups[[im]],
                    strategy = "posterior")
                )
            names(rtmp$gcvlist) <- names(level_sets)
        }
        print(Sys.time() - t0)
        
        cat("Saving results ... ")
        t0 <- Sys.time()
        saveRDS(rtmp[irsel], rfls[i.u])
        print(Sys.time() - t0)

        rm(stmodel)

    } ## end else fit
    
    results[[i.u]] <- rtmp
    cat("cpu.nfn:", rtmp$cpu, "\n")
    rm(rtmp)
    gc(reset = TRUE)
} ## end for stmodels

