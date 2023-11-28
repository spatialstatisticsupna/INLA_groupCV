## Data model for simulation case 3
##   y_i = a + v(t_i) + s(s_i) + u(s_i, t_i) + e_i

## data scenario design
##  v \sim RW1(\tau_v = 1)
##  s \sim N(0, Q(\theta)^{-1})
##   \theta = {range, \sigma}
##  u \sim N(0, | Q1 or Q2 )
##   Q1 is from model 102
##   Q2 is from model 121

if(!any(ls() == "sample.id")) {

    nd <- 1000 ## observations
    nt <- 20   ## time points

    Scenarios <- list(id = 1:3)
    Scenarios$label <- LETTERS[Scenarios$id]
    Scenarios$stcase <- c(0, 1, 2)
    
    rho.fixed <- 0.9

    ## models to fit
    ## Ma: a + v(t_i) + s(s_i) + e_i
    mm <- list(
        Ma = . ~ Intercept(1) +
            temporal(time, model = "ar1", 
                     hyper = list(theta1 = pprior,
                                  theta2 = list(
                                      initial = log((1 + rho.fixed)/
                                                    (1 - rho.fixed)),
                                  fixed = TRUE))) +
            spatial(cbind(xloc, yloc),
                    model = sspde, mapper = smapper)
    )
    ## Mb: a + v(t_i) + s(s_i) + u1(s_i, t_i) + e_i
    mm$Mb <- update(
        mm$Ma,
        . ~ . +
            spacetime(list(space = cbind(xloc, yloc),
                           time = time),
                      model = stmodel1))
    ## Mc: a + v(t_i) + s(s_i) + u2(s_i, t_i) + e_i
    mm$Mc <- update(
        mm$Ma,
        . ~ . +
            spacetime(list(space = cbind(xloc, yloc),
                           time = time),
                  model = stmodel2))
    mm
    (nM <- length(mm))
    
    ## temporal setup
    tmesh <- fm_mesh_1d(loc = 1:nt)
    
    ## scaling the temporal structure matrix for v
    Qt <- sparseMatrix(
        i = c(1:nt, 2:nt, 1:(nt-1)),
        j = c(1:nt, 1:(nt-1), 2:nt),
        x = c(rep(1 + rho.fixed^2, nt),
              rep(-rho.fixed, 2 * (nt-1))) / (1 - rho.fixed^2)
    )
    
### spatial setup
    bb <- c(0, 0, 2, 2)
    domain <- cbind(
        c(0, 1, 1, 0, 0) * bb[3],
        c(0, 0, 1, 1, 0) * bb[4]
    )
    
    smesh.sim <- fm_mesh_2d(
        loc.domain = domain,
        max.edge = c(0.1, .5),
        offset = c(0.5, 1),
        cutoff = 0.05
    )

    if(FALSE)
        plot(smesh)
    
    sspde <- inla.spde2.pcmatern(
        mesh = smesh.sim,
        alpha = 2,
        prior.range = c(2.0, NA),
        prior.sigma = c(1.0, NA),
        constr = TRUE
    ); Qs <- inla.spde2.precision(
           spde = sspde,
           theta = log(c(1, 1))
       ); rm(sspde)
    
    ## spacetime setup
    pars2 <- list(m1 = c(1, 10, 1),
                  m2 = c(1, 20, 1))
    Qst1 <- stModel.precision(smesh.sim, tmesh, "102", log(pars2[[1]]))
    Qst2 <- stModel.precision(smesh.sim, tmesh, "121", log(pars2[[2]]))
    c(length(Qst1@x), length(Qst2@x))
    
### sampling the random effects
    xt <- inla.qsample(
        n = nSim,
        Q = Qt)
    
    xs <- inla.qsample(nSim, Qs)
    
    xst1 <- inla.qsample(
        n = nSim, Q = Qst1)
    
    xst2 <- inla.qsample(
        n = nSim, Q = Qst2)
    
    summary(apply(xt, 2, sd))
    summary(apply(xs, 2, sd))
    
    summary(apply(xst1, 2, mean))
    summary(apply(xst1, 2, sd))

    summary(apply(xst2, 2, mean))
    summary(apply(xst2, 2, sd))
    
    rm(Qst1, Qst2); gc(reset = TRUE)
    
    smesh <- fm_mesh_2d(
        loc.domain = domain,
        max.edge = c(0.2, .7),
        offset = c(0.5, 1),
        cutoff = 0.1
    )
    smesh$n
    
    if(FALSE)
        plot(smesh)
    
    sspde <- inla.spde2.pcmatern(
        mesh = smesh,
        alpha = 2,
        prior.range = c(2.0, NA),
        prior.sigma = c(1.0, 0.05),
        constr = TRUE
    )
    smapper <- bru_mapper(smesh)
    
    cat("nt = ", nt, "smesh$n =", smesh.sim$n, 
        "nst =", nt*smesh$n, "and nd =", nd, "\n")
    
    stmodel1 <- stModel.define(
        smesh = smesh,
    tmesh = tmesh,
    model = "102",
    control.priors = list(
        prs = c(1, 0.00),
        prt = c(2, 0.05),
        psigma = c(1, 0.05)
    ),
    constr = TRUE
    )
    
    stmodel2 <- stModel.define(
        smesh = smesh,
        tmesh = tmesh,
        model = "121",
        control.priors = list(
            prs = c(1, 0.00),
            prt = c(4, 0.05),
            psigma = c(1, 0.05)
        ),
        constr = TRUE
    )
    
    rbind(stmodel1$f$cgeneric$data$matrices$x[1:2],
          stmodel2$f$cgeneric$data$matrices$x[1:2])

} else {
    
### simlate the outcome
    dataf <- data.frame(
        xloc = runif(nd, bb[1], bb[3]),
        yloc = runif(nd, bb[2], bb[4]),
        time = sample(1:nt, nd, replace = TRUE)
    )
    summary(dataf)
    
    As <- inla.spde.make.A(
        mesh = smesh.sim,
        loc = cbind(dataf$xloc, dataf$yloc)
    )
    Ast <- inla.spde.make.A(
        mesh = smesh.sim,
        loc = cbind(dataf$xloc, dataf$yloc),
        group = dataf$time,
        group.mesh = tmesh
    )
    
    eta.s <- xt[dataf$time, sample.id] +
        drop(As %*% xs[, sample.id])
    if(sc.id == 2)
        eta.s <- eta.s +
            drop(Ast %*% xst1[, sample.id])
    if(sc.id == 3)
        eta.s <- eta.s +
        drop(Ast %*% xst2[, sample.id])
    
    dataf$y <- eta.s + rnorm(nd, 0, 0.3)
    summary(dataf$y)
    
    lhood <- like(
        y ~ .,
        data = dataf,
        control.family = list(
            hyper = list(theta = pprior))
    )

}
