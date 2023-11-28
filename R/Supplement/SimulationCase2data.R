## Data model for simulation case 3
## y_i \sim Poisson( \exp(\eta_i) * E_i )
##  \eta_i = \alpha + v_t + s_i + d_{it}

## data scenario design
##  E_i \sim Gamma(shape = 4, scale = 0.001)
##     E(E_i)   = 4/0.001   = 4000
##     std(E_i) = 4/0.001^2 = 2000
##   \alpha = -5
##  v \sim RW1(\tau_v = 1)
##  s \sim N(0, Q(\theta)^{-1})
##     Q(\theta) = \tau(\kappa^2*I - R)
##             | n_i, if i = j
##      R_ij = |  -1, if i ~ j
##             |   0, otherwise
##   \theta = {\kappa, \tau} = {1, 1}
## d \sim N(0, | Q1 or Q4 )
##     Q1 = \tau_d Kronecker( I_{nt},       I_{ns})
##     Q4 = \tau_d Kronecker( R_{rw1_{nt}},      R)

if(!any(ls() == "sample.id")) {
    
    alpha <- -5
    d.fixed <- 1/3
    rho.fixed <- 0.9
    Scenarios <- list(id = 1:3)
    Scenarios$label <- LETTERS[Scenarios$id]
    Scenarios$stcase <- c(0, 1, 4)
    
    mm <- list(
        Ma = . ~ Intercept(1) +
            time(tt, model = "ar1", 
                 hyper = list(theta1 = pprior,
                              theta2 = list(
                                  initial = log((1 + rho.fixed)/
                                                (1 - rho.fixed)),
                                  fixed = TRUE))) +
            spatial(ii, model = "besagproper", graph = graph,
                    hyper = list(theta1 = pprior,
                                 theta2 = list(
                                     initial = log(d.fixed),
                                     fixed = TRUE))))
    mm$Mb <- update(
        mm$Ma,
        . ~ . +
            spacetime(it, model = "iid",
                      hyper = list(theta = pprior)))
    mm$Mc <- update(
        mm$Ma,
        . ~ . +
            spacetime(it, model = "generic0",
                      Cmatrix = Q4,
                      extraconstr = st4constr,
                      hyper = list(theta = pprior)))
    mm
    (nM <- length(mm))

### Data setup 
    map <- st_read(system.file(
        "shapes/sids.shp", package="spData")[1],
        quiet=TRUE)
    
    ## spatial setup
    nbl <- spdep::poly2nb(map)
    (ns <- length(nbl))
    nns <- spdep::card(nbl)
    graph <- sparseMatrix(
        i = rep(1:ns, nns),
        j = unlist(nbl[nns>0]), x = 1L)
    
    ## Scale the spatial structure matrix
    s1constr <- list(A = matrix(1, 1, ns), e = 0)
    Rs <- inla.scale.model(
        Q = Diagonal(ns, nns) - graph, 
        constr = s1constr
    )
    
    ## time setup
    nt <- 30
    t1constr <- list(A = matrix(1, 1, nt), e = 0)
    Rt <- inla.scale.model(
        Q = INLA:::inla.rw1(nt), 
        constr = t1constr
    )
    Qt <- sparseMatrix(
        i = c(1:nt, 2:nt, 1:(nt-1)),
        j = c(1:nt, 1:(nt-1), 2:nt),
        x = c(rep(1 + rho.fixed^2, nt),
              rep(-rho.fixed, 2 * (nt-1))) / (1 - rho.fixed^2)
    )
    
    ## spacetime
    nd <- nt * ns
    Q1 <- Diagonal(nd)
    Q4 <- kronecker(Rt, Rs)

    ## setup the constraints
    tconstr <- kronecker(diag(nt), matrix(1, 1, ns))
    sconstr <- kronecker(matrix(1,1,nt), diag(ns))
    st4constr <- list(A = rbind(tconstr[-1, ], sconstr),
                      e = rep(0, nt - 1 + ns))
    
### sampling the random effects
    xt <- inla.qsample(nSim, Q=Qt)
    xs <- inla.qsample(nSim, Diagonal(ns, d.fixed) + Rs)
    xst <- inla.qsample(nSim, Q4 + Diagonal(nd, 1e-7), constr = st4constr)
    
    summary(apply(xt, 2, sd))
    summary(apply(xs, 2, sd))
    summary(apply(xst, 2, sd))
    
    dataf <- list(
        ii = rep(1:ns, nt),
        tt = rep(1:nt, each = ns),
        it = 1:nd)
    
} else {
    
### simlate the outcome
    dataf$E <- rgamma(ns, 4, 0.001)[dataf$ii]
    eta.s <- alpha +
        xt[dataf$tt, sample.id] +
        xs[dataf$ii, sample.id]
    if(Scenario$id == 2)
        eta.s <- eta.s + rnorm(nd, 0, 1)
    if(Scenario$id == 3)
        eta.s <- eta.s + xst[, sample.id]
    
    dataf$y <- rpois(
        n = nd,
        lambda = exp(eta.s) * dataf$E
    )
    summary(dataf$y)
    
### likelihood model 
    lhood <- like(
        y ~ .,
        family = "poisson",
        data = dataf,
        E = E)

}
