## data model for simulation case 1
## y_i \sim Poisson( \exp(\eta_i) * E_i )
##  \eta_i = \alpha + \beta * X_i + s_i

## data scenario design
##  E_i \sim Gamma(shape = 4, scale = 0.001)
##     E(E_i)   = 4/0.001   = 4000
##     std(E_i) = 4/0.001^2 = 2000
##  x_i \sim U(0, 1)
##     E(x_i) = 0.5
##     V(x_i) = 1/12
##   \alpha = -5
##   \beta = 1.0
##  s \sim N(0, Q(\theta)^{-1})
##     Q(\theta) = \tau(\kappa^2*I - R)
##             | n_i, if i = j
##      R_ij = |  -1, if i ~ j
##             |   0, otherwise
##   \theta = {\kappa, \tau} = {1, 1}

if(!any(ls() == "sample.id")) {

    alpha <- -5
    d.fixed <- 1/3
    Scenarios <- list(id = 1:3)
    Scenarios$label <- LETTERS[Scenarios$id]
    Scenarios$beta <- c(1.0, 0.0, 1)
    Scenarios$include.s <- c(FALSE, TRUE, TRUE)

    mm <- list(
        Ma = . ~ Intercept(1) + x,
        Mb = . ~ Intercept(1) +
            f(i, model = "besagproper", graph = graph,
              hyper = list(theta1 = pprior,
                           theta2 = list(
                               initial = log(d.fixed),
                               fixed = TRUE))))
    mm$Mc <- update(mm$Mb, . ~ . + x)
    mm
    (nM <- length(mm))
    
### getting the map
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
    R <- Diagonal(ns, nns) - graph
    Q1 <- Diagonal(ns, d.fixed) + R
    
    exp(mean(log(diag(chol2inv(chol(Q1))))))
    
    nd <- ns
    dataf <- list(
        i = 1:nd)

## covariates
    xx <- matrix(runif(nSim * nd), nd)
    
    ## spatial 
    xs <- inla.qsample(nSim, Q1)
    
    summary(apply(xs, 2, sd))
    sd(runif(1e5))

} else {
    
### simulate the outcome
    dataf$E <- rgamma(nd, 4, 0.001)
    dataf$x <- xx[, sample.id]
    eta.s <- alpha +
        Scenario$beta * dataf$x +
        Scenario$include.s * xs[, sample.id]
    dataf$y <- rpois(
        n = nd,
        lambda = exp(eta.s) * dataf$E
    )
    summary(dataf$y)

    ## likelihood model

    lhood <- like(
        y ~ .,
        data = dataf,
        family = "poisson",
        E = E
    )

}
