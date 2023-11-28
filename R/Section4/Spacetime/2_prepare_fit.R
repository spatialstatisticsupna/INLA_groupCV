library(parallel)
library(ggplot2)
library(patchwork)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(fmesher)
library(INLA)
library(INLAspacetime)
library(inlabru)

inla.setOption(
    smtp = "pardiso",
    pardiso.lic = "~/.pardiso.lic",
    num.threads = "12:1", ## (1 + 1 + 3) + 1
    safe = FALSE
)

## local directory
wdir <- getwd()
wdir

## data dir
ddir <- file.path(wdir, "data/")
ddir

## directory for the output files
rdir <- file.path(wdir, "rdsfiles/")
rdir

if(!dir.exists(rdir))
    dir.create(rdir)

## get the UK map
map_ll <- ne_countries(
    scale = 50,
    country = c("ireland", "united kingdom"),
    returnclass = "sf")

## UTM
utm_crs <- "+proj=utm +zone=30 +datum=WGS84 +units=km"
map_utm <- st_transform(map_ll, utm_crs)

### load the data
dataf <- readRDS("data/fgselected.rds")
stations <- readRDS("data/fgstations.rds")

datalocs <- unique(cbind(dataf$xloc, dataf$yloc))
(nlocs <- nrow(stations))
stopifnot(nlocs == nrow(datalocs))
(ndata <- nrow(dataf))
(nt <- max(dataf$time))

## basic model components: \beta_0 + \beta_1 E(s) + r(s) + u(s, t)
comps <- ~ -1 + Intercept(1) + I(elevation/1000) +
    temporal(time, model = "rw1", scale.model = TRUE, hyper = list(theta = pcprec)) +
##    temporal(time, model = "ar1",
  ##           hyper = list(
    ##             theta1 = pcprec,
      ##           theta2 = list(
        ##             prior = "pccor1",
          ##           param = c(0, 0.9))),
            ## constr = TRUE) +
##    spatial(cbind(xloc, yloc),
  ##          model = sspde) + 
    field(list(space = cbind(xloc, yloc), 
               time = time),
          model = stmodel) 


## Define a temporal mesh, with each knot spaced by h,
## where tresol = 1 means one per day.
if(!any(ls() == "tresol"))
	tresol <- 1 # temporal mesh resolution in days
tmesh <- fm_mesh_1d(
  loc = seq(1, nt + tresol/2, tresol),
  degree = 1)

## Building mesh can play around with more values
if(!any(ls() == "sresol"))
    sresol <- 30 ## spatial mesh resolution in kilometers

smesh.fl <- file.path(
    rdir, paste0("smesh_sresol", sresol, ".rds"))
smesh.fl

if(all(file.exists(smesh.fl))) {

    smesh <- readRDS(smesh.fl)

} else {
    
    boundseg <- fm_sp2segment(map_utm)
    smesh <-  fm_mesh_2d(
        loc = datalocs, 
        boundary = boundseg, 
        offset = c(1, 3) * sresol,
        max.edge = c(1, 3) * sresol,
        cutoff = 0.5 * sresol)

    saveRDS(smesh, smesh.fl)
    
}

cat("nt =", tmesh$n, "and ns =", smesh$n, "\n")

## visualize the mesh
if(FALSE)
    ggplot() + theme_minimal() +
        geom_sf(data = map_utm) + 
        gg(smesh) +
        geom_sf(data = stations) +
        xlab("") + ylab("")

(sd(dataf$fg, na.rm = TRUE))
(mean(c(diff(range(dataf$xloc)),
        diff(range(dataf$yloc)))))

S0 <- 2    ## to use in P(sigma > S0) = 0.05
R0s <- 50  ## to use in P(spatial range < R0s) = 0.05
R0t <- 3   ## to use in P(temporal range < R0t) = 0.05
### (in the model 121 for u use 2*R0t)

sspde <- inla.spde2.pcmatern(
    mesh = smesh,
    alpha = 2,
    prior.range = c(R0s, 0.05),
    prior.sigma = c(S0, 0.05),
    constr = TRUE
)

## likelihood precision prior
pcprec <- list(
    prec = list(prior = "pcprec",
                param = c(S0, 0.05)))

## likelihood setup
lhood <- like(
  formula = fg ~ ., 
  family = "gaussian",
  control.family = list(
    hyper = pcprec),
  data = dataf)


## list for control.inla
ctrc <- list(waic = TRUE, dic = TRUE, cpo = TRUE)

## define m for gcv
if(!any(ls()=="level_sets"))
    level_sets <- c(m1 = 1, m3 = 3, m5 = 5, m10 = 10)

## models for u
models <- c("102", "121", "202", "220")
names(models) <- paste0('u', models)

irsel <- c(
    "cpu", "mode", "logfile", 
    paste0("summary.",
           c("fixed", "hyperpar", "random", "fitted.values")),
    "marginals.fixed", "internal.marginals.hyperpar",
    "mlik", "po", "cpo", "dic", "waic", "gcvlist")
