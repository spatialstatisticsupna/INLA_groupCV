library(parallel)
library(ggplot2)
library(patchwork)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(INLA)

## figures folder 
figures.dir <- here::here("LaTeX", "Figures/")
figures.path <- paste0(figures.dir, "Spacetime_")
figures.path

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

## load the mesh
smesh <- readRDS("rdsfiles/smesh_sresol30.rds")

datalocs <- unique(cbind(dataf$xloc, dataf$yloc))
(nlocs <- nrow(datalocs))

(ndata <- nrow(dataf))
(nt <- max(dataf$time))

bb <- matrix(st_bbox(map_utm), 2)
bb

rxy <- apply(bb, 1, diff)
rxy

models <- c(102, 121, 202, 220)

rfls <- paste0(
    "rdsfiles/",
    "fg_u", models, 
    "_n", ndata,
    "_ns", smesh$n,
    "_nt", nt,
    ".rds")
rfls

if(any(!file.exists(rfls)))
    stop("please run the 'fit4models.R' code!")

nm <- length(models)
results <- vector("list", nm)
names(results) <- paste0("u", models)
names(results)

for (k in 1:nm) {
### file with the result
    rfl <- rfls[k]
    if(file.exists(rfl)) {
        cat("Loading the previously saved result ... ")
        cat(rfl, "\n")
        results[[k]] <- readRDS(rfl)
    } 
} ## end for stmodels

sapply(results, function(r) r$cpu[["Total"]]) / 60
sapply(results, function(r) r$cpu[[length(r$cpu)]])

sapply(results, function(r) r$cpu[["Total"]]) / 
    sapply(results, function(r) r$cpu[[length(r$cpu)]])

sapply(results, function(r) unname(r$mode$theta))

(ab0.s <- range(sapply(results, function(r)
    range(r$summary.random$spatial$mean))))

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0), las = 1)
for(k in 1:nm) {
    hist(results[[k]]$summary.random$spatial$mean,
         round((ab0.s[1]-1)*10):round((ab0.s[2]+1)*10)/10, 
         xlab = '', ylab = '', main = '')
}

(ab0.st <- range(sapply(results, function(r)
    range(r$summary.random$field$mean))))

par(mfrow = c(2, 2), mar = c(3, 3, 0.5, 0.5), mgp = c(1.5, 0.5, 0), las = 1)
for(k in 1:nm) {
    hist(results[[k]]$summary.random$field$mean,
         round((ab0.st[1]-1)*10):round((ab0.st[2]+1)*10)/10,
         xlab = '', ylab = '', main = '')
}

sapply(results, function(r) mean(r$summary.random$field$mean))

ii.est <- which(!is.na(dataf$fg))
rstats <- data.frame(t(rbind(
    sapply(results, function(r)
        c(INLAspacetime::stats.inla(r, i = ii.est, y = dataf$fg,
                     fsummarize = sum)[c(5:6, 1:2, 7:8, 3:4)],
          sapply(r$gcvlist, function(r)
              -sum(log(r$cv[ii.est]))))))))
rstats$mse <- sqrt(rstats$mse)
names(rstats) <- toupper(gsub("mse", "smse", names(rstats)))
rownames(rstats) <- LETTERS[1:4]

print(round(rstats, 1))

xtable::xtable(rstats[, c(3, 10:12)])

nstt <- length(jjs <- 1:12)

snams <- toupper(colnames(rstats)[jjs])
snams
snams[snams=="SMSE"] <- 'sqrt(MSE)'
snams[snams=="M1"] <- "LGCPO, m=1"
snams[snams=="M3"] <- "LGCPO, m=3"
snams[snams=="M5"] <- "LGCPO, m=5"
snams[snams=="M10"] <- "LGCPO, m=10"
snams


int2user <- function(r) {
    r.m <- r$internal.marginals.hyperpar
    m <- list(e.sigma = inla.tmarginal(function(x) exp(-x/2), r.m[[1]]))
##    m$v.sigma <- inla.tmarginal(function(x) ### set it as second to show!
  ##      exp(-x/2), r.m[[5]])
    m$s.range <- inla.tmarginal(exp, r.m[[2]])
    m$s.sigma <- inla.tmarginal(exp, r.m[[3]])
    m$u.Rspace <- inla.tmarginal(exp, r.m[[4]]) 
    m$u.Rtime <- inla.tmarginal(exp, r.m[[5]])
    m$u.sigma <- inla.tmarginal(exp, r.m[[6]])
    return(m)
}

hmarginals <- lapply(results, int2user)
hmsummary <- lapply(hmarginals, sapply, function(m)
    unlist(inla.zmarginal(m, silent = TRUE)))
hmsummary[[1]]

sapply(results, function(x) x$summary.hy$mean[1])

jjpar <- 1:6; names(jjpar) <- names(hmarginals[[1]])
jjpar

mhpars <- sapply(jjpar, function(j) {
    sapply(hmsummary, function(mm) mm[1, j])
})

sapply(results, function(r)
    r$summary.fix$mean[2] / r$summary.fix$sd[2])

t(sapply(results, function(r) r$summary.fix[1, c(1, 2, 3, 5)]))
t(sapply(results, function(r) r$summary.fix[2, c(1, 2, 3, 5)]))

100 * cor(dataf$fg, 
          sapply(results, function(r) r$summary.fitted.values$mean[1:nrow(dataf)]),
          use = 'pair')

if(FALSE) {
   
    par(mfrow = c(2, 2), mar = c(2, 2, 0, 0), mgp = c(1.5, 0.5, 0), las = 1)
    for(i in 1:k) {
        plot(results[[k]]$summary.fitted.values$mean[1:nrow(dataf)] ,
             dataf$fg, xlab = '', ylab = '', bty = 'n', asp = 1, pch = 19)
        abline(0:1, col = 2)
    }

}

bb

sgrid <- list(
    x = seq(bb[1], bb[3], 2),
    y = seq(bb[2], bb[4], 2))
sgrid$nx <- length(sgrid$x)
sgrid$ny <- length(sgrid$y)
sgrid$loc <- as.matrix(expand.grid(sgrid[c("x", "y")]))
str(sgrid)

sgrid$over <- !is.na(over(
    SpatialPoints(sgrid$loc, proj4string = CRS(utm_crs)), 
    as_Spatial(map_utm))[, 1])

sgrid$proj <- inla.mesh.projector(
    mesh = smesh,
    loc = sgrid$loc[sgrid$over, ]
)

library(terra)
library(tidyterra)

ii.rsel <- 2

round(results[[ii.rsel]]$summary.hy, 2)

round(htab <- hmsummary[[ii.rsel]][c(1, 2), ], 2)
colnames(htab) <- paste0(
    "$",
    c(
        "\\sigma_e",
        "\\phi_t",
        "\\phi_s",
        "\\sigma_u",
        "\\rho_s",
        "\\sigma_s"
    ), "$")

print(
  xtable::xtable(htab, digits = 2, align = "lrrrrrr"),
  sanitize.text.function = function(x) x
)

zz <- matrix(NA, sgrid$nx, sgrid$ny)

zz[sgrid$over] <- as.vector(
    inla.mesh.project(
        sgrid$proj,
        results[[ii.rsel]]$summary.random$spatial$mean))

ab.s <- max(abs(range(zz[sgrid$over]))) * c(-1, 1)
ab.s

z <- rast(data.frame(
    x = rep(sgrid$x, times = sgrid$ny),
    y = rep(sgrid$y, each = sgrid$nx),
    u = as.vector(zz)),
    crs = utm_crs
    )

gg.smean <- ggplot() + theme_minimal() + xlab("") + ylab("") +
    geom_spatraster(data = z) +
    scale_fill_distiller(
        type = "div",
        palette = "RdBu",
        na.value = "transparent",
        limits = ab.s 
    ) +
    geom_sf(data = map_utm, fill = "transparent") +
    theme(legend.position = c(0.9, 0.8),
          legend.direction="vertical", ##horizontal",  
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) 

smean.fl <- paste0(figures.path, "smean.png")
smean.fl

png(smean.fl,
    pointsize = 20,
    width = 11,
    heigh = 20,
    units = "in",
    res = 300,
    type = "cairo"
    )

gg.smean

dev.off()

if(FALSE)
    system(paste("eog", smean.fl, "&"))

nt
length(tt1sel <- seq(1, nt, 1))
tt1sel

stprojected <- inla.mesh.project(
    sgrid$proj,
    matrix(results[[ii.rsel]]$summary.random$field$mean,
           smesh$n))

if(FALSE) {

    dim(stprojected)

    l0 <- locator()

    iig <- sapply(1:length(l0$x), function(i)
        which.min((l0$x[i] - sgrid$loc[sgrid$over, 1])^2 +
                  (l0$y[i] - sgrid$loc[sgrid$over, 2])^2))
    

    par(mfrow=c(1,1), mar = c(.1, .1, .1, .1))
    plot(SpatialPoints(datalocs), type = "n")
    plot(sgrid$loc[which(sgrid$over)[iig], ], asp = 1,
         xlim = bb[1, ], ylim = bb[2, ])
    INLAspacetime::stlines(t(stprojected[iig, ]),
                           SpatialPoints(sgrid$loc[which(sgrid$over)[iig], ]),
                           yscale = 0.25)
    
}

ab.st1 <- max(abs(range(stprojected))) * c(-1, 1)
ab.st1

Date0 <- as.Date("2021-06-30") + 1:nt
range(Date0)

stfigs <- vector("list", length(tt1sel))
k <- 0
for(tt in tt1sel) {
    k <- k + 1
    zz[sgrid$over] <- as.vector(stprojected[, tt])
    zt <- rast(data.frame(
        x = rep(sgrid$x, times = sgrid$ny),
        y = rep(sgrid$y, each = sgrid$nx),
        u = as.vector(zz)),
        crs = utm_crs
        )    
    stfigs[[k]] <-
        ggplot() + theme_minimal() + xlab("") + ylab("") +
        geom_spatraster(data = zt) +
        scale_fill_distiller(
            type = "div",
            palette = "RdBu",
            na.value = "transparent",
            limits = ab.st1
        ) +        
        theme(legend.position = c(0.9, 0.8),
              legend.direction="vertical", ##horizontal",  
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.ticks.y = element_blank()) + 
        labs(fill = format(Date0[tt], "%b, %d")) + 
        geom_sf(data = map_utm, fill = "transparent")
}

print(wrap_plots(stfigs[1:4], ncol = 4))

utt1.fl <- paste0(figures.path, "u1mean.png")

png(utt1.fl,
  pointsize = 10,
  width = 14,
  heigh = 20,
  units = "in",
  res = 100,
  type = "cairo"
)

print(wrap_plots(stfigs, ncol = 8))

dev.off()

