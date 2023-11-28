rm(list=ls())

library(INLAcomp)
library(sp)
library(ggplot2)
library(ggtext)
library(dplyr)
library(viridis)
library(INLA)
library(raster)
library(patchwork)
library(geodata)
library(parallel)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


### --- 0. Previous functions to use --- ####
#' Include a suffix to variable names
#'
#' Used with \code{\link[INLA]{inla.stack}} for building effect stacks.
#'
#' @param x a data.frame
#' @param suffix character string to be appended to variable names
#' @param makeNA logical. If \code{TRUE}, returns the suffixed data.frame filled
#'   with \code{NA}.
#' @export
suffix <- function(x, suffix, makeNA = FALSE) {
  z <- structure(x, names = paste(names(x), suffix, sep = '.'))
  if(makeNA) z[] <- NA
  z
}


### --- function to convert matrix to raster --- ####
rotate <- function(x)(apply(x,2,rev))

matrix_to_raster <- function(m, proj.grid.mat = proj.grid.mat) {
  raster(rotate(t(m)), 
         xmn = min(proj.grid.mat$lattice$loc[,1]), 
         xmx = max(proj.grid.mat$lattice$loc[,1]), 
         ymn = min(proj.grid.mat$lattice$loc[,2]), 
         ymx = max(proj.grid.mat$lattice$loc[,2]), 
         crs = proj4string(INLAcomp::polygon_IP))   
}

# Legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


### 
raster_spatial_effect <- function(spatial, boundary, 
                                  boundary1 = INLAcomp::polygon_IP,
                                  rast = TRUE,
                                  res = c(5000, 5000)) {
  
  if(rast == TRUE){
    res <- res(raster_predict)
  }
  proj <- inla.mesh.projector(mesh,
                              xlim = bbox(boundary)[1,],
                              ylim = bbox(boundary)[2,],
                              dims = c(round((bbox(boundary)[1,"max"] - bbox(boundary)[1,"min"])/res[1],0),
                                       round((bbox(boundary)[2,"max"]-bbox(boundary)[2,"min"])/res[2],0)))
  
  if(length(spatial) == mesh$n){
    field.proj = inla.mesh.project(proj, spatial[1:mesh$n])
    
    field1 <- matrix_to_raster(field.proj, proj.grid.mat = proj)
    field1 <- crop(field1, boundary)
    field1 <- mask(field1, boundary)
    #plot(field1)
  }else{
    res <- matrix(spatial, byrow = FALSE, ncol = k-1)
    total <- apply(res, 2,
                   function(x){
                     res1 <- inla.mesh.project(proj, x)
                     field1 <- matrix_to_raster(res1, proj.grid.mat = proj)
                     field1 <- crop(field1, boundary1)
                     field1 <- mask(field1, boundary1)
                     #plot(field1)
                     #field1
                   })
    field1 <- raster::stack(total)
  }
  field1
}


# Plotting rasters
plot_raster <- function(rast, cat, sc = FALSE, col1 = "G", ymin = 0, ymax = 1) {
  sp_eus2 <- fortify(INLAcomp::polygon_IP_low)
  rast_gc1 <- as(rast[[cat]], "SpatialPixelsDataFrame")
  rast_gc1 <- as.data.frame(rast_gc1)
  colnames(rast_gc1) <- c("value", "x", "y")
  
  gc1_pred <- ggplot() +
    geom_tile(data = rast_gc1, aes(x = x, y = y, fill = value)) +
    geom_polygon(data = INLAcomp::polygon_IP_low, aes(x = long, y = lat, group = group),
                 fill = NA, color = "gray20", size = 0.25) +
    coord_equal(ratio = 1) +
    ggtitle(cat) +
    theme(line = element_blank(), # remove axis lines ..
          axis.text = element_blank(),                       # .. tickmarks..
          axis.title = element_blank(),
          panel.background = element_blank(),
          legend.position = "bottom",
          legend.key.height = unit(0.5, "cm"),
          legend.key.width = unit(4, "cm"),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          plot.title = element_text(size = 17, hjust = 0.5)) +
    scale_fill_viridis(option=col1, direction = -1, limits = c(ymin, ymax)) 
  
  if(sc == TRUE){
    gc1_pred <- gc1_pred +
      ggsn::scalebar(sp_eus2, location = "bottomright",
                     dist = 200, transform = FALSE,
                     st.size = 3, st.dist = 0.05, height=0.03,
                     dist_unit = "km") 
    #ggsn::north(sp_eus2, location = "topright", scale = 0.2,
    #            symbol = 4)
    #ggsn::north2(sp_eus2, x = 0.65, y = 0.9, scale = 0.1, symbol = 1)
    # 
  }
  gc1_pred
}




### --- 0.1. Plotting groups --- ####
### Plotting neighbours
plot_neigh <- function(i, cv_mod, 
                       title = "Prior",
                       friends = friends_list,
                       polygon1 = NULL) {
  # cv_mod$groups[[i]]$idx %>% 
  #   cut(., c(1, 301, 602, 903)) %>%
  #   factor(labels = c("alr.gc1", "alr.gc2", "alr.gc3")) %>%
  #   data.frame(id = cv_mod$groups[[i]]$idx, .) -> groups
  
  index <- numeric()
  for (j in 1:length(friends[[i]])){
    index <- c(index, cv_mod$groups[[friends[[i]][j]]]$idx %in% ((n*(j-1)):(n*j)) %>%
      cv_mod$groups[[friends[[i]][j]]]$idx[.])
  }
  
  data_ext %>% .[index,] -> data_group
  
  data_ext[friends_list[[i]],]  ->  data_group1
  
  if(is.null(polygon1)==TRUE){
    data_group %>%
      ggplot(data = .) +
      geom_polygon(data = penin_2_fort, aes(x = long,
                                            y = lat,
                                            group = group),
                   colour = 'gray80', fill = 'gray100') +
      geom_point(data = data_ext,
                 #dplyr::filter(name == "alr.gc1"),
                 aes(x = x, y = y, group = alr_gc),
                 col = "gray80", size = 0.5) +
      tidyterra::geom_spatvector(data = ip_prov,
                                 colour = 'gray40', fill = 'gray100') +
      geom_point(aes(x = x,
                     y = y),
                 alpha = 1, 
                 size = 1.3,
                 shape = 3) +
      #           col = "blue", shape = 4) + 
      scale_color_viridis(option="viridis", direction = -1,
                          name = "Cor", limits = c(-1,1)) +
      geom_point(data =   data_group1,
                 aes(x = x,
                     y = y,
                     group = alr_gc), 
                 col = "blue", size = 2)+ 
      theme_bw() +
      theme(line = element_blank(),                          # remove axis lines ..
            axis.text = element_blank(),                       # .. tickmarks..
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            legend.position = "bottom",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15)) +
      #guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
      xlab("") +
      ylab("") +
      facet_wrap( ~ alr_gc ) -> p1
    
  }else{
        limits <- terra::ext(polygon1)[c(1:4)]
    data_group %>%
      ggplot(data = .) +
      tidyterra::geom_spatvector(data = polygon1,
                                 colour = 'gray40', fill = 'gray100') +
      geom_point(data = data_ext,
                 #dplyr::filter(name == "alr.gc1"),
                 aes(x = x, y = y, group = alr_gc),
                 col = "gray80", size = 0.5) +
      geom_point(aes(x = x,
                     y = y),
                 alpha = 1, shape = 3, size = 1.3) + 
      geom_point(data =   data_group1,
                 aes(x = x,
                     y = y,
                     group = alr_gc), 
                 col = "blue", size = 2) +
      xlim(limits[c(1,2)]) +
      ylim(limits[c(3,4)]) +
      theme_bw() +
      theme(line = element_blank(),                          # remove axis lines ..
            axis.text = element_blank(),                       # .. tickmarks..
            axis.title = element_blank(),
            plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            legend.position = "bottom",
            legend.key.height = unit(0.5, "cm"),
            legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 17),
            legend.text = element_text(size = 15),
            strip.text = element_text(size = 15)) +
      #guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
      xlab("") +
      ylab("") +
      facet_wrap( ~ alr_gc ) -> p1
  }
  p1
}


### --- 1. Reading the data --- ####
data <- INLAcomp::arabidopsis
head(data)

penin_2 <- INLAcomp::polygon_IP
plot(penin_2)

# We scale the covariates
bio_scaled <- data %>%
  dplyr::select(starts_with("bio")) %>%
  scale(.)

colnames(bio_scaled) <- paste0("sc_", colnames(bio_scaled ))
data <- cbind(data, bio_scaled)
data[, c("x", "y", paste0("sc_bio", c(1, 2, 3, 4, 8, 12, 15, 18)))]

# Aggregated to 10x10km by mean
bio1 <- raster("rasters/bio1.tif")
bio12 <- raster("rasters/bio12.tif")

# Scaled it
attr(bio_scaled,"scaled:center")["bio12"]
attr(bio_scaled,"scaled:scale")["bio12"]

sc_bio1 <- (bio1 - attr(bio_scaled,"scaled:center")["bio1"])/attr(bio_scaled,"scaled:scale")["bio1"]
sc_bio12 <- (bio12 - attr(bio_scaled,"scaled:center")["bio12"])/attr(bio_scaled,"scaled:scale")["bio12"]

par(mfrow=c(1,2))
plot(sc_bio1)
plot(sc_bio12)


## Boundaries: Spain
penin2_com <- geodata::gadm(country = "ESP", level = 1, path = getwd())
penin2_com <- project(penin2_com, "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")
penin2_com <- penin2_com[-which(penin2_com$NAME_1 == "Islas Baleares" | penin2_com$NAME_1 == "Islas Canarias"), ]

# Boundaries: Portugal
port_com <- geodata::gadm(country = "Portugal", level = 1, path = getwd())
port_com <- project(port_com, "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")

ip_prov <- rbind(penin2_com, port_com)
ip_prov <- crop(ip_prov, INLAcomp::polygon_IP)
penin_2_fort <- fortify(penin_2)


### --- 2. Descriptives --- ####
data_plot <- data %>%
  dplyr::select(x, y, gc1, gc2, gc3, gc4) %>%
  tidyr::pivot_longer(all_of(paste0("gc", 1:4)))


# Descriptive plot. Original compositional data
p1 <- ggplot(data_plot) +
  tidyterra::geom_spatvector(data = ip_prov,
                             colour = 'gray50',
                             fill = NA,
                             linewidth = 0.2) +
  geom_point(data = data_plot,
             aes(x = x,
                 y = y,
                 #size = value,
                 color = value),
             alpha = 1, size = 1) +
  scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(7, "Spectral"))) +
  #coord_fixed(ratio = 1) +
  facet_wrap(~ name, ncol = 4) +
  theme_bw() +
  theme(line = element_blank(),                          # remove axis lines ..
        axis.text = element_blank(),                       # .. tickmarks..
        axis.title = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  #guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  xlab("") +
  ylab("") 

print(p1)

# Selecting the reference category. The one with less variance
data %>% dplyr::select(gc1, gc2, gc3, gc4) %>%
  apply(., 2, function(x){var(log(x))})
# Also that they don't have so small values
data %>% 
  tidyr::pivot_longer(., cols = all_of(paste0("gc", 1:4)), 
                      names_to  = "gc",
                      values_to = "val_gc") %>%
  ggplot(data = .) +
  geom_density(aes(x = val_gc, fill = gc), alpha = 0.4) +
  theme_bw() +
  theme(legend.position = "bottom") +
  xlab("Membership probability to GC") +
  ylab("Density")

### Reference category is gc4
n <- dim(data)[1]
k <- 4
result_order <- c(paste0("gc", 1:k))
ref_cat <- paste0("gc", k)

names_gc <- c(result_order[!(result_order %in% ref_cat)])
result_order <- c(names_gc, ref_cat)
data <- cbind(data, alr=compositions::alr(data[, result_order]))

data[, c("x", "y", paste0("sc_bio", c(1, 2, 3, 4, 8, 12, 15, 18)), paste0("alr.gc", 1:3))]


# Descriptive map alr --- #
data_plot_alr <- data %>%
  dplyr::select(all_of(c("x", "y", paste0("alr.", names_gc)))) %>%
  tidyr::pivot_longer(all_of(paste0("alr.", names_gc)))

# Plotting the map
p2 <- ggplot(data_plot_alr) +
  tidyterra::geom_spatvector(data = ip_prov,
                             colour = 'gray50',
                             fill = NA,
                             linewidth = 0.2) +
  geom_point(data = data_plot_alr,
             aes(x = x,
                 y = y,
                 #size = value,
                 color = value),
             alpha = 1, size = 1) +
  #scale_color_viridis(option = "A", direction = -1)
  scale_color_gradientn(colours = (RColorBrewer::brewer.pal(7, "YlOrBr"))) +
  facet_wrap(~ name) +
  theme_bw() +
  theme(line = element_blank(),                          # remove axis lines ..
        axis.text = element_blank(),                       # .. tickmarks..
        axis.title = element_blank(),
        panel.background = element_blank(),
        legend.position = "bottom",
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        strip.text = element_text(size = 15)) +
  #guides(fill = guide_legend(title.position="top", title.hjust = 0.5)) +
  xlab("") +
  ylab("") 

graphics.off()
pdf("Fig7_CoDa_data.pdf",
   width =9, height = 7)
(p1)/(p2) +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(face = 'bold', size = 17))
dev.off()

graphics.off()
png("Fig7_CoDa_data.png",
   width = 1400, height = 900, res = 100)
(p1)/(p2) +
  plot_annotation(tag_levels = list(c('A', 'B'))) &
  theme(plot.tag = element_text(face = 'bold', size = 17))
dev.off()


### --- 3. Preparing the data for dealing with R-INLA --- ####
data_ext <- data %>% 
  tidyr::pivot_longer(., cols = all_of(paste0("alr.", names_gc)), 
                      names_to  = "alr_gc",
                      values_to = "val_alr_gc") %>%
  .[order(ordered(.$alr_gc)),]

head(data_ext[, c("x", "y", "alr_gc", "val_alr_gc")])


## -----------------------------------------------
data_ext$alr_gc <- ordered(data_ext$alr_gc)
k.group <- data_ext$alr_gc %>% as.numeric() #For group
k.repl  <- data_ext$alr_gc %>% as.numeric() #For replication
head(data.frame(k.group, k.repl))


# Boundaries and mesh
## -----------------------------------------------
penin_inla <- inla.sp2segment(INLAcomp::polygon_IP_low)
mesh <- inla.mesh.2d(boundary = penin_inla,
                     max.edge = c(40000, 100000),
                     cutoff   = 30000,
                     min.angle= 30,
                     offset   = c(10000, 300000))


## -----------------------------------------------
plot(penin_2, xlim = c(mesh$loc[,1] %>% range()),
     ylim = c(mesh$loc[,2] %>% range(.)))
plot(mesh, add = TRUE)
points(data[,c("x", "y")], pch = 20)


## SPDE
## -----------------------------------------------
size <- min(c(diff(range(penin_inla$loc[,1])), diff(range(penin_inla$loc[,2]))))
range0 <- size / 4 	# ~ default
spde <- inla.spde2.pcmatern(mesh        = mesh, 
                            prior.range = c(range0, 0.25), # P(range < range0) = 0.25
                            prior.sigma = c(1, 0.01)) # P(sigma > 1) = 0.01

## iNDEX for SPDE and CoDa
## -----------------------------------------------
iset <- inla.spde.make.index('i', n.spde = spde$n.spde,
                             n.repl  = k-1) #Replicating spatial effect


## -----------------------------------------------
iset2 <- inla.spde.make.index('i', n.spde = spde$n.spde)
iset_aux <- matrix(NA, ncol = k-1, nrow = k-1)
diag(iset_aux) <- rep(1, k-1)
iset2 <- kronecker(iset_aux, iset2$i)
colnames(iset2) <- paste0("iset.alr", 1:(k-1))

## Matrix A
## -----------------------------------------------
A.est <- inla.spde.make.A(mesh  = mesh,
                          loc   = cbind(data_ext$x, data_ext$y), 
                          repl  = k.repl) 
A.est2 <- A.est

## Index for the shared random effect
## -----------------------------------------------
id.z <- rep(1:(dim(data_ext)[1]/(k-1)), k-1)
id.z

## -----------------------------------------------
id.cat <- rep(1:(k-1), rep(n, k-1))


#Response
## -----------------------------------------------
names_alr <- paste0("alr.", names_gc)


1:length(names_alr) %>%
  lapply(., function(i){
    data_ext %>% 
      dplyr::filter(alr_gc == names_alr[i]) -> data_comp_i
    #Response
    y_alr <- matrix(ncol = names_alr %>% length(.), nrow = dim(data_comp_i)[1])
    y_alr[, i] <- data_comp_i$val_alr_gc
  }) -> y_alr

1:length(names_alr) %>%
  lapply(., function(i){
    y_aux <- data_ext %>%
      dplyr::select(val_alr_gc, alr_gc) %>%
      dplyr::filter(alr_gc == names_alr[i]) %>%
      dplyr::select(val_alr_gc) %>%
      as.matrix(.)
    aux_vec <- rep(NA, k-1)
    aux_vec[i] <- 1
    kronecker(aux_vec, y_aux)
  }) -> y_alr_list

y_alr <- do.call(cbind, y_alr_list)
y_alr[1:10,]


# Covariates
## -----------------------------------------------
variables <- c("intercept", data %>% 
                 dplyr::select(starts_with("sc_")) %>% 
                 colnames(.))
id.names <- paste0("id.", variables)
id.variables <- rep(id.cat, length(variables)) %>% 
  matrix(., ncol = length(variables), byrow = FALSE)
colnames(id.variables) <- id.names

# Inla.stack for estimation
## -----------------------------------------------
stk.est <- inla.stack(data    = list(resp = y_alr),
                      A       = list(A.est, 1, A.est2),
                      effects = list(c(iset),
                                     cbind(data_ext %>% 
                                             dplyr::select(starts_with("sc_bio")),
                                           id.z,
                                           id.variables, 
                                           intercept = 1),
                                     data.frame(iset2)),
                      tag     = 'est')


### --- 4. Fitting the differente models --- ####
### ----- 4.1. Type I --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeI <- resp ~ -1 + 
  intercept + 
  sc_bio1 +
  sc_bio12 +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeI <- inla(formula.typeI,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeI)

### ----- 4.2. Type II --- ####
formula.typeII <- resp ~ -1 + 
  f(id.intercept, intercept,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio1, sc_bio1, #BIO1
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio12, sc_bio12, #bio12
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeII <- inla(formula.typeII,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeII)


### ----- 4.3. Type III --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeIII <- resp ~ -1 + 
  intercept + 
  sc_bio1 +
  sc_bio12 +
  f(iset.alr1, model = spde) +
  f(iset.alr2, copy = "iset.alr1", 
    fixed = TRUE) +
  f(iset.alr3, copy = "iset.alr1", 
    fixed = TRUE) +  
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeIII <- inla(formula.typeIII,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeIII)


### ----- 4.4. Type IV --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeIV <- resp ~ -1 + 
  f(id.intercept, intercept,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio1, sc_bio1, #BIO1
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio12, sc_bio12, #bio12
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(iset.alr1, model = spde) +
  f(iset.alr2, copy = "iset.alr1", 
    fixed = TRUE) +
  f(iset.alr3, copy = "iset.alr1", 
    fixed = TRUE) + 
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))),
    constr = TRUE)

model.typeIV <- inla(formula.typeIV,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeIV)


### ----- 4.5. Type IV --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeV <- resp ~ -1 + 
  intercept + 
  sc_bio1 +
  sc_bio12 +
  f(iset.alr1, model = spde) +
  f(iset.alr2, copy = "iset.alr1", 
    fixed = FALSE) +
  f(iset.alr3, copy = "iset.alr1", 
    fixed = FALSE) +  
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeV <- inla(formula.typeV,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeV)


### ----- 4.6. Type VI --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeVI <- resp ~ -1 + 
    f(id.intercept, intercept,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE,
    constr = TRUE) +
  f(id.sc_bio1, sc_bio1, #bio1
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio12, sc_bio12, #bio12
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(iset.alr1, model = spde) +
  f(iset.alr2, copy = "iset.alr1",
    hyper = list(beta = list(fixed = FALSE, params = c(0, 1)))) +
  f(iset.alr3, copy = "iset.alr1",
         hyper = list(beta = list(fixed = FALSE, params = c(0, 1)))) +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))),
    constr = TRUE)

model.typeVI <- inla(formula.typeVI,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)
summary(model.typeVI)


### ----- 4.7. Model type VII --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeVII <- resp ~ -1 + 
  intercept + 
  sc_bio1 +
  sc_bio12 +
  f(i,
    model = spde,
    replicate = i.repl) +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeVII <- inla(formula.typeVII,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeVII)


### ----- 4.8. Model type VIII --- ####
list_prior <- rep(list(list(prior = "pc.prec", param = c(1, 0.01))), k-1)

formula.typeVIII <- resp ~ -1 + 
  f(id.intercept, intercept,
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio1, sc_bio1, #BIO1
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(id.sc_bio12, sc_bio12, #bio12
    model   = "iid",
    initial = log(1/1000),
    fixed   = TRUE) +
  f(i,
    model = spde,
    replicate = i.repl) +
  f(id.z,
    model = "iid",
    hyper = list(prec = list(prior = "pc.prec",
                             param = c(1, 0.01))), constr = TRUE)

model.typeVIII <- inla(formula.typeVIII,
                   family         = rep("gaussian", k - 1),
                   data           = inla.stack.data(stk.est),
                   control.compute = list(config = TRUE, 
                                          dic  = TRUE,
                                          waic = TRUE,
                                          cpo  = TRUE),
                   control.predictor = list(A = inla.stack.A(stk.est), 
                                            compute = TRUE),
                   control.family = list_prior,
                   verbose = FALSE)

summary(model.typeVIII)


### --- 5. Computing DIC, WAIC and predictive measures --- ####
model_list <- list("typeI"=model.typeI,
                   "typeII"=model.typeII,
                   "typeIII"=model.typeIII,
                   "typeIV"=model.typeIV,
                   "typeV"=model.typeV,
                   "typeVI"=model.typeVI,
                   "typeVII"=model.typeVII,
                   "typeVIII"=model.typeVIII)

# We select as reference model typeVIII
measures_fit <- pbapply::pblapply(model_list, function(mod1){
  xx <- inla.posterior.sample(1000, mod1)
  inf <- lapply(xx, INLAcomp::extract_lp_sigma)

  # DIC
  dic.mod1 <- INLAcomp::dic.mult(inf, y = data[, c(paste0("alr.gc", 1:(k-1)))])

  # WAIc
  waic.mod1 <- INLAcomp::waic.mult(inf, y = data[, c(paste0("alr.gc", 1:(k-1)))])
  
  cbind(dic.mod1, waic.mod1)
})
measures_fit <-  as.data.frame(do.call(rbind, measures_fit))

# Neighbours
friends_list <- 1:903 %>%
  lapply(., function(x){
    c(seq(x, 903, by = 301)[-1],
      rev(seq(x, 1, by = -301))) -> res
    res[order(res)]
  })

model_list_group <- model_list$typeVIII %>% list(.)

# Logarithmic Score (LS)
pbapply::pblapply(model_list_group, function(model.group){
  
  ## We create the groups with a general posterior group
  group.cpo <- INLA::inla.group.cv(result = model.group,
                                   num.level.sets = -1,
                                   strategy = "posterior",
                                   friends = friends_list)
  
  group.cpo3 <- INLA::inla.group.cv(result = model.group,
                                    num.level.sets = 3,
                                    strategy = "posterior",
                                    friends = friends_list)
  
  group.cpo5 <- INLA::inla.group.cv(result = model.group,
                                    num.level.sets = 5,
                                    strategy = "posterior",
                                    friends = friends_list)
  
  group.cpo10 <- INLA::inla.group.cv(result = model.group,
                                    num.level.sets = 10,
                                    strategy = "posterior",
                                    friends = friends_list)
  
  measures_cpo <- lapply(model_list, function(mod1){
    mod1.group.cv1 <- inla.group.cv(mod1, group.cv = group.cpo)
    mod1.group.cv3 <- inla.group.cv(mod1, group.cv = group.cpo3)
    mod1.group.cv5 <- inla.group.cv(mod1, group.cv = group.cpo5)
    mod1.group.cv10 <- inla.group.cv(mod1, group.cv = group.cpo10)
    
    measures <- data.frame(LS=-mean(log(mod1.group.cv1$cv)),
                           LS3=-mean(log(mod1.group.cv3$cv)),
                           LS5=-mean(log(mod1.group.cv5$cv)),
                           LS10=-mean(log(mod1.group.cv10$cv)))
    measures
  })
  
  measures_cpo <-  as.data.frame(do.call(rbind, measures_cpo))
  measures_cpo
}) -> measures_cpo


# Joining DIC, WAIC and LS
lapply(measures_cpo, function(xx){
  cbind(data.frame(DIC = measures_fit$dic,
                   WAIC = measures_fit$waic), 
        xx)
}) -> measures

round(measures[[1]],3)

measures[[1]] %>%
  dplyr::select(DIC, WAIC, LS, LS3, LS5, LS10)  %>%
  xtable::xtable(., digits = 3)


### --- 6. Plotting groups selected for comparing models --- ####
i1 <- 88

# Selecting a part around
epsilon  <- 8*10^4
epsilon <- 8*10^4
minext <- data_ext[i1, c("x", "y")] %>% as.numeric(.) - epsilon 
maxext <- data_ext[i1, c("x", "y")] %>% as.numeric(.) + epsilon
extent <- ext(minext[1], maxext[1], minext[2], maxext[2])
ip_i1 <- crop(ip_prov, extent)
plot(ip_i1)

model.typeVIII$group.cv10 <- inla.group.cv(model.typeVIII, 
                                           num.level.sets = 10, 
                                           friends = friends_list,
                                           strategy = "posterior")

graphics.off()
pdf("Fig8_CoDa_groups_m10.pdf", width = 10, height = 4)
plot_neigh(i = i1,
           cv_mod = model.typeVIII$group.cv10,
           friends = friends_list,
           title = i,
           polygon1 = ip_i1)
dev.off()
