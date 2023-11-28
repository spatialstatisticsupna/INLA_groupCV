rm(list=ls())
library(INLA)  # R-INLA stable version 23.09.09
library(RColorBrewer)
library(sf)
library(splines)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#########################
## Auxiliary functions ##
#########################
DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance,
             p.eff=x$dic$p.eff,
             DIC=x$dic$dic,
             WAIC=x$waic$waic)
}

CV <- function(x,O,E){
  LS <- -mean(log(x$cv))
  
  # Computation of E[Y_i | y_{I_i}] using trapezoid numerical integration #
  expectation = numeric(n)
  for (i in 1:n){
    mu = x$mean[i]
    sd = x$sd[i]
    xx = seq(mu-6*sd,mu+6*sd,length.out = 100)
    yy = E[i]*exp(xx)*dnorm(xx, mean=mu, sd=sd)
    expectation[i] = pracma::trapz(xx,yy)
  }
  MSPE <- mean((expectation-O)^2)
  
  data.frame(LS=LS, MSPE=MSPE)
}

marginal.basis <- function(x, p=3, q=10){

  x <- (x-min(x))/(max(x)-min(x))
  dist <- (max(x)-min(x))/q
  xl <- min(x)-dist*0.05
  xr <- max(x)+dist*0.05
  dx <- (xr-xl)/q
  knots <- seq(xl-p*dx, xr+p*dx, by=dx)
  
  B <- spline.des(knots,x,p+1)$design
  
  return(B)
}

Rten <- function(X1,X2){ ## Row-wise Kronecker product 
  one1 <- matrix(1,1,ncol(X1))
  one2 <- matrix(1,1,ncol(X2))
  kronecker(X1,one2)*kronecker(one1,X2)
}

sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"

lunif = "expression:
          a = 1;
          b = 1;
          beta = exp(theta)/(1+exp(theta));
          logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
          log_jacobian = log(beta*(1-beta));
          return(logdens+log_jacobian)"


###########################################
## Load dowry mortality data (year 2011) ##
###########################################
load("data_UP.Rdata")

data_UP$X1 <- as.numeric(scale(data_UP$x1)) # sex ratio
data_UP$X2 <- as.numeric(scale(data_UP$x2)) # per capita income 
data_UP$X3 <- as.numeric(scale(data_UP$x3)) # number of murders per 100000 inhabitants 

paleta <- brewer.pal(6,"Blues")
values <- c(-Inf,seq(-2,2),Inf)

tmap_mode("plot")
Mapa <- tm_shape(data_UP) + 
  tm_polygons(col=c("SMR",c("X1","X2","X3")), palette=paleta, title="", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position=0.2, legend.text.size=1, panel.label.size=1.5,
            panel.labels=c("SMR","Sex ratio","Per Capita Income","Murder rate"),
            panel.label.bg.color="lightskyblue", legend.outside=T, legend.outside.position="right", legend.frame=F) + 
  tm_facets(nrow=1, ncol=4)
tmap_save(Mapa, file="SpatialData_MapCovariates.pdf")


##############################
## Fit the models with INLA ##
##############################
n <- nrow(W)
Rs <- Diagonal(n,colSums(W))-W

data.INLA <- data.frame(O=data_UP$O, E=data_UP$E, ID.area=seq(1,n),
                        X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3)

## M1: Intercept + fixed effects 
#################################
formula1 <- O ~ 1 + X1 + X2 + X3

M1 <- inla(formula1, family="poisson", data=data.INLA, E=E,
           control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
           control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
summary(M1)


## M2: CAR model
#################
formula2 = O ~ 1 + 
  f(ID.area, model='bym2', graph=Rs, constr=TRUE,
    hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif)))

# formula2 = O ~ 1 + 
#   f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE,
#     hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),
#                phi=list(prior='pc', param=c(0.5,0.5))))

M2 <- inla(formula2, family="poisson", data=data.INLA, E=E,
           control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
           control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
summary(M2)


## M3: CAR + covariates 
########################
formula3 = O ~ 1 + X1 + X2 + X3 + 
  f(ID.area, model='bym2', graph=Rs, constr=TRUE,
    hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif)))

# formula3 = O ~ 1 + X1 + X2 + X3 +
#   f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE,
#     hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),
#                phi=list(prior='pc', param=c(0.5,0.5))))

M3 <- inla(formula3, family="poisson", data=data.INLA, E=E,
           control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
           control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
summary(M3)


## M4: P-spline model 
######################
B1 <- marginal.basis(x=coordinates(as(data_UP,"Spatial"))[,1], p=3, q=10)
k1 <- ncol(B1)

B2 <- marginal.basis(x=coordinates(as(data_UP,"Spatial"))[,2], p=3, q=10)
k2 <- ncol(B2)

Bs <- Rten(B2,B1)
ks <- ncol(Bs)

D1 <- diff(diag(k1),differences=1)
D2 <- diff(diag(k2),differences=1)
P1 <- t(D1)%*%D1
P2 <- t(D2)%*%D2

data.INLA <- list(O=data_UP$O, E=data_UP$E,
                  Intercept=c(1,rep(NA,ks)),
                  idx=c(NA,1:ks))

Ps <- list(as(kronecker(diag(k2),P1),"dgCMatrix"),
           as(kronecker(P2,diag(k1)),"dgCMatrix"))

formula4 <- O ~ -1 + Intercept +
  f(idx, model="generic3", Cmatrix=Ps, constr=TRUE, diagonal=1e-6,
    hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)))
  
# formula4 <- O ~ -1 + Intercept +
#   f(idx, model="generic0", Cmatrix=Ps, rankdef=1, constr=TRUE, hyper=list(prec=list(prior=sdunif)))
# 
# Ps <- kronecker(diag(k2),P1) + kronecker(P2,diag(k1))

M4 <- inla(formula4, family="poisson", data=data.INLA, E=E,
            control.predictor=list(compute=TRUE, A=cbind(rep(1,n), Bs), link=1, cdf=c(log(1))),
            control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
summary(M4)


## M5: P-spline + covariates 
#############################
data.INLA <- list(O=data_UP$O, E=data_UP$E, X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3,
                  Intercept=c(1,rep(NA,3+ks)),
                  beta1=c(rep(NA,1),1,rep(NA,2+ks)),
                  beta2=c(rep(NA,2),1,rep(NA,1+ks)),
                  beta3=c(rep(NA,3),1,rep(NA,ks)),
                  idx=c(rep(NA,4),1:ks))

formula5 <- O ~ -1 + Intercept + beta1 + beta2 + beta3 +
  f(idx, model="generic3", Cmatrix=Ps, constr=TRUE, diagonal=1e-6,
    hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)))

M5 <- inla(formula5, family="poisson", data=data.INLA, E=E,
            control.predictor=list(compute=TRUE, A=cbind(rep(1,n), X1, X2, X3, Bs), link=1, cdf=c(log(1))),
            control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
summary(M5)


###############################################
## Table: Model comparison for temporal data ##
###############################################
MODELS <- list(CAR=M2, `CAR+cov`=M3, Pspline=M4, `Pspline+cov`=M5)

Table.DIC <- round(do.call(rbind,lapply(MODELS, DIC)),1)

loocv <- lapply(MODELS, function(x) inla.group.cv(x, num.level.sets=-1))
Table.LOOCV <- round(do.call(rbind,lapply(loocv, function(x) CV(x, O=data_UP$O, E=data_UP$E))),3)

lgocv <- lapply(MODELS, function(x) inla.group.cv(x, num.level.sets=3, groups=NULL))
Table.LGOCV <- round(do.call(rbind,lapply(lgocv, function(x) CV(x, O=data_UP$O, E=data_UP$E))),3)

Table <- cbind(Table.DIC,Table.LOOCV,Table.LGOCV)
print(Table)


#######################
## Plot some results ##
#######################
data_UP$risks.M1 <- M1$summary.fitted.values$`0.5quant`
data_UP$risks.M2 <- M2$summary.fitted.values$`0.5quant`
data_UP$risks.M3 <- M3$summary.fitted.values$`0.5quant`
data_UP$risks.M4 <- M4$summary.fitted.values$`0.5quant`[1:n]
data_UP$risks.M5 <- M5$summary.fitted.values$`0.5quant`[1:n]

paleta <- brewer.pal(9,"YlOrRd")
values <- c(-Inf,0.25,0.5,0.75,1,1.25,1.5,1.75,2,Inf)

tmap_mode("plot")
Mapa.risks <- tm_shape(data_UP) + 
  tm_polygons(col=c("risks.M2","risks.M3","risks.M4","risks.M5"), palette=paleta, title="Relative risk", legend.show=T,
              legend.reverse=T, style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="", main.title.position=0.2, legend.text.size=1, panel.label.size=1.5,
            panel.labels=c("CAR model","CAR + covariates","Pspline model","Pspline + covariates"),
            panel.label.bg.color="lightskyblue", legend.outside=T, legend.outside.position="right", legend.frame=F) + 
  tm_facets(nrow=1, ncol=4)
tmap_save(Mapa.risks, file="SpatialData_MapRisks.pdf")

