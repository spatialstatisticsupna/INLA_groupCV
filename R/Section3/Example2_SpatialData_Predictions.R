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


###################################################################
## Compute predictions removing areas and its k-order neighbours ##
###################################################################
n <- nrow(W)
Rs <- Diagonal(n,colSums(W))-W

n.models <- length(unique(data_UP$dist))

k <- as.list(seq(3))
names(k) <- paste("k",seq(length(k)),sep="")


## Map of spatial neighbours for ID.area=1 ##
Map <- vector("list", length(k))

for(j in k){
  if(j==1) loc <- as.numeric(which(W[1,]==1))
  if(j>1){
    Wk <- W
    for(l in 2:j){
      Wk <- Wk%*%W
      diag(Wk) <- 0
      Wk[Wk>0] <- 1
    }
    loc <- as.numeric(which(Wk[1,]==1))
  }
  
  data_UP$ID.group <- rep(NA,n)
  data_UP$ID.group[1] <- 1
  data_UP$ID.group[loc] <- 2
  Map[[j]] <- tm_shape(data_UP) + tm_polygons(col="ID.group", legend.show=F, colorNA="lightgray")
}
Map <- tmap_arrange(Map[[1]],Map[[2]],Map[[3]], ncol=3)
tmap_save(Map, file="SpatialData_SimulationStudy.pdf", width=10, height=4)

#################################################################################################
compute <- FALSE
if(compute){
  
  ## M1: Intercept + fixed effects 
  #################################
  MRPE.M1 <- RRMSPE.M1 <- lapply(k, function(x) rep(NA,n.models))
  
  for(i in seq(n.models)){
    for(j in k){
      cat("GLM",sprintf("(%d of %d):",i,n.models),sprintf("k=%d",j),"\n")
      
      formula1 <- O ~ 1 + X1 + X2 + X3
      
      data.pred <- data.frame(O=data_UP$O, E=data_UP$E, ID.area=seq(1,n),
                              X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3)
      
      if(j==1) loc <- as.numeric(which(W[i,]==1))
      if(j>1){
        Wk <- W
        for(l in 2:j){
          Wk <- Wk%*%W
          diag(Wk) <- 0
          Wk[Wk>0] <- 1
        }
        loc <- as.numeric(which(Wk[i,]==1))
      }
      
      data_UP$ID <- seq(1,n)
      data_UP$ID <- rep(1,n)
      data_UP$ID[loc] <- 2
      tm_shape(data_UP) + tm_polygons(col="ID")
      
      data.pred$O[loc] <- NA
      
      model <- inla(formula1, family="poisson", data=data.pred, E=E,
                    control.predictor=list(compute=TRUE, link=1),
                    control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
      
      y <- data_UP$O[loc]
      y.pred <- model$summary.fitted.values$`0.5quant`[loc]*data_UP$E[loc]
      MRPE.M1[[j]][i] <- mean(abs((y-y.pred)/y))
      RRMSPE.M1[[j]][i] <- sqrt(mean(((y-y.pred)/y)^2))
    }
  }
  
  
  ## M2: CAR model
  #################
  MRPE.M2 <- RRMSPE.M2 <- lapply(k, function(x) rep(NA,n.models))
  
  for(i in seq(n.models)){
    for(j in k){
      cat("CAR",sprintf("(%d of %d):",i,n.models),sprintf("k=%d",j),"\n")
      
      formula2 = O ~ 1 + 
        f(ID.area, model='bym2', graph=Rs, constr=TRUE,
          hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif)))
      
      data.pred <- data.frame(O=data_UP$O, E=data_UP$E, ID.area=seq(1,n),
                              X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3)
      
      if(j==1) loc <- as.numeric(which(W[i,]==1))
      if(j>1){
        Wk <- W
        for(l in 2:j){
          Wk <- Wk%*%W
          diag(Wk) <- 0
          Wk[Wk>0] <- 1
        }
        loc <- as.numeric(which(Wk[i,]==1))
      }
      
      data.pred$O[loc] <- NA
      
      model <- inla(formula2, family="poisson", data=data.pred, E=E,
                    control.predictor=list(compute=TRUE, link=1),
                    control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
      
      y <- data_UP$O[loc]
      y.pred <- model$summary.fitted.values$`0.5quant`[loc]*data_UP$E[loc]
      MRPE.M2[[j]][i] <- mean(abs((y-y.pred)/y))
      RRMSPE.M2[[j]][i] <- sqrt(mean(((y-y.pred)/y)^2))
    }
  }
  
  
  ## M3: CAR + covariates
  ########################
  MRPE.M3 <- RRMSPE.M3 <- lapply(k, function(x) rep(NA,n.models))
  
  for(i in seq(n.models)){
    for(j in k){
      cat("CAR + covariates",sprintf("(%d of %d):",i,n.models),sprintf("k=%d",j),"\n")
      
      formula3 = O ~ 1 + X1 + X2 + X3 + 
        f(ID.area, model='bym2', graph=Rs, constr=TRUE,
          hyper=list(prec=list(prior=sdunif), phi=list(prior=lunif)))
      
      data.pred <- data.frame(O=data_UP$O, E=data_UP$E, ID.area=seq(1,n),
                              X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3)
      
      if(j==1) loc <- as.numeric(which(W[i,]==1))
      if(j>1){
        Wk <- W
        for(l in 2:j){
          Wk <- Wk%*%W
          diag(Wk) <- 0
          Wk[Wk>0] <- 1
        }
        loc <- as.numeric(which(Wk[i,]==1))
      }
      
      data.pred$O[loc] <- NA
      
      model <- inla(formula3, family="poisson", data=data.pred, E=E,
                    control.predictor=list(compute=TRUE, link=1),
                    control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
      
      y <- data_UP$O[loc]
      y.pred <- model$summary.fitted.values$`0.5quant`[loc]*data_UP$E[loc]
      MRPE.M3[[j]][i] <- mean(abs((y-y.pred)/y))
      RRMSPE.M3[[j]][i] <- sqrt(mean(((y-y.pred)/y)^2))
    }
  }
  
  
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
  
  Ps <- list(as(kronecker(diag(k2),P1),"dgCMatrix"),
             as(kronecker(P2,diag(k1)),"dgCMatrix"))
  
  MRPE.M4 <- RRMSPE.M4 <- lapply(k, function(x) rep(NA,n.models))
  
  for(i in seq(n.models)){
    for(j in k){
      cat("Pspline",sprintf("(%d of %d):",i,n.models),sprintf("k=%d",j),"\n")
      
      formula4 <- O ~ -1 + Intercept +
        f(idx, model="generic3", Cmatrix=Ps, constr=TRUE, diagonal=1e-6,
          hyper=list(prec1=list(prior=sdunif),prec2=list(prior=sdunif)))
      
      data.pred <- list(O=data_UP$O, E=data_UP$E,
                        Intercept=c(1,rep(NA,ks)),
                        idx=c(NA,1:ks))
      
      if(j==1) loc <- as.numeric(which(W[i,]==1))
      if(j>1){
        Wk <- W
        for(l in 2:j){
          Wk <- Wk%*%W
          diag(Wk) <- 0
          Wk[Wk>0] <- 1
        }
        loc <- as.numeric(which(Wk[i,]==1))
      }
      
      data.pred$O[loc] <- NA
      
      model <- inla(formula4, family="poisson", data=data.pred, E=E,
                    control.predictor=list(compute=TRUE, A=cbind(rep(1,n), Bs), link=1, cdf=c(log(1))),
                    control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
      
      y <- data_UP$O[loc]
      y.pred <- model$summary.fitted.values$`0.5quant`[loc]*data_UP$E[loc]
      MRPE.M4[[j]][i] <- mean(abs((y-y.pred)/y))
      RRMSPE.M4[[j]][i] <- sqrt(mean(((y-y.pred)/y)^2))
    }
  }
  
  
  ## M5: P-spline + covariates
  #############################
  MRPE.M5 <- RRMSPE.M5 <- lapply(k, function(x) rep(NA,n.models))
  
  for(i in seq(n.models)){
    for(j in k){
      cat("Pspline + covariates",sprintf("(%d of %d):",i,n.models),sprintf("k=%d",j),"\n")
      
      formula5 <- O ~ -1 + Intercept + beta1 + beta2 + beta3 +
        f(idx, model="generic3", Cmatrix=Ps, constr=TRUE, diagonal=1e-6,
          hyper=list(prec1=list(prior=sdunif), prec2=list(prior=sdunif)))
      
      data.INLA <- list(O=data_UP$O, E=data_UP$E, X1=data_UP$X1, X2=data_UP$X2, X3=data_UP$X3,
                        Intercept=c(1,rep(NA,3+ks)),
                        beta1=c(rep(NA,1),1,rep(NA,2+ks)),
                        beta2=c(rep(NA,2),1,rep(NA,1+ks)),
                        beta3=c(rep(NA,3),1,rep(NA,ks)),
                        idx=c(rep(NA,4),1:ks))
      
      if(j==1) loc <- as.numeric(which(W[i,]==1))
      if(j>1){
        Wk <- W
        for(l in 2:j){
          Wk <- Wk%*%W
          diag(Wk) <- 0
          Wk[Wk>0] <- 1
        }
        loc <- as.numeric(which(Wk[i,]==1))
      }
      
      data.pred$O[loc] <- NA
      
      model <- inla(formula5, family="poisson", data=data.INLA, E=E,
                    control.predictor=list(compute=TRUE, A=cbind(rep(1,n), X1, X2, X3, Bs), link=1, cdf=c(log(1))),
                    control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE))
      
      y <- data_UP$O[loc]
      y.pred <- model$summary.fitted.values$`0.5quant`[loc]*data_UP$E[loc]
      MRPE.M5[[j]][i] <- mean(abs((y-y.pred)/y))
      RRMSPE.M5[[j]][i] <- sqrt(mean(((y-y.pred)/y)^2))
    }
  }
  
  save(list=c("MRPE.M1","MRPE.M2","MRPE.M3","MRPE.M4","MRPE.M5",
              "RRMSPE.M1","RRMSPE.M2","RRMSPE.M3","RRMSPE.M4","RRMSPE.M5"),
       file="Example2_SpatialData_Preditions.Rdata")
}


#################################
## Compute predictive measures ##
#################################
load("Example2_SpatialData_Preditions.Rdata")

Table1 <- rbind(do.call(cbind,lapply(MRPE.M2, mean)),
                do.call(cbind,lapply(MRPE.M3, mean)),
                do.call(cbind,lapply(MRPE.M4, mean)),
                do.call(cbind,lapply(MRPE.M5, mean)))
rownames(Table1) <- c("CAR","CAR + cov","Pspline","Pspline + cov")
round(Table1,3)

Table2 <- rbind(do.call(cbind,lapply(RRMSPE.M2, mean)),
                do.call(cbind,lapply(RRMSPE.M3, mean)),
                do.call(cbind,lapply(RRMSPE.M4, mean)),
                do.call(cbind,lapply(RRMSPE.M5, mean)))
rownames(Table2) <- c("CAR","CAR + cov","Pspline","Pspline + cov")
round(Table2,3)
