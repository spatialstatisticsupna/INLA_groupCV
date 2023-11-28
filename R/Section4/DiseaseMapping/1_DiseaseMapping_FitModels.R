rm(list=ls())
library(dplyr)
library(fastDummies)
library(INLA)
library(spatialreg)
library(spdep)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))
wd <- getwd()
  
## Set the "compact" mode of INLA ##
inla.getOption()$inla.mode


###################################################################
## Pancreatic cancer data of England male population (2001-2020) ##
###################################################################
load("England_CancerData.Rdata")
str(data.MORT)
str(data.INCI)

S <- length(unique(data.MORT$Region))
T <- length(unique(data.MORT$Year))
T.from <- min(data.MORT$Year)
T.to <- max(data.MORT$Year)


## Join mortality and incidence data ##
CancerData <- rbind(data.MORT[order(data.MORT$Year,data.MORT$Region),],
                    data.INCI[order(data.INCI$Year,data.INCI$Region),])
CancerData$Site <- rep(c("Mortality","Incidence"),each=S*T)
rownames(CancerData) <- NULL
str(CancerData)


## Map of Clinical Commisioning Groups of England ##
carto.nb <- mat2listw(W, style="B")
plot(carto$geometry, axes=T)
plot(carto.nb, st_centroid(st_geometry(carto), of_largest_polygon=TRUE),
     pch=19, cex=0.5, col="red", add=TRUE)


##################################################
## Model 1 -> Space-time CAR model (univariate) ##
##################################################

## Define hyperprior distribution for the precision parameters ##
sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"


## Spatial and temporal structure matrices ##
Rs <- as(Diagonal(S,colSums(W))-W, "Matrix")

Dt <- diff(diag(T), differences=1)
Rt <- as(t(Dt)%*%Dt, "Matrix")


## Define INLA data ##
J <- 2
data.INLA <- data.frame(O=CancerData$Count, E=CancerData$E, ID.disease=rep(1:J,each=S*T))

intercepts <- dummy_cols(data.INLA$ID.disease)[,-1]
intercepts[intercepts==0] <- NA
colnames(intercepts) <- paste0("I",1:J)
data.INLA <- cbind(data.INLA, intercepts)

data.INLA$ID.area1 <- rep(1:S,T*J)*data.INLA$I1
data.INLA$ID.area2 <- rep(1:S,T*J)*data.INLA$I2

data.INLA$ID.year1 <- rep(rep(1:T,each=S),J)*data.INLA$I1
data.INLA$ID.year2 <- rep(rep(1:T,each=S),J)*data.INLA$I2

data.INLA$ID.area.year1 <- rep(1:(S*T),J)*data.INLA$I1
data.INLA$ID.area.year2 <- rep(1:(S*T),J)*data.INLA$I2


## TypeI interaction model ##
f.TypeI <- O ~ -1 + I1 + I2 + 
  f(ID.area1, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area2, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.year1, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year2, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year2, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))

M1.TypeI <- inla(f.TypeI, family="poisson", data=data.INLA, E=E,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M1.TypeI$loocv <- inla.group.cv(M1.TypeI, num.level.sets=-1)
M1.TypeI$lgocv.m3 <- inla.group.cv(M1.TypeI, num.level.sets=3)
M1.TypeI$lgocv.m5 <- inla.group.cv(M1.TypeI, num.level.sets=5)
M1.TypeI$lgocv.m10 <- inla.group.cv(M1.TypeI, num.level.sets=10)


## TypeII interaction model ##
R <- kronecker(Rt,Diagonal(S))
r.def <- S
A.constr <- kronecker(matrix(1,1,T),Diagonal(S))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeII <- O ~ -1 + I1 + I2 + 
  f(ID.area1, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area2, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.year1, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year2, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M1.TypeII <- inla(f.TypeII, family="poisson", data=data.INLA, E=E,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M1.TypeII$loocv <- inla.group.cv(M1.TypeII, num.level.sets=-1)
M1.TypeII$lgocv.m3 <- inla.group.cv(M1.TypeII, num.level.sets=3)
M1.TypeII$lgocv.m5 <- inla.group.cv(M1.TypeII, num.level.sets=5)
M1.TypeII$lgocv.m10 <- inla.group.cv(M1.TypeII, num.level.sets=10)


## TypeIII interaction model ##
R <- kronecker(Diagonal(T),Rs)
r.def <- T
A.constr <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ -1 + I1 + I2 + 
  f(ID.area1, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area2, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.year1, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year2, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M1.TypeIII <- inla(f.TypeIII, family="poisson", data=data.INLA, E=E,
                   control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M1.TypeIII$loocv <- inla.group.cv(M1.TypeIII, num.level.sets=-1)
M1.TypeIII$lgocv.m3 <- inla.group.cv(M1.TypeIII, num.level.sets=3)
M1.TypeIII$lgocv.m5 <- inla.group.cv(M1.TypeIII, num.level.sets=5)
M1.TypeIII$lgocv.m10 <- inla.group.cv(M1.TypeIII, num.level.sets=10)


## TypeIV interaction model ##
R <- kronecker(Rt,Rs)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T),Diagonal(S))
A2 <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ -1 + I1 + I2 + 
  f(ID.area1, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area2, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) +
  f(ID.year1, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year2, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M1.TypeIV <- inla(f.TypeIV, family="poisson", data=data.INLA, E=E,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M1.TypeIV$loocv <- inla.group.cv(M1.TypeIV, num.level.sets=-1)
M1.TypeIV$lgocv.m3 <- inla.group.cv(M1.TypeIV, num.level.sets=3)
M1.TypeIV$lgocv.m5 <- inla.group.cv(M1.TypeIV, num.level.sets=5)
M1.TypeIV$lgocv.m10 <- inla.group.cv(M1.TypeIV, num.level.sets=10)


## Save the models ##
setwd(wd)

M1 <- list(TypeI=M1.TypeI, TypeII=M1.TypeII, TypeIII=M1.TypeIII, TypeIV=M1.TypeIV)
save(M1, file="DiseaseMapping_UnivariateModels.Rdata")


####################################
## Model 2 -> Space-time M-models ##
####################################

## Download source code from GitHub repository ##
download.file(url="https://github.com/spatialstatisticsupna/BookChapter_STMmodels/archive/master.zip",
              destfile="BookChapter_STMmodels.zip")

unzip("BookChapter_STMmodels.zip")

setwd("BookChapter_STMmodels-main/R")
source("MCAR_INLA_ST_Model1.R")


## Fit the models ##
CancerData$ID.disease <- data.INLA$ID.disease

type <- list(TypeI="TypeI", TypeII="TypeII", TypeIII="TypeIII", TypeIV="TypeIV")

M2 <- lapply(type, function(x){
  model <- MCAR_INLA_ST(carto=carto, data=CancerData, O="Count", E="E",
                        ID.area="Region", ID.year="Year", ID.disease="ID.disease",
                        spatial="intrinsic", temporal="rw1", interaction=x)
  
  model$loocv <- inla.group.cv(model, num.level.sets=-1)
  model$lgocv.m3 <- inla.group.cv(model, num.level.sets=3)
  model$lgocv.m5 <- inla.group.cv(model, num.level.sets=5)
  model$lgocv.m10 <- inla.group.cv(model, num.level.sets=10)
  
  model
})


## Save the models ##
setwd(wd)
save(M2, file="DiseaseMapping_Mmodels.Rdata")


###################################################
## Model 3 -> Space-time shared component models ##
###################################################

## Define hyperprior distribution for the precision parameters ##
sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"


## Spatial and temporal structure matrices ##
Rs <- as(Diagonal(S,colSums(W))-W, "Matrix")

Dt <- diff(diag(T), differences=1)
Rt <- as(t(Dt)%*%Dt, "Matrix")


## Define INLA data ##
J <- 2
data.INLA <- data.frame(O=CancerData$Count, E=CancerData$E, ID.disease=rep(1:J,each=S*T))

intercepts <- dummy_cols(data.INLA$ID.disease)[,-1]
intercepts[intercepts==0] <- NA
colnames(intercepts) <- paste0("I",1:J)
data.INLA <- cbind(data.INLA, intercepts)

data.INLA$ID.area <- rep(1:S,T*J) + (data.INLA$ID.disease-1)*S
data.INLA$ID.year <- rep(rep(1:T,each=S),J) + (data.INLA$ID.disease-1)*T
data.INLA$ID.area.year1 <- c(1:(S*T),rep(NA,S*T))
data.INLA$ID.area.year2 <- c(rep(NA,S*T),1:(S*T))


## TypeI interaction model ##
f.TypeI <- O ~ - 1 + I1 + I2 + 
  f(ID.area, model='besag2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model='besag2', graph=Rt, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year2, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))

M3.TypeI <- inla(f.TypeI, family="poisson", data=data.INLA, E=E,
                 control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                 control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M3.TypeI$loocv <- inla.group.cv(M3.TypeI, num.level.sets=-1)
M3.TypeI$lgocv.m3 <- inla.group.cv(M3.TypeI, num.level.sets=3)
M3.TypeI$lgocv.m5 <- inla.group.cv(M3.TypeI, num.level.sets=5)
M3.TypeI$lgocv.m10 <- inla.group.cv(M3.TypeI, num.level.sets=10)


## TypeII interaction model ##
R <- kronecker(Rt,Diagonal(S))
r.def <- S
A.constr <- kronecker(matrix(1,1,T),Diagonal(S))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeII <- O ~ - 1 + I1 + I2 + 
  f(ID.area, model='besag2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model='besag2', graph=Rt, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M3.TypeII <- inla(f.TypeII, family="poisson", data=data.INLA, E=E,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M3.TypeII$loocv <- inla.group.cv(M3.TypeII, num.level.sets=-1)
M3.TypeII$lgocv.m3 <- inla.group.cv(M3.TypeII, num.level.sets=3)
M3.TypeII$lgocv.m5 <- inla.group.cv(M3.TypeII, num.level.sets=5)
M3.TypeII$lgocv.m10 <- inla.group.cv(M3.TypeII, num.level.sets=10)


## TypeIII interaction model ##
R <- kronecker(Diagonal(T),Rs)
r.def <- T
A.constr <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(A.constr[-1,],"matrix")

f.TypeIII <- O ~ - 1 + I1 + I2 + 
  f(ID.area, model='besag2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model='besag2', graph=Rt, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M3.TypeIII <- inla(f.TypeIII, family="poisson", data=data.INLA, E=E,
                   control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                   control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M3.TypeIII$loocv <- inla.group.cv(M3.TypeIII, num.level.sets=-1)
M3.TypeIII$lgocv.m3 <- inla.group.cv(M3.TypeIII, num.level.sets=3)
M3.TypeIII$lgocv.m5 <- inla.group.cv(M3.TypeIII, num.level.sets=5)
M3.TypeIII$lgocv.m10 <- inla.group.cv(M3.TypeIII, num.level.sets=10)


## TypeIV interaction model ##
R <- kronecker(Rt,Rs)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T),Diagonal(S))
A2 <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

f.TypeIV <- O ~ - 1 + I1 + I2 + 
  f(ID.area, model='besag2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.year, model='besag2', graph=Rt, constr=TRUE, hyper=list(prec=list(prior=sdunif))) + 
  f(ID.area.year1, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr)))) + 
  f(ID.area.year2, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)),
    extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))

M3.TypeIV <- inla(f.TypeIV, family="poisson", data=data.INLA, E=E,
                  control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                  control.compute=list(dic=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE))

M3.TypeIV$loocv <- inla.group.cv(M3.TypeIV, num.level.sets=-1)
M3.TypeIV$lgocv.m3 <- inla.group.cv(M3.TypeIV, num.level.sets=3)
M3.TypeIV$lgocv.m5 <- inla.group.cv(M3.TypeIV, num.level.sets=5)
M3.TypeIV$lgocv.m10 <- inla.group.cv(M3.TypeIV, num.level.sets=10)


## Save the models ##
setwd(wd)

M3 <- list(TypeI=M3.TypeI, TypeII=M3.TypeII, TypeIII=M3.TypeIII, TypeIV=M3.TypeIV)
save(M3, file="DiseaseMapping_SharedModels.Rdata")
