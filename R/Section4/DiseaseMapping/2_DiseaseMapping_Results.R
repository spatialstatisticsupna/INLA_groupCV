rm(list=ls())
library(INLA)
library(RColorBrewer)
library(tmap)

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
options(timeout=9999)


###################################
## Load previously fitted models ##
###################################
load("DiseaseMapping_UnivariateModels.Rdata")
load("DiseaseMapping_Mmodels.Rdata")
load("DiseaseMapping_SharedModels.Rdata")

## Available also from URL ##
# load(url("https://emi-sstcdapp.unavarra.es/INLA_groupCV/DiseaseMapping_UnivariateModels.Rdata"))
# load(url("https://emi-sstcdapp.unavarra.es/INLA_groupCV/DiseaseMapping_Mmodels.Rdata"))
# load(url("https://emi-sstcdapp.unavarra.es/INLA_groupCV/DiseaseMapping_SharedModels.Rdata"))


#########################
## Auxiliary functions ##
#########################
DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance,
             p.eff=x$dic$p.eff,
             DIC=x$dic$dic,
             WAIC=x$waic$waic)
}

sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"


########################################################################################
## TABLE 6. Pancreatic cancer data: model selection criteria (DIC) and predictive     ##
##          performance measures (LS) with automatic groups construction for m=3,5,10 ##
########################################################################################
load("England_CancerData.Rdata")

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


#######################
## Univariate models ##
#######################
Rs <- as(Diagonal(S,colSums(W))-W, "Matrix")
Dt <- diff(diag(T), differences=1)
Rt <- as(t(Dt)%*%Dt, "Matrix")

Table.M1 <- round(do.call(rbind, lapply(M1,DIC)),1)
Table.M1$LOOCV <- round(unlist(lapply(M1, function(x) -mean(log(x$loocv$cv)))),3)

M1.groups.m3 <- lapply(M1$TypeII$lgocv.m3$groups, function(x) x$idx)
M1.groups.m5 <- lapply(M1$TypeII$lgocv.m5$groups, function(x) x$idx)
M1.groups.m10 <- lapply(M1$TypeII$lgocv.m10$groups, function(x) x$idx)

lgocv.M1 <- data.frame(m3=rep(NA,4), m5=rep(NA,4), m10=rep(NA,4))
row.names(lgocv.M1) <- names(M1)


## Type I interaction ##
aux.m3 <- inla.group.cv(M1$TypeI, groups=M1.groups.m3)
aux.m5 <- inla.group.cv(M1$TypeI, groups=M1.groups.m5)
aux.m10 <- inla.group.cv(M1$TypeI, groups=M1.groups.m10)
lgocv.M1["TypeI",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))


## Type II interaction ##
R <- kronecker(Rt,Diagonal(S))
r.def <- S
A.constr <- kronecker(matrix(1,1,T),Diagonal(S))
A.constr <- as(A.constr[-1,],"matrix")

aux.m3 <- inla.group.cv(M1$TypeII, groups=M1.groups.m3)
aux.m5 <- inla.group.cv(M1$TypeII, groups=M1.groups.m5)
aux.m10 <- inla.group.cv(M1$TypeII, groups=M1.groups.m10)
lgocv.M1["TypeII",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))


## Type III interaction ##
R <- kronecker(Diagonal(T),Rs)
r.def <- T
A.constr <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(A.constr[-1,],"matrix")

aux.m3 <- inla.group.cv(M1$TypeIII, groups=M1.groups.m3)
aux.m5 <- inla.group.cv(M1$TypeIII, groups=M1.groups.m5)
aux.m10 <- inla.group.cv(M1$TypeIII, groups=M1.groups.m10)
lgocv.M1["TypeIII",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))


## Type IV interaction ##
R <- kronecker(Rt,Rs)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T),Diagonal(S))
A2 <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

aux.m3 <- inla.group.cv(M1$TypeIV, groups=M1.groups.m3)
aux.m5 <- inla.group.cv(M1$TypeIV, groups=M1.groups.m5)
aux.m10 <- inla.group.cv(M1$TypeIV, groups=M1.groups.m10)
lgocv.M1["TypeIV",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))


##############
## M-models ##
##############
Table.M2 <- round(do.call(rbind, lapply(M2,DIC)),1)
Table.M2$LOOCV <- round(unlist(lapply(M2, function(x) -mean(log(x$loocv$cv)))),3)

M2.groups.m3 <- lapply(M2$TypeII$lgocv.m3$groups, function(x) x$idx)
M2.groups.m5 <- lapply(M2$TypeII$lgocv.m5$groups, function(x) x$idx)
M2.groups.m10 <- lapply(M2$TypeII$lgocv.m10$groups, function(x) x$idx)

aux.m3 <- lapply(M2, function(x) inla.group.cv(x, groups=M2.groups.m3))
aux.m5 <- lapply(M2, function(x) inla.group.cv(x, groups=M2.groups.m5))
aux.m10 <- lapply(M2, function(x) inla.group.cv(x, groups=M2.groups.m10))

Table.M2$m3 <- unlist(lapply(aux.m3, function(x) -mean(log(x$cv))))
Table.M2$m5 <- unlist(lapply(aux.m5, function(x) -mean(log(x$cv))))
Table.M2$m10 <- unlist(lapply(aux.m10, function(x) -mean(log(x$cv))))


#############################
## Shared-component models ##
#############################
Table.M3 <- round(do.call(rbind, lapply(M3,DIC)),1)
Table.M3$LOOCV <- round(unlist(lapply(M3, function(x) -mean(log(x$loocv$cv)))),3)

M3.groups.m3 <- lapply(M3$TypeII$lgocv.m3$groups, function(x) x$idx)
M3.groups.m5 <- lapply(M3$TypeII$lgocv.m5$groups, function(x) x$idx)
M3.groups.m10 <- lapply(M3$TypeII$lgocv.m10$groups, function(x) x$idx)

lgocv.M3 <- data.frame(m3=rep(NA,4), m5=rep(NA,4), m10=rep(NA,4))
row.names(lgocv.M3) <- names(M3)


## Type I interaction ##
aux.m3 <- inla.group.cv(M3$TypeI, groups=M3.groups.m3)
aux.m5 <- inla.group.cv(M3$TypeI, groups=M3.groups.m5)
aux.m10 <- inla.group.cv(M3$TypeI, groups=M3.groups.m10)
lgocv.M3["TypeI",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))


## Type II interaction ##
R <- kronecker(Rt,Diagonal(S))
r.def <- S
A.constr <- kronecker(matrix(1,1,T),Diagonal(S))
A.constr <- as(A.constr[-1,],"matrix")

aux.m3 <- inla.group.cv(M3$TypeII, groups=M3.groups.m3)
aux.m5 <- inla.group.cv(M3$TypeII, groups=M3.groups.m5)
aux.m10 <- inla.group.cv(M3$TypeII, groups=M3.groups.m10)
lgocv.M3["TypeII",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))



## Type III interaction ##
R <- kronecker(Diagonal(T),Rs)
r.def <- T
A.constr <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(A.constr[-1,],"matrix")

aux.m3 <- inla.group.cv(M3$TypeIII, groups=M3.groups.m3)
aux.m5 <- inla.group.cv(M3$TypeIII, groups=M3.groups.m5)
aux.m10 <- inla.group.cv(M3$TypeIII, groups=M3.groups.m10)
lgocv.M3["TypeIII",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))



## Type IV interaction ##
R <- kronecker(Rt,Rs)
r.def <- S+T-1
A1 <- kronecker(matrix(1,1,T),Diagonal(S))
A2 <- kronecker(Diagonal(T),matrix(1,1,S))
A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")

aux.m3 <- inla.group.cv(M3$TypeIV, groups=M3.groups.m3)
aux.m5 <- inla.group.cv(M3$TypeIV, groups=M3.groups.m5)
aux.m10 <- inla.group.cv(M3$TypeIV, groups=M3.groups.m10)
lgocv.M3["TypeIV",] <- unlist(lapply(list(aux.m3,aux.m5,aux.m10), function(x) -mean(log(x$cv))))



#######################
## Print the results ##
#######################
Table.M1 <- cbind(Table.M1,lgocv.M1)
print(round(Table.M1[,c("DIC","m3","m5","m10")],3))

print(round(Table.M2[,c("DIC","m3","m5","m10")],3))

Table.M3 <- cbind(Table.M3,lgocv.M3)
print(round(Table.M3[,c("DIC","m3","m5","m10")],3))


######################################################################################################
## FIGURE A1. Top: posterior median estimates $\exp(\phi_i)$ and posterior exceedence probabilities ##
##            $P(\exp(\phi_i)>1 | \boldsymbol{y})$ for the shared spatial random effect (M3).       ##
##            Bottom: Posterior median estimates and 95\% credible intervals of the temporal        ##
##            shared random effect (M3).                                                            ##
######################################################################################################
Model <- M3$TypeII

####################
## Spatial effect ##
####################
spatial.pattern <- lapply(Model$marginals.random$ID.area[1:S], function(x) inla.tmarginal(function(y) exp(y), x))

carto$spatial.mean <- unlist(lapply(spatial.pattern, function(x) inla.qmarginal(0.5,x)))
carto$spatial.prob <- unlist(lapply(spatial.pattern, function(x) 1-inla.pmarginal(1,x)))

Map1 <- tm_shape(carto) +
  tm_polygons(col="spatial.mean", palette=brewer.pal(5,"YlOrRd"),
              title="Risk", legend.show=T, legend.reverse=T,
              style="fixed", breaks=c(-Inf,0.91,0.95,1,1.05,1.1,Inf), interval.closure="left") + 
  tm_layout(legend.title.size=1.5, legend.text.size=1)

Map2 <- tm_shape(carto) +
  tm_polygons(col="spatial.prob", palette=brewer.pal(6,"Blues")[-1],
              title="Prob", legend.show=T, legend.reverse=T,
              style="fixed", breaks=c(0,0.1,0.2,0.8,0.9,1), interval.closure="left",
              labels=c("[0-0.1)","[0.1-0.2)","[0.2-0.8)","[0.8-0.9)","[0.9-1]")) + 
  tm_layout(legend.title.size=1.5, legend.text.size=1)

Map <- tmap_arrange(Map1, Map2, nrow=1)

tmap_mode("plot")
tmap_save(Map, file="FigA1_SpatialPattern.pdf", width=12, height=6)


#####################
## Temporal effect ##
#####################
temporal.pattern <- lapply(Model$marginals.random$ID.year[1:T], function(x) inla.tmarginal(function(y) exp(y), x))

temporal.mean <- unlist(lapply(temporal.pattern, function(x) inla.qmarginal(0.5, x)))
q1 <- unlist(lapply(temporal.pattern, function(x) inla.qmarginal(0.025,x)))
q2 <- unlist(lapply(temporal.pattern, function(x) inla.qmarginal(0.975,x)))

graphics.off()
pdf("FigA1_TemporalPattern.pdf", height=6)

x <- 1:T
plot(range(x), c(0.9, 1.1), type="n", xlab="",ylab="", xaxt="n", main="", cex.lab=1.3, cex.axis=1.3, cex.main=1.5)
axis(1, at=round(seq(1,T,3)), labels=round(seq(T.from,T.to,3)), las=0, cex.axis=1.3)
X.Vec <- c(x, tail(x, 1), rev(x), x[1])
Y.Vec <- c(q1, tail(q2, 1), rev(q2), q1[1])
polygon(X.Vec, Y.Vec, col = "grey", border = NA)
lines(temporal.mean)
abline(h=1, lty=2)

dev.off()


##############################################################################################
## FIGURE A2. Posterior median estimates of pancreatic cancer mortality and incidence risks ##
##            for the shared-component model (M3) with Type II interaction effect.          ##
##############################################################################################
Model <- M3$TypeII

RR1 <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==1],S,T,byrow=F)
RR2 <- matrix(Model$summary.fitted.values$`0.5quant`[Model$.args$data$ID.disease==2],S,T,byrow=F)
colnames(RR1) <- colnames(RR2) <- paste("Year",seq(T.from,T.to),sep=".")

paleta <- brewer.pal(8,"RdYlGn")[8:1]
values <- c(-Inf,0.83,0.91,0.95,1,1.05,1.1,1.20,Inf)


####################
## Mortality data ##
####################
carto.aux <- cbind(carto,RR1)

Map.RR1 <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",round(seq(T.from,T.to,length.out=6)),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Pancreatic cancer mortality data", main.title.position=0.2, panel.label.size=1.5,
            panel.labels=paste("Year",round(seq(T.from,T.to,length.out=6))),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_save(Map.RR1, file="FigA2_MortalityRR.pdf")


####################
## Incidence data ##
####################
carto.aux <- cbind(carto,RR2)

Map.RR2 <- tm_shape(carto.aux) +
  tm_polygons(col=paste("Year",round(seq(T.from,T.to,length.out=6)),sep="."), id="name", palette=paleta,
              title="", legend.show=T, legend.reverse=T,
              style="fixed", breaks=values, interval.closure="left") +
  tm_layout(main.title="Pancreatic cancer incidence data", main.title.position=0.2, panel.label.size=1.5,
            panel.labels=paste("Year",round(seq(T.from,T.to,length.out=6))),
            legend.outside=T, legend.outside.position="right", legend.frame=F,
            legend.outside.size=0.2, outer.margins=c(0.02,0.01,0.02,0.01)) + 
  tm_facets(nrow=2, ncol=3)

tmap_save(Map.RR2, file="FigA2_IncidenceRR.pdf")
