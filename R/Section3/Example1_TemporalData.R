rm(list=ls())
library(INLA)  # R-INLA stable version 23.09.09

setwd(dirname(rstudioapi::getSourceEditorContext()$path))


#########################
## Auxiliary functions ##
#########################
DIC <- function(x){
  data.frame(mean.deviance=x$dic$mean.deviance,
             p.eff=x$dic$p.eff,
             dic=x$dic$dic,
             waic=x$waic$waic)
}

CV <- function(x){
  LS <- -mean(log(x$cv))
  
  # Computation of E[Y_i | y_{I_i}] using trapezoid numerical integration #
  expectation = numeric(n)
  for (i in 1:n){
    mu = x$mean[i]
    sd = x$sd[i]
    xx = seq(mu-6*sd,mu+6*sd,length.out = 100)
    yy = xx*dnorm(xx, mean=mu, sd=sd)
    expectation[i] = pracma::trapz(xx,yy)
  }
  MSPE <- mean((expectation-temperature)^2)
  
  data.frame(LS=LS, MSPE=MSPE)
}

sdunif="expression:
  logdens=-log_precision/2;
  return(logdens)"


###################
## Simulate data ##
###################
set.seed(123)
n <- 2000
id <- 1:n
rho <- 0.9
var_marginal <- 1
sd_marginal <- sqrt(var_marginal)
var_additive <- (1-rho^2)*var_marginal
sd_additive <- sqrt(var_additive)
radiation <- arima.sim(model = list(ar = rho),n = n,sd = sd_additive)
sd.error <- 0.05
temperature <- 30 + radiation + rnorm(n, sd=sd.error)
radiation_sq <- sign(radiation)*abs(radiation)^2


##############################
## Fit the models with INLA ##
##############################
formula1 <- temperature ~ 1 + radiation
formula2 <- temperature ~ 1 + f(time, model='ar1')

M1 <- inla(formula=formula1, data=list(temperature=temperature, radiation=radiation_sq),
           control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE))

M2 <- inla(formula=formula2, data=list(temperature=temperature, time=id),
           control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE),
           control.family=list(hyper=list(prec=list(initial=log(1/sd.error), fixed=TRUE))))


#######################
## Plot some results ##
#######################
graphics.off()
pdf("TemporalData_tempVSradiation.pdf", width=10)
par(cex.axis=1.5, cex.lab=1.5)
plot(radiation_sq, temperature, ylim=c(25,35), xlab="Radiation", ylab="Temperature", cex.axis=1.5)
points(sort(M1$.args$data$radiation), sort(M1$summary.fitted.values$`0.5quant`), type="l", col="blue", lwd=2)
points(sort(M1$.args$data$radiation), sort(M2$summary.fitted.values$`0.5quant`), type="l", col="red", lwd=2)
legend("topleft", col=c("blue","red"), legend=c("Linear model","AR1 model"), lty=1, bty="n", cex=1.5, lwd=2)
dev.off()

graphics.off()
pdf("TemporalData_ModelFitting.pdf", width=16)
par(mfrow=c(1,2), cex.main=2, cex.axis=1.5, cex.lab=1.5)
plot(temperature, ylim=c(25,35), type="p", pch=19,  main="Linear model")
lines(M1$summary.fitted.values$`0.5quant`, col="blue")
plot(temperature, ylim=c(25,35), type="p", pch=19, main="AR1 model")
lines(M2$summary.fitted.values$`0.5quant`, col="red")
dev.off()


#################################################
## Table 1: Model comparison for temporal data ##
#################################################
Table.DIC <- round(do.call(rbind,lapply(list(Linear=M1, AR1=M2), DIC)),1)

## Leave-one-out cross validation ##
loocv1 <- inla.group.cv(result=M1, num.level.sets=-1)
loocv2 <- inla.group.cv(result=M2, num.level.sets=-1)
Table.LOOCV <- round(do.call(rbind,lapply(list(Linear=loocv1, AR1=loocv2), CV)),3)

## Leave-group-out cross validation ##
## The groups derived from the AR1 model are used as a reference ## 
m <- 3

# Based on prior correlations #
lgocv2.prior <- inla.group.cv(result=M2, num.level.sets=m, strategy="prior")
groups <- lapply(lgocv2.prior$groups, function(x) x$idx)
lgocv1.prior <- inla.group.cv(result=M1, groups=groups)
Table.LGOCV.prior <- round(do.call(rbind,lapply(list(Linear.prior=lgocv1.prior, AR1.prior=lgocv2.prior), CV)),3)

# Based on posterior correlations #
lgocv2.post <- inla.group.cv(result=M2, num.level.sets=m, strategy="posterior")
groups <- lapply(lgocv2.post$groups, function(x) x$idx)
lgocv1.post <- inla.group.cv(result=M1, groups=groups)
Table.LGOCV.post <- round(do.call(rbind,lapply(list(Linear.post=lgocv1.post, AR1.post=lgocv2.post), CV)),3)


## Print the table ##
Table1 <- cbind(Table.DIC,Table.LOOCV,Table.LGOCV.prior,Table.LGOCV.post)
print(Table1)
