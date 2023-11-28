rm(list=ls())
library(INLA)  # R-INLA stable version 23.09.09

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

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


#########################################
## Compute t-steps forward predictions ##
#########################################
n.pred <- 100
n.models <- 500

compute <- FALSE

if(compute){
  
  MAPE <- list(Linear=rep(NA,n.models), AR1=rep(NA,n.models))
  RMSPE <- list(Linear=rep(NA,n.models), AR1=rep(NA,n.models))
  
  
  ## Linear models 
  #################
  models <- vector("list",n.models)
  
  for(i in seq(0,n.models-1)){
    cat("Linear models",i+1,"of",n.models,"\n")
    
    formula1 = temperature ~ 1 + radiation
    
    data.pred <- list(temperature=c(temperature[seq(1,n-n.pred-i)], rep(NA,n.pred)),
                      radiation=radiation_sq[seq(1,n-i)])
    nn <- length(data.pred$temperature)
    
    models[[i+1]] <- inla(formula=formula1, data=data.pred,
                          control.compute=list(cpo=TRUE, dic=TRUE, waic=TRUE, config=TRUE))
    
    y <- temperature[seq(nn-n.pred+1,nn)]
    y.pred <- models[[i+1]]$summary.fitted.values$`0.5quant`[seq(nn-n.pred+1,nn)]
    MAPE$Linear[i+1] <- mean(abs(y-y.pred))
    RMSPE$Linear[i+1] <- sqrt(mean((y-y.pred)^2))
  }
  Linear.model <- models[[1]]
  rm(models)
  gc()
  
  
  ## AR1 models 
  ###############
  models <- vector("list",n.models)
  
  for(i in seq(0,n.models-1)){
    cat("AR1 models",i+1,"of",n.models,"\n")
    
    formula2 = temperature ~ 1 + f(time, model='ar1')
    
    data.pred <- list(temperature=c(temperature[seq(1,n-n.pred-i)], rep(NA,n.pred)),
                      radiation=radiation_sq[seq(1,n-i)])
    nn <- length(data.pred$temperature)
    data.pred$time <- 1:nn
    
    models[[i+1]] <- inla(formula=formula2, data=data.pred,
                          control.compute=list(cpo=TRUE,  dic=TRUE, waic=TRUE, config=TRUE),
                          control.family=list(hyper=list(prec=list(initial=log(1/sd.error),fixed=TRUE))))
    
    y <- temperature[seq(nn-n.pred+1,nn)]
    y.pred <- models[[i+1]]$summary.fitted.values$`0.5quant`[seq(nn-n.pred+1,nn)]
    MAPE$AR1[i+1] <- mean(abs(y-y.pred))
    RMSPE$AR1[i+1] <- sqrt(mean((y-y.pred)^2))
  }
  AR1.model <- models[[1]]
  rm(models)
  gc()
  
  save(list=c("Linear.model","AR1.model","MAPE","RMSPE"),
       file="Example1_TemporalData_Preditions.Rdata")
}


#################################
## Compute predictive measures ##
#################################
load("Example1_TemporalData_Preditions.Rdata")

Table <- data.frame(MAPE=do.call(rbind,lapply(MAPE, mean, na.rm=T)),
                    RMSPE=do.call(rbind,lapply(RMSPE, mean, na.rm=T)))
round(Table,3)


##########################################################
## Plot of the predictions for the last 100 time points ##
##########################################################
graphics.off()
pdf("TemporalData_Predictions.pdf", width=16)

par(mfrow=c(1,2), cex.main=2, cex.axis=1.5, cex.lab=1.5)

plot(temperature, xlim=c(1500,2000), ylim=c(26,34), type="p", pch=19, cex=0.3, main="Linear model", cex.main=2)
lines(Linear.model$summary.fitted.values$`0.5quant`, col="red")
lines(Linear.model$summary.fitted.values$mean[seq(1,n-n.pred)], col="blue")
abline(v=n-n.pred, lty=2)

plot(temperature, xlim=c(1500,2000), ylim=c(26,34), type="p", pch=19, cex=0.3, main="AR1 model", cex.main=2)
lines(AR1.model$summary.fitted.values$`0.5quant`, col="red")
lines(AR1.model$summary.fitted.values$mean[seq(1,n-n.pred)], col="blue")
abline(v=n-n.pred, lty=2)

dev.off()
