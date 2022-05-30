### MCMC Analysis

inits=list(list(mu1=0.5,tau1=0.2,a2=1.9,b2=1,int=0.01,sigma=0.1),list(mu1=0.1,tau1=0.2,a2=1.9,b2=1,int=0.01,sigma=0.1),list(mu1=0.8,tau1=0.2,a2=1.9,b2=1,int=0.01,sigma=0.1))

n.adapt=100
n.update=500
n.iter=100

# setwd('')
dnorm_model_jags <- jags.model("CRU_model_dnorm_shapeMATMAP_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt) # logistic MAT, quadratic logistic MAP
dnorm_model_jags <- jags.model("CRU_model_dnorm_shapeMATMAP_logistic_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt) # logistic MAT, MAP
dnorm_model_jags <- jags.model("CRU_model_dnorm_shapeMAT_lognormal_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt) # logistic MAT, MAP
dnorm_model_jags <- jags.model("CRU_model_dnorm_shapeMAT_normal_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt) # logistic MAT, MAP

inits=list(list(mu1=0.5,tau1=0.2,a2=-0.1,b2=1.8,c2=2,d2=0.4,int=0.01,sigma=0.1),list(mu1=0.1,tau1=0.2,a2=-0.2,b2=1.2,c2=2,d2=0.4,int=0.01,sigma=0.1),list(mu1=0.8,tau1=0.2,a2=-0.2,b2=1.2,c2=2,d2=0.5,int=0.01,sigma=0.1))
dnorm_model_jags <- jags.model("dnorm_shapeMAT_normal_MAP_sh_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt) # logistic MAT, MAP
update(dnorm_model_jags,n.iter=n.update)
dnorm_model_coda=coda.samples(dnorm_model_jags,variable.names=c("mu1","tau1","a2","b2","c2","d2","int","sigma"),n.iter=n.iter,thin=1)
dnorm_model_jags_pred=jags.samples(dnorm_model_jags,variable.names=c("y_pred_MAT","y_pred_MAP"),n.iter=n.iter,thin=1)

plot(dnorm_model_coda)
xyplot(dnorm_model_coda)
densityplot(dnorm_model_coda)

####################
#################### # Get predictions
####################

fix_jags_pred <- dnorm_model_jags_pred
fix_jags <- dnorm_model_jags
fix_coda <- dnorm_model_coda

### MAT
y_pred_MAT <- summary(fix_jags_pred$y_pred_MAT, quantile,c(0.025,0.5,0.975))$stat
lowerfit_MAT<-y_pred_MAT[1,]
medianfit_MAT<-y_pred_MAT[2,]
upperfit_MAT<-y_pred_MAT[3,]

colnames(y_pred_MAT)<-MAT_x
dat_pred_MAT<-t(y_pred_MAT)
dat_pred_MAT<-as.data.frame(dat_pred_MAT)
dat_pred_MAT$MAT_x<-rownames(dat_pred_MAT)

### MAP
y_pred_MAP <- summary(fix_jags_pred$y_pred_MAP, quantile,c(0.025,0.5,0.975))$stat
lowerfit_MAP<-y_pred_MAP[1,]
medianfit_MAP<-y_pred_MAP[2,]
upperfit_MAP<-y_pred_MAP[3,]

colnames(y_pred_MAP)<-MAP_x
dat_pred_MAP<-t(y_pred_MAP)
dat_pred_MAP<-as.data.frame(dat_pred_MAP)
dat_pred_MAP$MAP_x<-rownames(dat_pred_MAP)

####################
#################### # Save parameter values
####################

all_model_stat<-cbind(summary(fix_coda)$statistics,summary(fix_coda)$quantiles)
all_model_stat<-as.data.frame(all_model_stat)
#setwd("")
write.csv(all_model_stat,file="CRU_all_dnorm_model_stat_logistic_160208.csv")

#########################################################################################
##################### MS MS MS MS MS MS MS MS MS MS MS MS MS MS PLOT#####################
#########################################################################################

#setwd("")
dat_ba_mat_fix<-read.csv("CRU_dat_ba_mat_fix_160119.csv",header=T)
dat_ba_map_fix<-read.csv("CRU_dat_ba_map_fix_160119.csv",header=T)
dat_ba_mat_fix<-dat_ba_mat_fix[,c(-1)]
dat_ba_map_fix<-dat_ba_map_fix[,c(-1)]

BA_clm$fMAT<-floor(BA_clm$MAT)
BA_clm$fMAP<-floor(BA_clm$MAP/10)*10
ba_mean_mat_fix<-tapply(BA_clm$perc_fix_all,BA_clm$fMAT,mean)
ba_mat_fix<-as.data.frame(ba_mean_mat_fix)
ba_mat_fix$fMAT<-as.double(rownames(ba_mat_fix))
ba_mean_map_fix<-tapply(BA_clm$perc_fix_all,BA_clm$fMAP,mean)
ba_map_fix<-as.data.frame(ba_mean_map_fix)
ba_map_fix$fMAP<-as.double(rownames(ba_map_fix))

#setwd("")
pdf("CRU_Fig1_dnorm_MATfit_logistic_160208.pdf",width=4.5,height=4.5)

plot(ba_mat_fix$ba_mean_mat*100~ba_mat_fix$fMAT,cex=1,pch=19,col="grey",xlab=expression(paste("Mean Annual Temperature (",degree,"C)")),ylab="N Fixer Basal Area (%)",ylim=c(0,80),main="Temperature",xlim=c(-5,30))
points(medianfit_MAT*100~MAT_x,typ="l",lwd=2)
points(lowerfit_MAT*100~MAT_x,col="red",typ="l",lwd=2,lty=2)
points(upperfit_MAT*100~MAT_x,col="red",typ="l",lty=2,lwd=2)
legend("topright",col=c("grey","black","red"),legend=c("Mean BA% of Data","Best Fit of Data","95% CI of Best Fit"),lty=c(NA,1,2),pch=c(19,NA,NA),bty='n',lwd=2)

#setwd("")
pdf("CRU_Fig2_dnorm_logistic_MAPfit_160208.pdf",width=4.5,height=4.5)

plot(ba_map_fix$ba_mean_map*100~ba_map_fix$fMAP,cex=1,pch=19,col="grey",main="Precipitation",xlab="Mean Annual Precipitation (mm)",ylab="N Fixer Basal Area (%)",ylim=c(0,30),xlim=c(0,420))
points(medianfit_MAP*100~MAP_x,typ="l",lwd=2)
points(lowerfit_MAP*100~MAP_x,col="red",typ="l",lwd=2,lty=2)
points(upperfit_MAP*100~MAP_x,col="red",typ="l",lty=2,lwd=2)
legend("topright",col=c("grey","black","red"),legend=c("Mean BA% of Data","Best Fit of Data","95% CI of Best Fit"),lty=c(NA,1,2),pch=c(19,NA,NA),bty='n',lwd=2)

dev.off()