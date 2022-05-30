### MCMC Analysis

inits=list(list(a1=2,b1=1,a2=1.9,int=0.01,sigma=0.1),list(a1=2,b1=1,a2=-2,int=0.01,sigma=0.1),list(a1=2.1,b1=1.1,a2=-1.8,int=0.1,sigma=0.1))

n.adapt=100
n.update=500
n.iter=100

#setwd("")
dnorm_model_jags <- jags.model("CRU_model_dnorm_shapeMAT_160208.R",data=data,inits,n.chains=length(inits),n.adapt=n.adapt)
update(dnorm_model_jags,n.iter=n.update)
dnorm_model_coda=coda.samples(dnorm_model_jags,variable.names=c("a1","a2","int","sigma"),n.iter=n.iter,thin=1)
dnorm_model_jags_pred=jags.samples(dnorm_model_jags,variable.names=c("y_pred_MAT","y_pred_MAP"),n.iter=n.iter,thin=1)

plot(dnorm_model_coda)
xyplot(dnorm_model_coda)
densityplot(dnorm_model_coda)
