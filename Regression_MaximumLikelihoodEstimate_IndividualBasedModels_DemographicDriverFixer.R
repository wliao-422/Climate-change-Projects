# R code created to analyze the demographic drivers of nitrogen-fixer distribution across succession in U.S.
# Author: WL

################################################################
######################       Libraries       ###################
################################################################

library(bbmle) # mle2
library(MASS) # mvrnorm

################################################################
######################    Working Directory   ##################
################################################################

# setwd("")

################################################################
######################       Original Data    ##################
################################################################

# Species: species.csv
# FIA_fixerspcd_merged_noCeltis.csv
	# GRIN report: GRIN_nodulation_reports_Fixer0.6.csv
	# Huss Dannel report: Actinorhizal_spcd_FIA_HussD.csv
# DATA SHARE POSSIBLE UPON REQUEST

################################################################
######################       Species Info     ##################
################################################################

# setwd("")
sp<-read.csv("species.csv",header=T)

robinia<-sp_spcd[sp_spcd$GENUS=="Robinia",]
cercocarpus<-sp_spcd[sp_spcd$GENUS=="Cercocarpus",]
alnus<-sp_spcd[sp_spcd$GENUS=="Alnus",]
albizia<-sp_spcd[sp_spcd$GENUS=="Albizia",]
elaeagnus<-sp_spcd[sp_spcd$GENUS=="Elaeagnus",]
prosopis<-sp_spcd[sp_spcd$GENUS=="Prosopis",]

################################################################
#####################  Assign Fixing Status   ##################
################################################################


# setwd("")
sp<-read.csv("species.csv",header=T)
grin_sp<-read.csv("GRIN_nodulation_reports_greater0.6.csv",header=T)

# Merge
sp_spcd<-sp[,1:4]
grin_sp_FIA<-merge(sp_spcd,grin_sp,by.x="GENUS",by.y="Row.Labels",all.x=FALSE,all.y=TRUE)

#setwd("")
#write.csv(grin_sp_FIA,file="GRIN_nodulation_reports_greater0.6_spcdMerged.csv")


####################################################################################################
################################              Growth Analysis         ##############################
####################################################################################################

######################################################
##################  Rate Calculation  ################
######################################################

# setwd("")
t<-read.csv("tree_small.csv",header=TRUE)

ft<-t[,c(2,4,5,7)] #pcn, tcn, spcd, dbh
ts<-read.csv("tree_ts.csv",header=TRUE) # Load in tree time series data
ts1<-merge(ts,ft,by.x="tcn1",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(ts1)[c(10,11,12)]<-c("pcn1","spcd1","dbh1")
ts2<-merge(ts1,ft,by.x="tcn2",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(ts2)[c(13,14,15)]<-c("pcn2","spcd2","dbh2")
ts3<-merge(ts2,ft,by.x="tcn3",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(ts3)[c(16:18)]<-c("pcn3","spcd3","dbh3")
rm(ts2)
ts4<-merge(ts3,ft,by.x="tcn4",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(ts4)[c(19:21)]<-c("pcn4","spcd4","dbh4")
rm(ts3)
rm(ft)

tbind<-ts4[!is.na(ts4$dbh1)|!is.na(ts4$dbh2)|!is.na(ts4$dbh3)|!is.na(ts4$dbh4),]

#setwd("")
#write.csv(tbind,file="FIAGRSTAT_all_tree_ts_150707.csv")

tb12<-tbind[,c(3,4,6,7,c(10:15))]
tb23<-tbind[,c(2,3,7,8,c(13:18))]
tb34<-tbind[,c(1,2,8,9,c(16:21))]

tb12[tb12$mt1<0,]$mt1<-NA
tb12[tb12$mt2<0,]$mt2<-NA

tb23[tb23$mt2<0,]$mt2<-NA
tb23[tb23$mt3<0,]$mt3<-NA

tb34[tb34$mt3<0,]$mt3<-NA
tb34[tb34$mt4<0,]$mt4<-NA

# get rid of NA's for dbh (diameter) and mt(measurement time)
tb.12<-tb12[!is.na(tb12$mt1)&!is.na(tb12$mt2)&!is.na(tb12$dbh1)&!is.na(tb12$dbh2),]
tb.23<-tb23[!is.na(tb23$mt2)&!is.na(tb23$mt3)&!is.na(tb23$dbh2)&!is.na(tb23$dbh3),]
tb.34<-tb34[!is.na(tb34$mt3)&!is.na(tb34$mt4)&!is.na(tb34$dbh3)&!is.na(tb34$dbh4),]

colnames(tb.23)<-c("tcn2","tcn1","mt1","mt2","pcn1","spcd1","dbh1","pcn2","spcd2","dbh2")
colnames(tb.34)<-c("tcn2","tcn1","mt1","mt2","pcn1","spcd1","dbh1","pcn2","spcd2","dbh2")
tb1234<-rbind(tb.12,tb.23,tb.34)

# setwd("")
# write.csv(tb1234,file="FIAGRSTAT_all_tree_ts_rbind_150707.csv")

######### Step2: Assign fixing status

fixer<-read.csv("FIA_fixerspcd_merged_noCeltis.csv",header=T) 
tree_spcd<-merge(tb1234,fixer,by.x="spcd1",by.y="SPCD",all.x=TRUE,all.y=FALSE)
tree_spcd[is.na(tree_spcd$FIX),]$FIX<-0

#setwd("")
#write.csv(tree_spcd,file="FIAGRSTAT_tree_plot_fixerstatus_150707.csv")

############ Step3: Get the plots with fixer present

#setwd("")
p<-read.csv("plot.csv",header=T) # Load in plot data
fp<-p[p$natforsamp==1&p$propsamp==1&p$ntlog==0,] # Screen plot data by following criteria: naturally regenerating, all trees sampled (100% sampling), not affected by logging
ffp<-fp[,c(2,8)] # pcn (plot number), std(stand age)
gdata<-merge(ffp,tree_spcd,by.x="pcn",by.y="pcn1",all.x=FALSE,all.y=TRUE) 

gdata<-gdata[gdata$stdage>=0&!is.na(gdata$stdage),]
gdata_fp_pcn<-unique(gdata[gdata$FIX==1,]$pcn) #unique plot number with fixer present
gdata_fp<-data.frame(pcn=gdata_fp_pcn,fixer_present=1)
gdata_fp1<-merge(gdata_fp,gdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
gdata_fp2<-gdata_fp1[!is.na(gdata_fp1$fixer_present),]

gdata<-gdata_fp2[gdata_fp2$dbh2>0 & gdata_fp2$dbh1>0,]

#setwd("")
#write.csv(gdata,file="FIAGRSTAT_tree_plotswithfixer_150707.csv")

######### Step4: Compute relative growth rates of individuals for each year

t<-gdata
t$lg<-log(t$dbh2)-log(t$dbh1)

	# Initialzize positivized relative growth rates
t$lgp<-log(t$dbh2)-log(t$dbh1)

	#Compute "positivized" relative growth rates following approach in Condit et al. 2006 Science (Eqn.S5)
MDL<-0.05 # Minimum diameter measurement is 1
hMDL<-MDL/2
t[t$lg<=0&!is.na(t$lg),]$lgp<-log(t[t$lg<=0&!is.na(t$lg),]$dbh1+hMDL)-log(t[t$lg<=0&!is.na(t$lg),]$dbh1)

t$LGR<-t$lg/(t$mt2-t$mt1)
t$LGRp<-t$lgp/(t$mt2-t$mt1)
t$t<-t$mt2-t$mt1

tree<-data.frame(pcn=t$pcn,tcn=t$tcn1,stdage=t$stdage,year=t$mt1,t=t$t,spcd=t$spcd1,FIX=t$FIX,dbh=t$dbh1,LGR=t$LGR,LGRp=t$LGRp)

#setwd("")
#write.csv(tree,file="FIAGRSTAT_all_tree_log_plot_LGR_FixerPresent_150707.csv")

	#Get rid of the data above 99th percentile and below 1st percentile of growth rate (uncorrected)
gdata<-tree[tree$LGR>quantile(tree$LGR,0.01) & tree$LGR<quantile(tree$LGR,0.99),]
gdata_all<-gdata

#setwd("")
#write.csv(gdata,file="FIAGRSTAT_final_data_for_fit_150707.csv") 

	########## Here: print out the no. of individuals for each species

gdata$count<-1
spcd_gr<-tapply(gdata[gdata$FIX==1,]$count,gdata[gdata$FIX==1,]$spcd,sum)
spcd_GR<-data.frame(spcd=rownames(spcd_gr),count=spcd_gr)
spcd_GR<-merge(spcd_GR,fixer,by.x="spcd",by.y="SPCD",all.x=TRUE,all.y=FALSE)
spcd_GR_rank<-spcd_GR[order(-spcd_GR$count),]

#setwd("")
#write.csv(spcd_GR,file="FIAGRSTAT_FixerInfo_150707.csv")

######################################################
##################  Mean Calculation  ################
######################################################

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

u_stdage_all<-u_stdage
u_meanLGR_fix_all<-u_meanLGR_fix
u_meanLGR_non_all<-u_meanLGR_non

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig1_GRmean_150707.pdf",width=4,height=4)

plot(u_stdage_all,u_meanLGR_fix_all*100,col="red",xlab="Stand Age (Year)",ylab="Growth Rate (% diameter growth/year)", main="Growth",pch=19)
points(u_stdage_all,u_meanLGR_non_all*100,col="blue",pch=17)
legend("topright",col=c('red','blue'),legend=c('Fixer','Non-Fixer'),pch=c(19,17),bty='n')

dev.off()

######################################################
##################    STAT Analysis   ################
######################################################

library(bbmle)

MDL<-0.05
lMDL<-function(D,T){
	log((log(D+MDL)-log(D))/T)
}

################### Saturating ###################

###### For N Fixers ######

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood	
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)

	# Plot Fit with Initial Values
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)

	# Plot Fits with Initial Values
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Sigmoidal ###################
###### For N Fixers ######

	# Functions
gMLF_2_f<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_f<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_f(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_f(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-7,-5,32,-0.5,-7)

	# Plot Fits with Initial Values


	# MLE
fit_logGR_2_f<-mle2(normNNLF_2_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5],k=pv_f[6]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_2_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-7,-5,32,-0.5,-7)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Sigmoidal Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],pv_f[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_2_f)
curve(exp(gMLF_2_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],pv_f_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_2_n<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_n<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_n(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_n(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.1)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_n<-mle2(normNNLF_2_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5],k=pv_n[6]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_2_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.1)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Sigmoidal Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_n[2],pv_n[3],pv_n[4],pv_n[5],pv_n[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_2_n)
curve(exp(gMLF_2_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],pv_n_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Ricker ###################
###### For N Fixers ######

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_2_f),AIC(fit_logGR_4_f)))

# [1] 558.5292
# [1] 183.4808
# [1] 0

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_2_n),AIC(fit_logGR_4_n)))

# [1] 0
# [1] 8.822845
# [1] 11.56831

#setwd("")
pv_f_gr_all<-coef(fit_logGR_4_f)
pv_n_gr_all<-coef(fit_logGR_1_n)
write.csv(pv_f_gr_all,file="pv_f_gr_all_updated.csv")
write.csv(pv_n_gr_all,file="pv_n_gr_all_updated.csv")

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

# How many times you want to run it (draws from the multivariate normal)
ndraws<-1000
# Create an x vector
xs<-seq(0, 250, 1)
xs_all<-seq(0,250,1)
# Multi-variate normal draw from parameter means and variance-covariance matrices
# from an object that came from an mle2 fit
vmat_fix<-mvrnorm(ndraws, mu=coef(fit_logGR_4_f)[c(2:5)],Sigma=vcov(fit_logGR_4_f)[c(2:5),c(2:5)])
vmat_non<-mvrnorm(ndraws, mu=coef(fit_logGR_1_n)[c(2:5)],Sigma=vcov(fit_logGR_1_n)[c(2:5),c(2:5)])
# New array of y values, x long for each draw
ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))
# For each draw, fit the likelihood function using mle2
# Then write into the new y array for that draw
for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_4_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_1_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}
civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))
# For each x, write in the upper and lower CIs for quantiles corresponding to lower bound and upper bound
# Set to 2.5 and 97.5 for 95% CI
lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_all<-civec_fix
civec_non_all<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig1_GRfit_wCI_Updated_150928.pdf",width=4.5,height=4.5)

agevec_noOG_all<-sort(unique(gdata_all$stdage))
plot(u_stdage_all,u_meanLGR_fix_all*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (% diameter growth/years)", main="Growth",cex=1,pch=1,ylim=c(0,4))
points(u_stdage_all,u_meanLGR_non_all*100,col="blue",cex=1,pch=2,ylim=c(0,4))

# pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_gr_all[2],pv_f_gr_all[3],pv_f_gr_all[4],pv_f_gr_all[5],exp(mean(log(gdata_all$dbh))),x))*100,from=min(agevec_noOG_all),to=max(agevec_noOG_all),lwd=3,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
points(civec_fix_all[1,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_all[2,]~xs,lwd=2,lty=2,typ="l",col="red")

# pv_n_fit<-coef(fit_logGR_5_n)
curve(exp(gMLF_1_n(pv_n_gr_all[2],pv_n_gr_all[3],pv_n_gr_all[4],pv_n_gr_all[5],exp(mean(log(gdata_all$dbh))),x))*100,from=min(agevec_noOG_all),to=max(agevec_noOG_all),lwd=3,lty=1,col="blue",add=TRUE)
points(civec_non_all[1,]~xs_all,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_all[2,]~xs_all,lwd=2,lty=2,typ="l",col="blue")

legend("topright", bty="n",legend=c("N fixers","Non-fixers "),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2),cex=1.2)

dev.off()

######################################################
############    Robinia Growth Analysis  #############
######################################################

gdata.original<-gdata_all
gdata<-gdata.original

# Get the plots with robinia present

gdata_rob_pcn<-unique(gdata[gdata$spcd==901,]$pcn) #Robinia pseudoacacia present
gdata_rob<-data.frame(pcn=gdata_rob_pcn,rob_present=1)
gdata_rob1<-merge(gdata_rob,gdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
gdata_rob2<-gdata_rob1[!is.na(gdata_rob1$rob_present),]

gdata_robinia<-gdata_rob2[(gdata_rob2$FIX==1 & gdata_rob2$spcd==901) | gdata_rob2$FIX==0,]
gdata<-gdata_robinia

#setwd("")
write.csv(gdata_robinia,file="FIA_GrowthData_RobiniaPseudoacacia_150709.csv")

######################################################
##################  Mean Calculation  ################
######################################################

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

u_stdage_robinia<-u_stdage
u_meanLGR_fix_robinia<-u_meanLGR_fix
u_meanLGR_non_robinia<-u_meanLGR_non

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig2_GRmean_Robinia_150709.pdf",width=4,height=4)

plot(u_stdage_robinia,u_meanLGR_fix_robinia*100,col="red",xlab="Stand Age (Year)",ylab="Growth Rate (% diameter growth/year)", main=expression(italic(Robinia~pseudoacacia)),pch=19)
points(u_stdage_robinia,u_meanLGR_non_robinia*100,col="blue",pch=17)
legend("topright",col=c('red','blue'),legend=c(expression(italic(Robinia~pseudoacacia)),'Non-Fixer'),pch=c(19,17),bty='n')

dev.off()

######################################################
##################    STAT Analysis   ################
######################################################

MDL<-0.05
lMDL<-function(D,T){
	log((log(D+MDL)-log(D))/T)
}

################### Saturating ###################

###### For N Fixers ######

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)

	# Plot Fit with Initial Values
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)

	# Plot Fits with Initial Values
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Sigmoidal ###################
###### For N Fixers ######

	# Functions
gMLF_2_f<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_f<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_f(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_f(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-7,-5,32,-0.5,-7)

	# Plot Fits with Initial Values


	# MLE
fit_logGR_2_f<-mle2(normNNLF_2_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5],k=pv_f[6]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_2_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-7,-5,32,-0.5,-7)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Sigmoidal Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],pv_f[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_2_f)
curve(exp(gMLF_2_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],pv_f_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_2_n<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_n<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_n(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_n(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.1)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_n<-mle2(normNNLF_2_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5],k=pv_n[6]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_2_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.1)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Sigmoidal Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_n[2],pv_n[3],pv_n[4],pv_n[5],pv_n[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_2_n)
curve(exp(gMLF_2_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],pv_n_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Ricker ###################
###### For N Fixers ######

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_2_f),AIC(fit_logGR_4_f)))

# [1] 3.201478
# [1] 3.079871
# [1] 0

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_2_n),AIC(fit_logGR_4_n)))

# [1] 55.89526
# [1] 26.47894
# [1] 0

#setwd("")
pv_f_fit_robinia<-coef(fit_logGR_4_f)
pv_n_fit_robinia<-coef(fit_logGR_4_n)
write.csv(pv_f_fit_robinia,file="pv_f_gr_robinia_updated.csv")
write.csv(pv_n_fit_robinia,file="pv_n_gr_robinia_updated.csv")

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

# How many times you want to run it (draws from the multivariate normal)
ndraws<-1000
# Create an x vector
agevec_noOG<-sort(unique(gdata$stdage))
xs<-seq(min(agevec_noOG), max(agevec_noOG), 1)
xs_robinia<-seq(min(agevec_noOG), max(agevec_noOG), 1)
# Multi-variate normal draw from parameter means and variance-covariance matrices
# from an object that came from an mle2 fit
vmat_fix<-mvrnorm(ndraws, mu=coef(fit_logGR_4_f)[c(2:5)],Sigma=vcov(fit_logGR_4_f)[c(2:5),c(2:5)])
vmat_non<-mvrnorm(ndraws, mu=coef(fit_logGR_4_n)[c(2:5)],Sigma=vcov(fit_logGR_4_n)[c(2:5),c(2:5)])
# New array of y values, x long for each draw
ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))
# For each draw, fit the likelihood function using mle2
# Then write into the new y array for that draw
for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_4_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_4_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}
civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))
# For each x, write in the upper and lower CIs for quantiles corresponding to lower bound and upper bound
# Set to 2.5 and 97.5 for 95% CI
lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

pv_f_fit_robinia<-coef(fit_logGR_4_f)
pv_n_fit_robinia<-coef(fit_logGR_4_n)

civec_fix_robinia<-civec_fix
civec_non_robinia<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig2_GRfit_wCI_Robinia_Updated_150928.pdf",width=4.5,height=4.5)

agevec_noOG_robinia<-sort(unique(gdata_robinia$stdage))
plot(u_stdage_robinia,u_meanLGR_fix_robinia*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (% diameter growth/year)", main=expression(italic(Robinia~pseudoacacia)),cex=1,pch=1,ylim=c(0,4))
points(u_stdage_robinia,u_meanLGR_non_robinia*100,col="blue",cex=1,pch=2,ylim=c(0,4))

# pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit_robinia[2],pv_f_fit_robinia[3],pv_f_fit_robinia[4],pv_f_fit_robinia[5],exp(mean(log(gdata_robinia$dbh))),x))*100,from=min(agevec_noOG_robinia),to=max(agevec_noOG_robinia),lwd=3,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
points(civec_fix_robinia[1,]~xs_robinia,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_robinia[2,]~xs_robinia,lwd=2,lty=2,typ="l",col="red")

# pv_n_fit<-coef(fit_logGR_5_n)
curve(exp(gMLF_4_n(pv_n_fit_robinia[2],pv_n_fit_robinia[3],pv_n_fit_robinia[4],pv_n_fit_robinia[5],exp(mean(log(gdata_robinia$dbh))),x))*100,from=min(agevec_noOG_robinia),to=max(agevec_noOG_robinia),lwd=3,lty=1,col="blue",add=TRUE)
points(civec_non_robinia[1,]~xs_robinia,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_robinia[2,]~xs_robinia,lwd=2,lty=2,typ="l",col="blue")

legend("topright", bty="n",legend=c(expression(italic(Robinia~pseudoacacia)),"Non-fixers "),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2),cex=1.2)

dev.off()

######################################################
############    Alnus   Growth Analysis  #############
######################################################

gdata.original<-gdata_all
gdata<-gdata.original

# Get the plots with Alnus rubrum present

gdata_alnr_pcn<-unique(gdata[gdata$spcd==351,]$pcn) #Alnus rubrum present
gdata_alnr<-data.frame(pcn=gdata_alnr_pcn,alnr_present=1)
gdata_alnr1<-merge(gdata_alnr,gdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
gdata_alnr2<-gdata_alnr1[!is.na(gdata_alnr1$alnr_present),]

gdata_alnus<-gdata_alnr2[(gdata_alnr2$FIX==1 & gdata_alnr2$spcd==351) | gdata_alnr2$FIX==0,]
gdata<-gdata_alnus

#setwd("")
write.csv(gdata_alnus,file="FIA_GrowthData_AlnusRubra_150709.csv")

######################################################
##################  Mean Calculation  ################
######################################################

gdata<-gdata_alnus

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="Mean Growth Rate",pch=1)
points(u_stdage,u_meanLGR_non*100,col="blue",pch=2)
legend("topright",col=c('red','blue'),legend=c('Fixer','Non-Fixer'),pch=c(1,2),bty='n')

u_stdage_alnus<-u_stdage
u_meanLGR_fix_alnus<-u_meanLGR_fix
u_meanLGR_non_alnus<-u_meanLGR_non

######################################################
##################    STAT Analysis   ################
######################################################

library(bbmle)

MDL<-0.05
lMDL<-function(D,T){
	log((log(D+MDL)-log(D))/T)
}

################### Saturating ###################

###### For N Fixers ######

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.5,-4,30,-0.2)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.5,-4,30,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-6,-3.9,25,-0.3)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-6,-3.9,25,-0.3)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Sigmoidal ###################

###### For N Fixers ######

	# Functions
gMLF_2_f<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_f<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_f(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_f(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.2,-3.6,15,-0.2,-0.1)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_f<-mle2(normNNLF_2_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5],k=pv_f[6]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_2_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.2,-3.6,15,-0.2,-0.1)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Sigmoidal Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],pv_f[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_2_f)
curve(exp(gMLF_2_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],pv_f_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_2_n<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_n<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_n(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_n(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.4)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_n<-mle2(normNNLF_2_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5],k=pv_n[6]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_2_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-4.9,-3.6,20,-0.2,-0.4)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Sigmoidal Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_n[2],pv_n[3],pv_n[4],pv_n[5],pv_n[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_2_n)
curve(exp(gMLF_2_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],pv_n_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Ricker ###################
###### For N Fixers ######

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-7,-4,-0.3,0.01)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_2_f),AIC(fit_logGR_4_f)))

# [1] 0
# [1] 2.153264
# [1] 0.1535219

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_2_n),AIC(fit_logGR_4_n)))

# [1] 1.681435
# [1] 6.535772
# [1] 0

#setwd("")
pv_f_fit_alnus<-coef(fit_logGR_1_f)
pv_n_fit_alnus<-coef(fit_logGR_4_n)
write.csv(pv_f_fit_alnus,file="pv_f_fit_alnus_updated.csv")
write.csv(pv_n_fit_alnus,file="pv_n_fit_alnus_updated.csv")

agevec_noOG_alnus<-sort(unique(gdata_alnus$stdage))

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

ndraws<-10000
xs<-seq(min(agevec_noOG_alnus),max(agevec_noOG_alnus), 1)
xs_alnus<-seq(min(agevec_noOG_alnus),max(agevec_noOG_alnus), 1)
vmat_fix<-mvrnorm(ndraws, mu=coef(fit_logGR_1_f)[c(2:5)],Sigma=vcov(fit_logGR_1_f)[c(2:5),c(2:5)])
vmat_non<-mvrnorm(ndraws, mu=coef(fit_logGR_4_n)[c(2:5)],Sigma=vcov(fit_logGR_4_n)[c(2:5),c(2:5)])

ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_1_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_4_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975

for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_alnus<-civec_fix
civec_non_alnus<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig3_GRfit_wCI_Alnus_Updated_150928.pdf",width=4.5,height=4.5)

agevec_noOG_alnus<-sort(unique(gdata_alnus$stdage))
plot(u_stdage_alnus,u_meanLGR_fix_alnus*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (%diameter growth/year)", main=expression(italic(Alnus~rubra)),cex=1,pch=1,ylim=c(0,4))
points(u_stdage_alnus,u_meanLGR_non_alnus*100,col="blue",cex=1,pch=2,ylim=c(0,4))

# pv_f_fit_alnus<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit_alnus[2],pv_f_fit_alnus[3],pv_f_fit_alnus[4],pv_f_fit_alnus[5],exp(mean(log(gdata_alnus$dbh))),x))*100,from=min(agevec_noOG_alnus),to=max(agevec_noOG_alnus),lwd=3,lty=1,col="red",add=TRUE)

# pv_n_fit_alnus<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit_alnus[2],pv_n_fit_alnus[3],pv_n_fit_alnus[4],pv_n_fit_alnus[5],exp(mean(log(gdata_alnus$dbh))),x))*100,from=min(agevec_noOG_alnus),to=max(agevec_noOG_alnus),lwd=3,lty=1,col="blue",add=TRUE)

legend("topright", bty="n",legend=c(expression(italic(Alnus~rubra)),"Non-fixers"),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2))

points(civec_fix_alnus[1,]~xs_alnus,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_alnus[2,]~xs_alnus,lwd=2,lty=2,typ="l",col="red")
points(civec_non_alnus[1,]~xs_alnus,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_alnus[2,]~xs_alnus,lwd=2,lty=2,typ="l",col="blue")

dev.off()

######################################################
############ Cercocarpus Growth Analysis #############
######################################################

gdata.original<-gdata_all
gdata<-gdata.original

# Get the plots with cercocarpus present

gdata_crc_pcn<-unique(gdata[gdata$spcd==475,]$pcn) #Cercocarpus ledifolius present
gdata_crc<-data.frame(pcn=gdata_crc_pcn,crc_present=1)
gdata_crc1<-merge(gdata_crc,gdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
gdata_crc2<-gdata_crc1[!is.na(gdata_crc1$crc_present),]

gdata_cercocarpus<-gdata_crc2[(gdata_crc2$FIX==1 & gdata_crc2$spcd==475) | gdata_crc2$FIX==0,]
gdata<-gdata_cercocarpus

#setwd("")
write.csv(gdata_cercocarpus,file="FIA_GrowthData_CercocarpusLedifolius_150709.csv")

######################################################
##################  Mean Calculation  ################
######################################################

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="Mean Growth Rate",pch=1)
points(u_stdage,u_meanLGR_non*100,col="blue",pch=2)
legend("topright",col=c('red','blue'),legend=c('Fixer','Non-Fixer'),pch=c(1,2),bty='n')

u_stdage_cercocarpus<-u_stdage
u_meanLGR_fix_cercocarpus<-u_meanLGR_fix
u_meanLGR_non_cercocarpus<-u_meanLGR_non

######################################################
##################    STAT Analysis   ################
######################################################

library(bbmle)

MDL<-0.05
lMDL<-function(D,T){
	log((log(D+MDL)-log(D))/T)
}

################### Saturating ###################

###### For N Fixers ######

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6,-4,30,-0.2)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6,-4,30,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5.5,-3.9,45,-0.59)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Sigmoidal ###################
###### For N Fixers ######

	# Functions
gMLF_2_f<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_f<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_f(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_f(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.2,-3.6,15,-0.2,-0.1)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_f<-mle2(normNNLF_2_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5],k=pv_f[6]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_2_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.2,-3.6,15,-0.2,-0.1)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Sigmoidal Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],pv_f[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_2_f)
curve(exp(gMLF_2_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],pv_f_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_2_n<-function(a,b,c,d,k,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D)))
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

	# Negative loglikelihood
normNNLF_2_n<-function(sd,a,b,c,d,k,G,D,A,T){
	up<-gMLF_2_n(a,b,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_2_n(a,b,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1.3,-5.6,-5,-7,-1,-0.3)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_2_n<-mle2(normNNLF_2_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5],k=pv_n[6]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_2_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1.3,-5.6,-5,17,-0.98,-0.3)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Sigmoidal Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_2_f(pv_n[2],pv_n[3],pv_n[4],pv_n[5],pv_n[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_2_n)
curve(exp(gMLF_2_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],pv_n_fit[6],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################### Ricker ###################

###### For N Fixers ######

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

###### For Non Fixers ######

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1.5,-6,-2.3,-01,0.001)

	# Plot Fits with Initial Values
	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1.5,-6,-2.3,-01,0.001)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_2_f),AIC(fit_logGR_4_f)))

# [1] 3.421301
# [1] 0
# [1] 16.96502

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_2_n),AIC(fit_logGR_4_n)))

# [1] 0
# [1] 42.20651
# [1] 40.20651

#setwd("")
pv_f_fit_cercocarpus<-coef(fit_logGR_2_f)
pv_n_fit_cercocarpus<-coef(fit_logGR_1_n)
write.csv(pv_f_fit_cercocarpus,file="pv_f_gr_cercocarpus.csv")
write.csv(pv_n_fit_cercocarpus,file="pv_n_gr_cercocarpus.csv")

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

ndraws<-10000
agevec_noOG_cercocarpus<-sort(unique(gdata_cercocarpus$stdage))
xs_cercocarpus<-seq(min(agevec_noOG_cercocarpus),max(agevec_noOG_cercocarpus), 1)
xs<-seq(min(agevec_noOG_cercocarpus),max(agevec_noOG_cercocarpus), 1)
vmat_fix<-mvrnorm(ndraws, mu=coef(fit_logGR_2_f)[c(2:6)],Sigma=vcov(fit_logGR_2_f)[c(2:6),c(2:6)])
vmat_non<-mvrnorm(ndraws, mu=coef(fit_logGR_1_n)[c(2:5)],Sigma=vcov(fit_logGR_1_n)[c(2:5),c(2:5)])

ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_2_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],vmat_fix[i,5],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_1_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975

for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_cercocarpus<-civec_fix
civec_non_cercocarpus<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig4_GRfit_wCI_Cercocarpus_Updated_150928.pdf",width=4.5,height=4.5)

agevec_noOG_cercocarpus<-sort(unique(gdata_cercocarpus$stdage))
plot(u_stdage_cercocarpus,u_meanLGR_fix_cercocarpus*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (%diameter growth/year)", main=expression(italic(Cercocarpus~ledifolius)),cex=1,pch=1,ylim=c(0,20))
points(u_stdage_cercocarpus,u_meanLGR_non_cercocarpus*100,col="blue",cex=1,pch=2,ylim=c(0,4))

# pv_f_fit_cercocarpus<-coef(cerc_fix_best)
curve(exp(gMLF_2_f(pv_f_fit_cercocarpus[2],pv_f_fit_cercocarpus[3],pv_f_fit_cercocarpus[4],pv_f_fit_cercocarpus[5],pv_f_fit_cercocarpus[6],exp(mean(log(gdata_cercocarpus$dbh))),x))*100,from=min(agevec_noOG_cercocarpus),to=max(agevec_noOG_cercocarpus),lwd=3,lty=1,col="red",add=TRUE)

# pv_n_fit_cercocarpus<-coef(cerc_non_best)
curve(exp(gMLF_1_n(pv_n_fit_cercocarpus[2],pv_n_fit_cercocarpus[3],pv_n_fit_cercocarpus[4],pv_n_fit_cercocarpus[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG_cercocarpus),to=max(agevec_noOG_cercocarpus),lwd=3,lty=1,col="blue",add=TRUE)

legend("topright", bty="n",legend=c(expression(italic(Cercocarpus~ledifolius)),"Non-fixers"),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2))

points(civec_fix_cercocarpus[1,]~xs_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_cercocarpus[2,]~xs_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_non_cercocarpus[1,]~xs_cercocarpus,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_cercocarpus[2,]~xs_cercocarpus,lwd=2,lty=2,typ="l",col="blue")

dev.off()




#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################




####################################################################################################
################################           Mortality Analysis         ##############################
####################################################################################################

######################################################
##################  Rate Calculation  ################
######################################################

#setwd("")

t<-read.csv("tree_small.csv",header=TRUE) # Load in tree data
ft<-t[,c(2,4,5,6,7)] #pcn, tcn, spcd, statuscd, dbh
t.ts<-read.csv("tree_ts.csv",header=TRUE) # Load in tree time series data
rm(t)

ts1<-t.ts[,c(2,3,6,7)]
ts2<-t.ts[,c(3,4,7,8)]
ts3<-t.ts[,c(4,5,8,9)]

colnames(ts1)<-c("tcn0","tcnT","mt0","mtT")
colnames(ts2)<-c("tcn0","tcnT","mt0","mtT")
colnames(ts3)<-c("tcn0","tcnT","mt0","mtT")
ts12<-rbind(ts1,ts2)
ts123<-rbind(ts12,ts3)
rts<-ts123[!is.na(ts123$tcn0)&!is.na(ts123$tcnT),]

mtts1<-merge(rts,ft,by.x="tcn0",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(mtts1)[5:8]<-c("pcn0","spcd0","statuscd0","dbh0")

fft<-ft[,c(2,4)]
mtts2<-merge(mtts1,fft,by.x="tcnT",by.y="tcn",all.x=TRUE,all.y=FALSE)
colnames(mtts2)[9]<-"statuscdT"

mtts<-mtts2
#setwd("")
write.csv(mtts,file="FIAMRSTAT_merged_ts_tree_status_150709.csv")

mtts<-mtts[mtts$statuscd0==1 & (mtts$statuscdT==1|mtts$statuscdT==2),] # alive=1, dead=2
mtts$N<-NA
mtts$S<-NA
mtts[mtts$statuscd0==1,]$N<-1
mtts[mtts$statuscd0==1&mtts$statuscdT==1,]$S<-1 # survived
mtts[mtts$statuscd0==1&mtts$statuscdT==2,]$S<-0 # dead

#setwd("")
write.csv(mtts,file="FIAMRSTAT_merged_ts_tree_NS_150709.csv")
rm(t)

#setwd("")
p<-read.csv("plot.csv",header=T)
fp<-p[p$natforsamp==1&p$propsamp==1&p$ntlog==0,]
ffp<-fp[,c(2,8)]
mdat<-merge(ffp,mtts,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

#setwd("")
fixer<-read.csv("FIA_fixerspcd_merged_noCeltis.csv",header=T) 
mdat_fixer<-merge(mdat,fixer,by.x="spcd0",by.y="SPCD",all.x=TRUE,all.y=FALSE)
mdat_fixer[is.na(mdat_fixer$FIX),]$FIX<-0

mdat<-mdat_fixer

mdat<-mdat[mdat$stdage>=0&!is.na(mdat$stdage),]
mdat_fp_pcn<-unique(mdat[mdat$FIX==1,]$pcn) #fixer present
mdat_fp<-data.frame(pcn=mdat_fp_pcn,fixer_present=1)
mdat_fp1<-merge(mdat_fp,mdat,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
mdat_fp2<-mdat_fp1[!is.na(mdat_fp1$fixer_present),]

mdat<-mdat_fp2

mdat$mt<-NA
mdat[mdat$mtT>0,]$mt<-mdat[mdat$mtT>0,]$mtT-mdat[mdat$mtT>0,]$mt0

mdata<-mdat[!is.na(mdat$mt),]

#setwd("")
write.csv(mdata,file="FIAMRSTAT_morality_data_stat_final_150709.csv")

########## Here: print out the no. of individuals for each species

mdata$count<-1
spcd_mr<-tapply(mdata[mdata$FIX==1,]$count,mdata[mdata$FIX==1,]$spcd0,sum)
spcd_MR<-data.frame(spcd=rownames(spcd_mr),count=spcd_mr)
spcd_MR<-merge(spcd_MR,fixer,by.x="spcd",by.y="SPCD",all.x=FALSE,all.y=FALSE)
spcd_MR<-spcd_MR[order(-spcd_MR$count),]

#setwd("")
write.csv(spcd_MR,file="FIAMRSTAT_FixerInfo_150709.csv")

mdata.original<-mdata

######################################################
##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1
######################################################

#Now: Fixers

mdata<-mdata.original
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_all_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata.original
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_all_v1<-favg
navg_mr<-favg

######################################################
##################  Mean Calculation  ################ Version 2: plot with S=0, omit
######################################################

#Now: Fixers

mdata<-mdata.original
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

pAgeNS_fix_omit<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_all_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata.original
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

pAgeNS_non_omit<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_all_v2<-favg
navg_mr<-favg

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig3_MRmean_150714.pdf",width=9,height=4.5)

par(mfrow=c(1,2))

plot(fxavg_mr_all_v1$meanMR*100~fxavg_mr_all_v1$stdage,xlim=c(0,280),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main="Mortality - Version 1")
points(navg_mr_all_v1$meanMR*100~navg_mr_all_v1$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c("N Fixers","Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

plot(fxavg_mr_all_v2$meanMR*100~fxavg_mr_all_v2$stdage,xlim=c(0,280),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main="Mortality - Version 2")
points(navg_mr_all_v2$meanMR*100~navg_mr_all_v2$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c("N Fixers","Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

dev.off()

######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)

#setwd("")
write.csv(mdata_plotInfo,file="FIA_mdata_all_plotInfo_S0_150728.csv")

######################################################
##################    STAT Analysis   ################
######################################################

library(bbmle)
mdata<-mdata.original

################## Saturating

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	}

# fixers

# MLE Fit

pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
x<-0:300
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f_new

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit

pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
x<-0:300
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n_new

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

################## Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# fixers

# MLE Fit

pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}


fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit
x<-0:300
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f_new

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit

pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
x<-0:300
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n_new

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
#[1] 2.639268
#[1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
#[1] 2.665175
#[1] 0

#setwd("")
pv_f_mr_all<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_all<-coef(fit_ricker_surv_mi_d_t_n)
write.csv(pv_f_mr_all,file="pv_f_mr_all.csv")
write.csv(pv_n_mr_all,file="pv_n_mr_all.csv")

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

agevec_mr_all<-unique(mdata$stdage)

ndraws<-1000
xs<-seq(min(agevec_mr_all),max(agevec_mr_all),1)
xs_mr_all<-seq(min(agevec_mr_all),max(agevec_mr_all),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_all,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_all,Sigma=vcov(fit_ricker_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata$dbh0))),xs,exp(mean(log(mdata$dbh0))),1)*100
	ys_non[i,]<-fn_ricker_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata$dbh0))),xs,exp(mean(log(mdata$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_all<-civec_fix
civec_non_mr_all<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

mdata<-mdata.original

#setwd("")
pdf("Fig5_MRfit_wCI_version1_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_all_v1$meanMR*100~fxavg_mr_all_v1$stdage,xlim=c(0,250),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Mortality")
points(navg_mr_all_v1$meanMR*100~navg_mr_all_v1$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_all[1],pv_f_mr_all[2],pv_f_mr_all[3],pv_f_mr_all[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)*100,from=min(agevec_mr_all),to=max(agevec_mr_all),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_all[1],pv_n_mr_all[2],pv_n_mr_all[3],pv_n_mr_all[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)*100,from=min(agevec_mr_all),to=max(agevec_mr_all),lwd=3,col="blue",lty=1,add=TRUE)

### CI for fixers
points(civec_fix_mr_all[1,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_all[2,]~xs_mr_all,lwd=2,lty=2,typ="l",col="red")

### CI for non fixers
points(civec_non_mr_all[1,]~xs_mr_all,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_all[2,]~xs_mr_all,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("N Fixers","Non Fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

mdata<-mdata.original

#setwd("")
pdf("Fig5_MRfit_wCI_version2_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_all_v2$meanMR*100~fxavg_mr_all_v2$stdage,xlim=c(0,250),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Mortality")
points(navg_mr_all_v2$meanMR*100~navg_mr_all_v2$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_all[1],pv_f_mr_all[2],pv_f_mr_all[3],pv_f_mr_all[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)*100,from=min(agevec_mr_all),to=max(agevec_mr_all),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_all[1],pv_n_mr_all[2],pv_n_mr_all[3],pv_n_mr_all[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)*100,from=min(agevec_mr_all),to=max(agevec_mr_all),lwd=3,col="blue",lty=1,add=TRUE)

### CI for fixers
points(civec_fix_mr_all[1,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_all[2,]~xs_mr_all,lwd=2,lty=2,typ="l",col="red")

### CI for non fixers
points(civec_non_mr_all[1,]~xs_mr_all,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_all[2,]~xs_mr_all,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("N fixers","Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

######################################################
##########    Robinia Mortality Analysis  ############
######################################################

#setwd("")
mdata.original<-read.csv("FIAMRSTAT_morality_data_stat_final_150709.csv",header=T)
mdata.original<-mdata.original[,c(-1)]

mdata<-mdata.original

# Get the plots with robinia present

mdata_rob_pcn<-unique(mdata[mdata$spcd0==901,]$pcn) #Robinia pseudoacacia present
mdata_rob<-data.frame(pcn=mdata_rob_pcn,rob_present=1)
mdata_rob1<-merge(mdata_rob,mdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
mdata_rob2<-mdata_rob1[!is.na(mdata_rob1$rob_present),]

mdata_robinia<-mdata_rob2[(mdata_rob2$FIX==1 & mdata_rob2$spcd==901) | mdata_rob2$FIX==0,]
mdata<-mdata_robinia

#setwd("")
write.csv(mdata_robinia,file="FIA_MRdataforfit_Robinia_150723.csv")

######################################################
##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1
######################################################

#Now: Fixers

mdata<-mdata_robinia
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_robinia_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_robinia
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_robinia_v1<-favg
navg_mr<-favg

######################################################
##################  Mean Calculation  ################ Version 2: plot with S=0, omit
######################################################

#Now: Fixers

mdata<-mdata_robinia
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_robinia_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_robinia
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_robinia_v2<-favg
navg_mr<-favg


######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)
mdata_plotInfo

#setwd("")
write.csv(mdata_plotInfo,file="FIA_mdata_robinia_plotInfo_S0_150728.csv")

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# pdf("Fig4_MRmean_Robinia_150714.pdf",width=9,height=4.5)

par(mfrow=c(1,2))

plot(fxavg_mr_robinia_v1$meanMR*100~fxavg_mr_robinia_v1$stdage,xlim=c(0,200),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Robinia~pseudoacacia)))
points(navg_mr_robinia_v1$meanMR*100~navg_mr_robinia_v1$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c("N Fixers","Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

plot(fxavg_mr_robinia_v2$meanMR*100~fxavg_mr_robinia_v2$stdage,xlim=c(0,200),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Robinia~pseudoacacia)))
points(navg_mr_robinia_v2$meanMR*100~navg_mr_robinia_v2$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c("N Fixers","Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

# dev.off()

######################################################
##################    STAT Analysis   ################
######################################################

library(bbmle)
mdata<-mdata_robinia

################ Saturating

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
}

# fixers

# MLE Fit

pv<-c(-4.5,-10.2,  -3.2,-5.9,  10,2,  -0.1,-1.5) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f_new

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit

pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

################ Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-2.5,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}


fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-2.5,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit

pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
########    Model Comparisons - AIC   ################
######################################################
deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
# [1] 6.697303
# [1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
# [1] 0
# [1] 4.164615

#setwd("")
pv_f_mr_robinia<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_robinia<-coef(fit_sat_surv_mi_d_t_n)
write.csv(pv_f_mr_robinia,file="pv_f_mr_robinia_updated.csv")
write.csv(pv_n_mr_robinia,file="pv_n_mr_robinia_updated.csv")

######################################################
############    CI Computation   #####################
######################################################
library(MASS)
agevec_mr_robinia<-unique(mdata_robinia$stdage)

ndraws<-1000
xs<-seq(min(agevec_mr_robinia),max(agevec_mr_robinia),1)
xs_mr_robinia<-seq(min(agevec_mr_robinia),max(agevec_mr_robinia),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_robinia,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_robinia,Sigma=vcov(fit_sat_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata_robinia$dbh0))),xs,exp(mean(log(mdata_robinia$dbh0))),1)*100
	ys_non[i,]<-fn_sat_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata_robinia$dbh0))),xs,exp(mean(log(mdata_robinia$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_robinia<-civec_fix
civec_non_mr_robinia<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig6_MRfit_wCI_version1_Updated_150731.pdf",width=4.5,height=4.5)

plot(fxavg_mr_robinia_v1$meanMR*100~fxavg_mr_robinia_v1$stdage,xlim=c(0,200),ylim=c(-2,100),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (% dead/year)",main=expression(italic(Robinia~pseudoacacia)))
points(navg_mr_robinia_v1$meanMR*100~navg_mr_robinia_v1$stdage,xlim=c(0,200),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_robinia[1],pv_f_mr_robinia[2],pv_f_mr_robinia[3],pv_f_mr_robinia[4],exp(mean(log(mdata_robinia$dbh0))),x,exp(mean(log(mdata_robinia$dbh0))),1)*100,from=min(agevec_mr_robinia),to=max(agevec_mr_robinia),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_sat_mort_d_t(pv_n_mr_robinia[1],pv_n_mr_robinia[2],pv_n_mr_robinia[3],pv_n_mr_robinia[4],exp(mean(log(mdata_robinia$dbh0))),x,exp(mean(log(mdata_robinia$dbh0))),1)*100,from=min(agevec_mr_robinia),to=max(agevec_mr_robinia),lwd=3,col="blue",lty=1,add=TRUE)

### CI for fixers
points(civec_fix_mr_robinia[1,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="red",xlim=c(3,173))
points(civec_fix_mr_robinia[2,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="red",xlim=c(3,173))

### CI for non fixers
points(civec_non_mr_robinia[1,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="blue",xlim=c(3,173))
points(civec_non_mr_robinia[2,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="blue",xlim=c(3,173))

legend('topright',c(expression(italic(Robinia~pseudoacacia)),"Non Fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

#setwd("")
pdf("Fig6_MRfit_wCI_version2_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_robinia_v2$meanMR*100~fxavg_mr_robinia_v2$stdage,xlim=c(0,200),ylim=c(-2,20),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Robinia~pseudoacacia)))
points(navg_mr_robinia_v2$meanMR*100~navg_mr_robinia_v2$stdage,xlim=c(0,200),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_robinia[1],pv_f_mr_robinia[2],pv_f_mr_robinia[3],pv_f_mr_robinia[4],exp(mean(log(mdata_robinia$dbh0))),x,exp(mean(log(mdata_robinia$dbh0))),1)*100,from=min(agevec_mr_robinia),to=max(agevec_mr_robinia),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_sat_mort_d_t(pv_n_mr_robinia[1],pv_n_mr_robinia[2],pv_n_mr_robinia[3],pv_n_mr_robinia[4],exp(mean(log(mdata_robinia$dbh0))),x,exp(mean(log(mdata_robinia$dbh0))),1)*100,from=min(agevec_mr_robinia),to=max(agevec_mr_robinia),lwd=3,col="blue",lty=1,add=TRUE)

### CI for fixers
points(civec_fix_mr_robinia[1,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="red",xlim=c(3,173))
points(civec_fix_mr_robinia[2,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="red",xlim=c(3,173))

### CI for non fixers
points(civec_non_mr_robinia[1,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="blue",xlim=c(3,173))
points(civec_non_mr_robinia[2,]~xs_mr_robinia,lwd=2,lty=2,typ="l",col="blue",xlim=c(3,173))

legend('topright',c(expression(italic(Robinia~pseudoacacia)),"Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

######################################################
########  Cercocarpus Mortality Analysis  ############
######################################################

mdata<-mdata.original

# Get the plots with cercocarpus present

mdata_crc_pcn<-unique(mdata[mdata$spcd==475,]$pcn) #Cercocarpus ledifolius present
mdata_crc<-data.frame(pcn=mdata_crc_pcn,crc_present=1)
mdata_crc1<-merge(mdata_crc,mdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
mdata_crc2<-mdata_crc1[!is.na(mdata_crc1$crc_present),]

mdata_cercocarpus<-mdata_crc2[(mdata_crc2$FIX==1 & mdata_crc2$spcd==475) | mdata_crc2$FIX==0,]
mdata<-mdata_cercocarpus

#setwd("")
write.csv(mdata_cercocarpus,file="FIAMR_mdata_cercocarpus_dataforfit_150729.csv")

######################################################
##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1
######################################################

#Now: Fixers

mdata<-mdata_cercocarpus
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_cercocarpus_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_cercocarpus
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_cercocarpus_v1<-favg
navg_mr<-favg

######################################################
##################  Mean Calculation  ################ Version 2: plot with S=0, omit
######################################################

#Now: Fixers

mdata<-mdata_cercocarpus
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_cercocarpus_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_cercocarpus
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_cercocarpus_v2<-favg
navg_mr<-favg

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig5_MRmean_Cercocarpus_150714.pdf",width=9,height=4.5)

par(mfrow=c(1,2))

plot(fxavg_mr_cercocarpus_v1$meanMR*100~fxavg_mr_cercocarpus_v1$stdage,xlim=c(0,250),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Cercocarpus~ledifolius)))
points(navg_mr_cercocarpus_v1$meanMR*100~navg_mr_cercocarpus_v1$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c(expression(italic(Cercocarpus~ledifolius)),"Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

plot(fxavg_mr_cercocarpus_v2$meanMR*100~fxavg_mr_cercocarpus_v2$stdage,xlim=c(0,250),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Cercocarpus~ledifolius)))
points(navg_mr_cercocarpus_v2$meanMR*100~navg_mr_cercocarpus_v2$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c(expression(italic(Cercocarpus~ledifolius)),"Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

dev.off()

######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)
mdata_plotInfo

#setwd("")
write.csv(mdata_plotInfo,file="FIA_mdata_cercocarpus_plotInfo_S0_150728.csv")

######################################################
##################    STAT Analysis   ################
######################################################

mdata<-mdata_cercocarpus

################ Saturating

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	}

# MLE Fit

pv<-c(-4.5,-5.2,-3.2,-6.0,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-4.5,-5.2,-3.2,-6.0,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,1))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

################ Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# MLE Fit

pv<-c(-5.6,-5.6, -2.9,-7.3,  -0.1,-0.1,  0.08,0.1) #7.3
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}


fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit

pv<-c(-5.6,-5.6, -2.9,-7.3,  -0.1,-0.1,  0.08,0.1) #7.3
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit

pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
# [1] 16.96953
# [1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
# [1] 6.854258
# [1] 0

#setwd("")
pv_f_mr_cercocarpus<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_cercocarpus<-coef(fit_ricker_surv_mi_d_t_n)
write.csv(pv_f_mr_cercocarpus,file="pv_f_mr_cercocarpus_updated.csv")
write.csv(pv_n_mr_cercocarpus,file="pv_n_mr_cercocarpus_updated.csv")

######################################################
############    CI Computation   #####################
######################################################

agevec_mr_cercocarpus<-unique(mdata_cercocarpus$stdage)

ndraws<-1000
xs_mr_cercocarpus<-seq(min(agevec_mr_cercocarpus),max(agevec_mr_cercocarpus),1)
xs<-seq(min(agevec_mr_cercocarpus),max(agevec_mr_cercocarpus),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_cercocarpus,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_cercocarpus,Sigma=vcov(fit_ricker_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata_cercocarpus$dbh0))),xs,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100
	ys_non[i,]<-fn_ricker_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata_cercocarpus$dbh0))),xs,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_cercocarpus<-civec_fix
civec_non_mr_cercocarpus<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig7_MRfit_wCI_Cercocarpus_version1_Updated_170728.pdf",width=4.5,height=4.5)

plot(fxavg_mr_cercocarpus_v1$meanMR*100~fxavg_mr_cercocarpus_v1$stdage,xlim=c(min(agevec_mr_cercocarpus),max(agevec_mr_cercocarpus)),col="red",cex=1,pch=1,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Cercocarpus~ledifolius)),ylim=c(0,100))
points(navg_mr_cercocarpus_v1$meanMR*100~navg_mr_cercocarpus_v1$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_cercocarpus[1],pv_f_mr_cercocarpus[2],pv_f_mr_cercocarpus[3],pv_f_mr_cercocarpus[4],exp(mean(log(mdata_cercocarpus$dbh0))),x,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_cercocarpus[1],pv_n_mr_cercocarpus[2],pv_n_mr_cercocarpus[3],pv_n_mr_cercocarpus[4],exp(mean(log(mdata_cercocarpus$dbh0))),x,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100,from=min(agevec_mr_cercocarpus),to=max(agevec_mr_cercocarpus),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_cercocarpus[1,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_cercocarpus[2,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_cercocarpus[1,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_cercocarpus[2,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c(expression(italic(Cercocarpus~ledifolius)),"Non Fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

#setwd("")
pdf("Fig7_MRfit_wCI_Cercocarpus_version2_Updated_170928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_cercocarpus_v2$meanMR*100~fxavg_mr_cercocarpus_v2$stdage,xlim=c(min(agevec_mr_cercocarpus),max(agevec_mr_cercocarpus)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Cercocarpus~ledifolius)),ylim=c(0,15))
points(navg_mr_cercocarpus_v2$meanMR*100~navg_mr_cercocarpus_v2$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_cercocarpus[1],pv_f_mr_cercocarpus[2],pv_f_mr_cercocarpus[3],pv_f_mr_cercocarpus[4],exp(mean(log(mdata_cercocarpus$dbh0))),x,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_cercocarpus[1],pv_n_mr_cercocarpus[2],pv_n_mr_cercocarpus[3],pv_n_mr_cercocarpus[4],exp(mean(log(mdata_cercocarpus$dbh0))),x,exp(mean(log(mdata_cercocarpus$dbh0))),1)*100,from=min(agevec_mr_cercocarpus),to=max(agevec_mr_cercocarpus),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_cercocarpus[1,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_cercocarpus[2,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_cercocarpus[1,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_cercocarpus[2,]~xs_mr_cercocarpus,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c(expression(italic(Cercocarpus~ledifolius)),"Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

######################################################
###########  Alnus Mortality Analysis  ###############
######################################################

mdata<-mdata.original

# Get the plots with Alnus rubrum present

mdata_alnr_pcn<-unique(mdata[mdata$spcd0==351,]$pcn) #Alnus rubrum present
mdata_alnr<-data.frame(pcn=mdata_alnr_pcn,alnr_present=1)
mdata_alnr1<-merge(mdata_alnr,mdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
mdata_alnr2<-mdata_alnr1[!is.na(mdata_alnr1$alnr_present),]

mdata_alnus<-mdata_alnr2[(mdata_alnr2$FIX==1 & mdata_alnr2$spcd==351) | mdata_alnr2$FIX==0,]
mdata<-mdata_alnus

#setwd("")
write.csv(mdata_alnus,file="FIA_MR_mdata_alnus_dataforfit_150729.csv")

######################################################
##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1
######################################################

#Now: Fixers

mdata<-mdata_alnus
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_alnus_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_alnus
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_alnus_v1<-favg
navg_mr<-favg

######################################################
##################  Mean Calculation  ################ Version 2: plot with S=0, omit
######################################################

#Now: Fixers

mdata<-mdata_alnus
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_alnus_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_alnus
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_alnus_v2<-favg
navg_mr<-favg

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig6_MRmean_Alnus_150715.pdf",width=9,height=4.5)

par(mfrow=c(1,2))

plot(fxavg_mr_alnus_v1$meanMR*100~fxavg_mr_alnus_v1$stdage,xlim=c(0,200),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Alnus~rubra)))
points(navg_mr_alnus_v1$meanMR*100~navg_mr_alnus_v1$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c(expression(italic(Alnus~rubra)),"Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

plot(fxavg_mr_alnus_v2$meanMR*100~fxavg_mr_alnus_v2$stdage,xlim=c(0,200),col="red",cex=0.8,pch=19,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Alnus~rubra)))
points(navg_mr_alnus_v2$meanMR*100~navg_mr_alnus_v2$stdage,xlim=c(0,280),col="blue",cex=0.8,pch=17)
legend('topright',c(expression(italic(Alnus~rubra)),"Non Fixers"),col=c("red","blue"),bty="n",pch=c(19,17))

dev.off()

######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)
mdata_plotInfo

#setwd("")
write.csv(mdata_plotInfo,file="FIA_mdata_alnus_plotInfo_S0_150728.csv")

######################################################
##################    STAT Analysis   ################
######################################################

mdata<-mdata.original

################ Saturating

# fixers

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	}

# MLE Fit

x<-0:300
pv<-c(-4.5,-5.2,-3.2,-6.0,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-4.5,-5.2,-3.2,-6.0,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

################ Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# fixers

# MLE Fit
pv<-c(-5.6,-6.6,-2.9,-3.9,-0.1,-0.1,0.08,0.08)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}


fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-5.6,-6.6,-2.9,-3.9,-0.1,-0.1,0.08,0.08)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
########    Model Comparisons - AIC   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
# [1] 42.22902
# [1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
# [1] 2.665175
# [1] 0

#setwd("")
pv_f_mr_alnus<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_alnus<-coef(fit_ricker_surv_mi_d_t_n)
write.csv(pv_f_mr_alnus,file="pv_f_mr_alnus_updated.csv")
write.csv(pv_n_mr_alnus,file="pv_n_mr_alnus_updated.csv")

######################################################
############    CI Computation   #####################
######################################################

agevec_mr_alnus<-unique(mdata_alnus$stdage)

ndraws<-10000
xs_mr_alnus<-seq(min(agevec_mr_alnus),max(agevec_mr_alnus),1)
xs<-seq(min(agevec_mr_alnus),max(agevec_mr_alnus),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_alnus,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_alnus,Sigma=vcov(fit_ricker_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata_alnus$dbh0))),xs,exp(mean(log(mdata_alnus$dbh0))),1)*100
	ys_non[i,]<-fn_ricker_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata_alnus$dbh0))),xs,exp(mean(log(mdata_alnus$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_alnus<-civec_fix
civec_non_mr_alnus<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig8_MRfit_wCI_aLNUS_version1_Updated_150728.pdf",width=4.5,height=4.5)

plot(fxavg_mr_alnus_v1$meanMR*100~fxavg_mr_alnus_v1$stdage,xlim=c(min(agevec_mr_alnus),max(agevec_mr_alnus)),col="red",cex=1,pch=1,xlab="Stand Age (Year)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Alnus~rubra)),ylim=c(0,100))
points(navg_mr_alnus_v1$meanMR*100~navg_mr_alnus_v1$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_alnus[1],pv_f_mr_alnus[2],pv_f_mr_alnus[3],pv_f_mr_alnus[4],exp(mean(log(mdata_alnus$dbh0))),x,exp(mean(log(mdata_alnus$dbh0))),1)*100,from=min(agevec_mr_alnus),to=max(agevec_mr_alnus),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_alnus[1],pv_n_mr_alnus[2],pv_n_mr_alnus[3],pv_n_mr_alnus[4],exp(mean(log(mdata_alnus$dbh0))),x,exp(mean(log(mdata_alnus$dbh0))),1)*100,from=min(agevec_mr_alnus),to=max(agevec_mr_alnus),lwd=3,col="blue",lty=1,add=TRUE)

legend('topright',c(expression(italic(Alnus~rubra)),"Non Fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

points(civec_fix[1,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_fix[2,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_non[1,]~xs,lwd=2,lty=2,typ="l",col="blue")
points(civec_non[2,]~xs,lwd=2,lty=2,typ="l",col="blue")

dev.off()

#setwd("")
pdf("Fig8_MRfit_wCI_aLNUS_version2_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_alnus_v2$meanMR*100~fxavg_mr_alnus_v2$stdage,xlim=c(min(agevec_mr_alnus),max(agevec_mr_alnus)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main=expression(italic(Alnus~rubra)),ylim=c(0,15))
points(navg_mr_alnus_v2$meanMR*100~navg_mr_alnus_v2$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_alnus[1],pv_f_mr_alnus[2],pv_f_mr_alnus[3],pv_f_mr_alnus[4],exp(mean(log(mdata_alnus$dbh0))),x,exp(mean(log(mdata_alnus$dbh0))),1)*100,from=min(agevec_mr_alnus),to=max(agevec_mr_alnus),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_alnus[1],pv_n_mr_alnus[2],pv_n_mr_alnus[3],pv_n_mr_alnus[4],exp(mean(log(mdata_alnus$dbh0))),x,exp(mean(log(mdata_alnus$dbh0))),1)*100,from=min(agevec_mr_alnus),to=max(agevec_mr_alnus),lwd=3,col="blue",lty=1,add=TRUE)

legend('topright',c(expression(italic(Alnus~rubra)),"Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

points(civec_fix_mr_alnus[1,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_alnus[2,]~xs,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_alnus[1,]~xs,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_alnus[2,]~xs,lwd=2,lty=2,typ="l",col="blue")

dev.off()




#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################




####################################################################################################
################################         Recruitment Analysis         ##############################
####################################################################################################

#setwd("")
t.ts<-read.csv("tree_ts.csv",header=T)
tcn1<-t.ts$tcn1[!is.na(t.ts$tcn1)]
tcn1<-unique(tcn1)
tcn<-tcn1
t_1m<-data.frame(tcn=tcn,msID=1) #tree_first measurement # Here: I got the trees that are measured for the first time.

p.ts<-read.csv("plot_ts.csv",header=T)
p_rec<-data.frame(pcn=unique(c(p.ts$pcn2,p.ts$pcn3,p.ts$pcn4)),pID=0) #not first measurement plot=0
p_rec<-na.omit(p_rec) # Plots that are not measured for the first time.

t<-read.csv("tree_small.csv",header=T)
ft<-t[,c(2,4)]
rm(t)

tp<-merge(ft,t_1m,by.x="tcn",by.y="tcn",all.x=FALSE,all.y=FALSE) # This gives me the information of first time measurement trees and their corresponding plot information

mtp<-merge(tp,p_rec,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE) 
# This gives me the trees that are measured for the first time, and in the plots that are not measured for the first time --> recruited (not present for the first plot measurement)

#setwd("")
write.csv(mtp,file="FIARRSTAT_recruited_tree_150715.csv")

############# Step2: Assign status for each individual tree (N,S,R)

#setwd("")
mtts<-read.csv("FIAMRSTAT_merged_ts_tree_status_150709.csv",header=T)
mtts<-mtts[,c(-1)]

rec<-data.frame(tcn=mtp$tcn,R=1)
mR<-merge(rec,mtts,by.x="tcn",by.y="tcn0",all.x=TRUE,all.y=TRUE) # This assigns the ones that are recruited as R=1
mR$N<-NA
mR$S<-NA
mR[!is.na(mR$statuscd0)& mR$statuscd0==1,]$N<-1 # Alive
mR[!is.na(mR$statuscd0)&!is.na(mR$statuscdT)&mR$statuscd0==1&mR$statuscdT==2,]$S<-0
mR[mR$statuscd0==1&mR$statuscdT==1&!is.na(mR$statuscd0)&!is.na(mR$statuscdT),]$S<-1

#setwd("")
write.csv(mR,file="FIARRSTAT_merged_ts_tree_NS_150715.csv")

############ Step3: Assign fixing status to species and trees

#setwd("")
p<-read.csv("plot.csv",header=T)
fp<-p[p$natforsamp==1&p$propsamp==1&p$ntlog==0,]
ffp<-fp[,c(2,8)]
mdat<-merge(ffp,mR,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

#setwd("")
fixer<-read.csv("FIA_fixerspcd_merged_noCeltis.csv",header=T) 
rdata<-merge(mdat,fixer,by.x="spcd0",by.y="SPCD",all.x=TRUE,all.y=FALSE)
rdata[is.na(rdata$FIX),]$FIX<-0

#setwd("")
write.csv(rdata,file="FIARRSTAT_recruitment_data_stat_final_150715.csv")

############ Step4: Get the plots with fixer present

rdata<-rdata[rdata$stdage>=0&!is.na(rdata$stdage),]
rdata_fp_pcn<-unique(rdata[rdata$FIX==1,]$pcn) #fixer present
rdata_fp<-data.frame(pcn=rdata_fp_pcn,fixer_present=1)
rdata_fp1<-merge(rdata_fp,rdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
rdata_fp2<-rdata_fp1[!is.na(rdata_fp1$fixer_present),]

rdata.original<-rdata
rdata<-rdata_fp2

#setwd("")
write.csv(rdata_fp2,file="FIARRSTAT_recruitment_data_stat_final_fixerpresent_150715.csv")

########## Here: print out the no. of individuals for each species

rdata$count<-1
spcd_rr<-tapply(rdata[rdata$FIX==1,]$count,rdata[rdata$FIX==1,]$spcd0,sum)
spcd_RR<-data.frame(spcd=rownames(spcd_rr),count=spcd_rr)
spcd_RR<-merge(spcd_RR,fixer,by.x="spcd",by.y="SPCD",all.x=FALSE,all.y=FALSE)
spcd_RR<-spcd_RR[order(-spcd_RR$count),]

#setwd("")
write.csv(spcd_RR,file="FIARRSTAT_FixerInfo_150715.csv")

############ Step5: Now compute for number of trees alive at time0, and number of trees recruited at timeT

#setwd("")
pts<-read.csv("plot_ts.csv",header=T)
p1<-data.frame(pcn0=pts$pcn1,pcnT=pts$pcn2)
p2<-data.frame(pcn0=pts$pcn2,pcnT=pts$pcn3)
p3<-data.frame(pcn0=pts$pcn3,pcnT=pts$pcn4)
pcn<-rbind(p1,p2,p3)
pcn<-pcn[!is.na(pcn$pcn0)&!is.na(pcn$pcnT),]

pcn0<-data.frame(pcn0=pcn$pcn0)
pcnT<-data.frame(pcnT=pcn$pcnT)

rdata0<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcn0,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=FALSE)
colnames(rdata0)<-c("pcn0","spcd_0","stdage_0","tcn_0","R_0","dbh_0","N_0","FIX_0")
rdataT<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcnT,by.x="pcn",by.y="pcnT",all.x=FALSE,all.y=FALSE)
colnames(rdataT)<-c("pcnT","spcd_T","stdage_T","tcn_T","R_T","dbh_T","N_T","FIX_T")

sum.na<-function(x){
	x<-x[!is.na(x)]
	sum(x)
}

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

rdata.original<-rdata

# For pcn0: time 0

rdata<-rdata0

nonN<-tapply(rdata[rdata$FIX_0==0,]$N_0,rdata[rdata$FIX_0==0,]$pcn0,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_0==0,]$dbh_0,rdata[rdata$FIX_0==0,]$pcn0,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na) # for time0, calculate mean diameter for all trees

fixN<-tapply(rdata[rdata$FIX_0==1,]$N_0,rdata[rdata$FIX_0==1,]$pcn0,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_0==1,]$dbh_0,rdata[rdata$FIX_0==1,]$pcn0,mean.na)
fixmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na)

n1<-data.frame(pcn0=as.numeric(rownames(nonN)),N=nonN)
n4<-data.frame(pcn0=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcn0=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcn0=as.numeric(rownames(fixN)),N=fixN)
f4<-data.frame(pcn0=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcn0=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)

fix_0<-mf2
fix_0$FIX<-1
non_0<-mn2
non_0$FIX<-0
final_0<-rbind(fix_0,non_0)
final_0<-final_0[!is.na(final_0$N)&!is.na(final_0$pcn0)&!is.na(final_0$mDBH)&!is.na(final_0$mDBHbt)&!is.na(final_0$FIX),]

# For pcnT: time T

rdata<-rdataT

nonR<-tapply(rdata[rdata$FIX_T==0,]$R_T,rdata[rdata$FIX_T==0,]$pcnT,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_T==0,]$dbh_T,rdata[rdata$FIX_T==0,]$pcnT,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na) # for time0, calculate mean diameter for all trees

fixR<-tapply(rdata[rdata$FIX_T==1,]$R_T,rdata[rdata$FIX_T==1,]$pcnT,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_T==1,]$dbh_T,rdata[rdata$FIX_T==1,]$pcnT,mean.na)
fixmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na)

n1<-data.frame(pcnT=as.numeric(rownames(nonR)),R=nonR)
n4<-data.frame(pcnT=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcnT=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcnT=as.numeric(rownames(fixR)),R=fixR)
f4<-data.frame(pcnT=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcnT=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)

fix_T<-mf2
fix_T$FIX<-1
non_T<-mn2
non_T$FIX<-0
final_T<-rbind(fix_T,non_T)
final_T<-final_T[!is.na(final_T$R)&!is.na(final_T$pcnT)&!is.na(final_T$mDBH)&!is.na(final_T$mDBHbt)&!is.na(final_T$FIX),]

pcn1<-merge(pcn,final_0,by.x="pcn0",by.y="pcn0",all.x=FALSE,all.y=FALSE)
pcn2<-merge(pcn1,final_T[,c(1,2,5)],by=c("pcnT","FIX"),all.x=FALSE,all.y=FALSE)

final<-pcn2

fp<-p[,c(1,2,5,6,8)]
MRdat<-merge(fp,final,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

mt<-data.frame(pcn=p$pcn,mt=p$meastime)
MRdat_mt0<-merge(MRdat,mt,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mt0)[12]<-"mt0"
MRdat_mtT<-merge(MRdat_mt0,mt,by.x="pcnT",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mtT)[13]<-"mtT"

MRdat_mtT$mt<-MRdat_mtT$mtT-MRdat_mtT$mt0
MRdat<-MRdat_mtT
MRdat$RR<-NA
MRdat<-MRdat[MRdat$N>0,]
MRdat$RR<-MRdat$R/(MRdat$N*MRdat$mt)

MRdat<-MRdat[MRdat$stdage>=0,]

#setwd("")
write.csv(MRdat,"FIARR_Recruitment_data_for_fit_final_fixerpresent_150715.csv")

######################################################
##############    Mean Calculations   ################
######################################################

avg<-tapply(MRdat[MRdat$FIX==1,]$RR,MRdat[MRdat$FIX==1,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
fxavg_rr<-favg
fxavg_rr_all<-favg
plot(fxavg_rr$meanRR~fxavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Fixers")

avg<-tapply(MRdat[MRdat$FIX==0,]$RR,MRdat[MRdat$FIX==0,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
navg_rr<-favg
navg_rr_all<-favg
plot(navg_rr$meanRR~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Non Fixers")

plot(navg_rr$meanRR*100~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",lwd=1,col="blue",main="U.S. Trees",pch=2,ylim=c(0,100))
points(fxavg_rr$meanRR*100~fxavg_rr$stdage,lwd=1,col="red",pch=1)
legend("topright",c("N Fixers","Non Fixers"),col=c("red","blue"),pch=c(1,2),bty="n")

######################################################
##################    STAT Analysis   ################
######################################################

MRdat.original<-MRdat

################ Saturating

fn1_rec_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
	y <- y*t
} # Saturating with diameter

# fixers

# MLE
pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_f<-pv

NLL2_rec_d_t_f <- function(af,bf,cf,df,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	yf <- fn1_rec_d_t(af, bf, cf, df, D[F==1], A[F==1], Dhat, t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit2_rec_d_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_f)

fit2_rec_d_c_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_f)

# non-fixers

pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_n<-pv

NLL2_rec_d_t_n <- function(a0,b0,c0,d0,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	y0 <- fn1_rec_d_t(a0, b0, c0, d0, D[F==0], A[F==0], Dhat, t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit2_rec_d_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_n)

fit2_rec_d_c_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_n)

################ Ricker

fn4_rec_d_t<-function(a,c,d,k,D,A,Dhat,t){
		a<-exp(a)
	c<-exp(c)
y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	y<-y*t
	} #Right-skewed with diameter effect

# fixers

# MLE

pv<-c(-4,-8.6,-3.7,-1.6,-0.6,-0.2,0.2,0.18)
pv_f<-pv

NLL5_rec_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	yf<-fn4_rec_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit5_rec_d_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_f)

fit5_rec_d_c_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_f)

# non-fixers

# MLE

NLL5_rec_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	y0<-fn4_rec_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit5_rec_d_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_n)

fit5_rec_d_c_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_n)

######################################################
################   Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

# fixer
deltaAICs(c(AIC(fit2_rec_d_t_f),AIC(fit2_rec_d_c_t_f),AIC(fit5_rec_d_t_f),AIC(fit5_rec_d_c_t_f)))

# [1] 10.69596
# [1] 0.2586587
# [1] 11.41815
# [1] 0

# non-fixers
deltaAICs(c(AIC(fit2_rec_d_t_n),AIC(fit2_rec_d_c_t_n),AIC(fit5_rec_d_t_n),AIC(fit5_rec_d_c_t_n)))
# [1] 0.08704712
# [1] 0 
# [1] 215.1468
# [1] 209.3612

#setwd("")
pv_f_rr_all<-coef(fit5_rec_d_c_t_f)
pv_n_rr_all<-coef(fit2_rec_d_c_t_n)
write.csv(pv_f_rr_all,file="pv_f_rr_all.csv")
write.csv(pv_n_rr_all,file="pv_n_rr_all.csv")

######################################################
#################    CI Computation   ################
######################################################

agevec_rr_all<-unique(MRdat$stdage)
mD<-mean(MRdat$mDBH)
mD_fix<-mean(MRdat[MRdat$FIX==1,]$mDBH)
mD_non<-mean(MRdat[MRdat$FIX==0,]$mDBH)

ndraws<-10000
xs_rr_all<-seq(min(agevec_rr_all),max(agevec_rr_all),1)
xs<-seq(min(agevec_rr_all),max(agevec_rr_all),1)

vmat_fix_rr<-mvrnorm(ndraws,mu=pv_f_rr_all,Sigma=vcov(fit5_rec_d_c_t_f)) 
vmat_non_rr<-mvrnorm(ndraws,mu=pv_n_rr_all,Sigma=vcov(fit2_rec_d_c_t_n))

ys_fix_rr<-array(dim=c(ndraws,length(xs)))
ys_non_rr<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix_rr[i,]<-fn4_rec_d_t(vmat_fix_rr[i,1],vmat_fix_rr[i,2],vmat_fix_rr[i,3],vmat_fix_rr[i,4],mD_fix,xs,mD,1)*100
	ys_non_rr[i,]<-fn1_rec_d_t(vmat_non_rr[i,1],vmat_non_rr[i,2],vmat_non_rr[i,3],vmat_non_rr[i,4],mD_non,xs,mD,1)*100
}

civec_fix_rr<-array(dim=c(2,length(xs)))
civec_non_rr<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix_rr[,j]<-quantile(ys_fix_rr[,j],c(lowb,upb),na.rm=TRUE)
	civec_non_rr[,j]<-quantile(ys_non_rr[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_rr_all<-civec_fix_rr
civec_non_rr_all<-civec_non_rr

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig9_RRfit_wCI_150928.pdf",width=4.5,height=4.5)

x<-c(0,180)
y<-c(0,100)
plot(y~x,type="n",xlim=c(0,180),y=c(0,100),xlab="Stand Age (Years)",ylab="Recruitment Rate (% Recruited/year)",main="Recruitment")
points(navg_rr$meanRR*100~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",col="blue",main="U.S. Trees",pch=2,cex=0.8)
points(fxavg_rr$meanRR*100~fxavg_rr$stdage,col="red",pch=1,cex=0.8)

# fixers
	# Ricker
curve(fn4_rec_d_t(pv_f_rr_all[1],pv_f_rr_all[2],pv_f_rr_all[3],pv_f_rr_all[4],mD_fix,x,mD,1)*100,from=min(agevec_rr_all),to=max(agevec_rr_all),lwd=3,lty=1,add=TRUE,col="red")

# non-fixers
	# Saturating
curve(fn1_rec_d_t(pv_n_rr_all[1],pv_n_rr_all[2],pv_n_rr_all[3],pv_n_rr_all[4],mD_non,x,mD,1)*100,from=min(agevec_rr_all),to=max(agevec_rr_all),lwd=3,col="blue",lty=1,add=TRUE,xlab="Stand Age",ylab="Estimated Relative Recruitment Rate (%)")

legend("topright",c("N fixers","Non-fixers"),lty=1,bty="n",cex=1.2,lwd=2,col=c("red","blue"),pch=c(1,2))

points(civec_fix_rr_all[1,]~xs_rr_all,lwd=2,typ="l",lty=2,col="red")
points(civec_fix_rr_all[2,]~xs_rr_all,lwd=2,typ="l",lty=2,col="red")

points(civec_non_rr_all[1,]~xs_rr_all,lwd=2,typ="l",lty=2,col="blue") # This CI is really wierd. Need to troubleshoot.
points(civec_non_rr_all[2,]~xs_rr_all,lwd=2,typ="l",lty=2,col="blue")

dev.off()

######################################################
#########    Robinia Recruitment Analysis  ###########
######################################################

# rdata_fp2 ---> rdata with fixer presents (plots)

#setwd("")
rdata_fp2<-read.csv("FIARRSTAT_recruitment_data_stat_final_fixerpresent_150715.csv",header=T)
rdata_fp2<-rdata_fp2[,c(-1)]

rdata<-rdata_fp2

rdata_rob_pcn<-unique(rdata[rdata$spcd==901,]$pcn) #Robinia pseudoacacia present
rdata_rob<-data.frame(pcn=rdata_rob_pcn,rob_present=1)
rdata_rob1<-merge(rdata_rob,rdata,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)
rdata_rob2<-rdata_rob1[!is.na(rdata_rob1$rob_present),]

rdata_robinia<-rdata_rob2[(rdata_rob2$FIX==1 & rdata_rob2$spcd==901) | rdata_rob2$FIX==0,]
rdata<-rdata_robinia

#setwd("")
pts<-read.csv("plot_ts.csv",header=T)
p1<-data.frame(pcn0=pts$pcn1,pcnT=pts$pcn2)
p2<-data.frame(pcn0=pts$pcn2,pcnT=pts$pcn3)
p3<-data.frame(pcn0=pts$pcn3,pcnT=pts$pcn4)
pcn<-rbind(p1,p2,p3)
pcn<-pcn[!is.na(pcn$pcn0)&!is.na(pcn$pcnT),]

pcn0<-data.frame(pcn0=pcn$pcn0)
pcnT<-data.frame(pcnT=pcn$pcnT)

rdata0<-merge(rdata_robinia[,c(1,4,5,6,7,12,14,20)],pcn0,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=FALSE)
colnames(rdata0)<-c("pcn0","spcd_0","stdage_0","tcn_0","R_0","dbh_0","N_0","FIX_0")
rdataT<-merge(rdata_robinia[,c(1,4,5,6,7,12,14,20)],pcnT,by.x="pcn",by.y="pcnT",all.x=FALSE,all.y=FALSE)
colnames(rdataT)<-c("pcnT","spcd_T","stdage_T","tcn_T","R_T","dbh_T","N_T","FIX_T")

sum.na<-function(x){
	x<-x[!is.na(x)]
	sum(x)
}

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

rdata.original<-rdata
# For pcn0: time 0

rdata<-rdata0

nonN<-tapply(rdata[rdata$FIX_0==0,]$N_0,rdata[rdata$FIX_0==0,]$pcn0,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_0==0,]$dbh_0,rdata[rdata$FIX_0==0,]$pcn0,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na) # for time0, calculate mean diameter for all trees

fixN<-tapply(rdata[rdata$FIX_0==1,]$N_0,rdata[rdata$FIX_0==1,]$pcn0,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_0==1,]$dbh_0,rdata[rdata$FIX_0==1,]$pcn0,mean.na)
fixmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na)

n1<-data.frame(pcn0=as.numeric(rownames(nonN)),N=nonN)
n4<-data.frame(pcn0=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcn0=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcn0=as.numeric(rownames(fixN)),N=fixN)
f4<-data.frame(pcn0=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcn0=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)

fix_0<-mf2
fix_0$FIX<-1
non_0<-mn2
non_0$FIX<-0
final_0<-rbind(fix_0,non_0)
final_0<-final_0[!is.na(final_0$N)&!is.na(final_0$pcn0)&!is.na(final_0$mDBH)&!is.na(final_0$mDBHbt)&!is.na(final_0$FIX),]

# For pcnT: time T

rdata<-rdataT

nonR<-tapply(rdata[rdata$FIX_T==0,]$R_T,rdata[rdata$FIX_T==0,]$pcnT,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_T==0,]$dbh_T,rdata[rdata$FIX_T==0,]$pcnT,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na) # for time0, calculate mean diameter for all trees

fixR<-tapply(rdata[rdata$FIX_T==1,]$R_T,rdata[rdata$FIX_T==1,]$pcnT,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_T==1,]$dbh_T,rdata[rdata$FIX_T==1,]$pcnT,mean.na)
fixmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na)

n1<-data.frame(pcnT=as.numeric(rownames(nonR)),R=nonR)
n4<-data.frame(pcnT=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcnT=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcnT=as.numeric(rownames(fixR)),R=fixR)
f4<-data.frame(pcnT=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcnT=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)

fix_T<-mf2
fix_T$FIX<-1
non_T<-mn2
non_T$FIX<-0
final_T<-rbind(fix_T,non_T)
final_T<-final_T[!is.na(final_T$R)&!is.na(final_T$pcnT)&!is.na(final_T$mDBH)&!is.na(final_T$mDBHbt)&!is.na(final_T$FIX),]

pcn1<-merge(pcn,final_0,by.x="pcn0",by.y="pcn0",all.x=FALSE,all.y=FALSE)
pcn2<-merge(pcn1,final_T[,c(1,2,5)],by=c("pcnT","FIX"),all.x=FALSE,all.y=FALSE)

final<-pcn2

#setwd("")

p<-read.csv("plot.csv",header=T) # Load in plot data
fp<-p[,c(1,2,5,6,8)]
MRdat<-merge(fp,final,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

mt<-data.frame(pcn=p$pcn,mt=p$meastime)
MRdat_mt0<-merge(MRdat,mt,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mt0)[12]<-"mt0"
MRdat_mtT<-merge(MRdat_mt0,mt,by.x="pcnT",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mtT)[13]<-"mtT"

MRdat_mtT$mt<-MRdat_mtT$mtT-MRdat_mtT$mt0
MRdat<-MRdat_mtT
MRdat$RR<-NA
MRdat<-MRdat[MRdat$N>0,]
MRdat$RR<-MRdat$R/(MRdat$N*MRdat$mt)

MRdat_robinia<-MRdat

#setwd("")
write.csv(MRdat_robinia,file="FIA_RRdataforfit_Robinia_150723.csv")

######################################################
###############    Mean Calculations   ###############
######################################################

avg<-tapply(MRdat[MRdat$FIX==1,]$RR,MRdat[MRdat$FIX==1,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
fxavg_rr<-favg
fxavg_rr_robinia<-favg
plot(fxavg_rr$meanRR~fxavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Fixers")

avg<-tapply(MRdat[MRdat$FIX==0,]$RR,MRdat[MRdat$FIX==0,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
navg_rr<-favg
navg_rr_robinia<-favg
plot(navg_rr$meanRR~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Non Fixers")

plot(navg_rr$meanRR*100~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",lwd=1,col="blue",main="U.S. Trees",pch=2,ylim=c(0,100))
points(fxavg_rr$meanRR*100~fxavg_rr$stdage,lwd=1,col="red",pch=1)
legend("topright",c("N Fixers","Non Fixers"),col=c("red","blue"),pch=c(1,2),bty="n")

######################################################
##################    STAT Analysis   ################
######################################################

MRdat.original<-MRdat

################ Saturating

fn1_rec_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
	y <- y*t
} # Saturating with diameter

# fixers

# MLE
pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_f<-pv

NLL2_rec_d_t_f <- function(af,bf,cf,df,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	yf <- fn1_rec_d_t(af, bf, cf, df, D[F==1], A[F==1], Dhat, t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit2_rec_d_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_f)

fit2_rec_d_c_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_f)

# non-fixers

pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_n<-pv

NLL2_rec_d_t_n <- function(a0,b0,c0,d0,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	y0 <- fn1_rec_d_t(a0, b0, c0, d0, D[F==0], A[F==0], Dhat, t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit2_rec_d_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_n)

fit2_rec_d_c_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_n)

################ Ricker

fn4_rec_d_t<-function(a,c,d,k,D,A,Dhat,t){
		a<-exp(a)
	c<-exp(c)
y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	y<-y*t
	} #Right-skewed with diameter effect

# fixers

# MLE

pv<-c(-4,-8.6,-3.7,-1.6,-0.6,-0.2,0.2,0.18)
pv_f<-pv

NLL5_rec_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	yf<-fn4_rec_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit5_rec_d_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_f)

fit5_rec_d_c_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_f)

# non-fixers

# MLE

NLL5_rec_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	y0<-fn4_rec_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit5_rec_d_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_n)

fit5_rec_d_c_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_n)

######################################################
################   Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

# fixer
deltaAICs(c(AIC(fit2_rec_d_t_f),AIC(fit2_rec_d_c_t_f),AIC(fit5_rec_d_t_f),AIC(fit5_rec_d_c_t_f)))

# [1] 9.81149
# [1] 0.2780213
# [1] 10.4186
# [1] 0

# non-fixers
deltaAICs(c(AIC(fit2_rec_d_t_n),AIC(fit2_rec_d_c_t_n),AIC(fit5_rec_d_t_n),AIC(fit5_rec_d_c_t_n)))
# [1] 0.9851506
# [1] 0
# [1] 329.4753
# [1] 321.6797

#setwd("")
pv_f_rr_robinia<-coef(fit5_rec_d_c_t_f)
pv_n_rr_robinia<-coef(fit2_rec_d_c_t_n)
write.csv(pv_f_rr_robinia,file="pv_f_rr_robinia.csv")
write.csv(pv_n_rr_robinia,file="pv_n_rr_robinia.csv")

######################################################
#################    CI Computation   ################
######################################################

agevec_rr_robinia<-unique(MRdat_robinia$stdage)
mD<-mean(MRdat_robinia$mDBH)
mD_fix<-mean(MRdat_robinia[MRdat_robinia$FIX==1,]$mDBH)
mD_non<-mean(MRdat_robinia[MRdat_robinia$FIX==0,]$mDBH)

ndraws<-10000
xs_rr_robinia<-seq(min(agevec_rr_robinia),max(agevec_rr_robinia),1)
xs<-seq(min(agevec_rr_robinia),max(agevec_rr_robinia),1)

vmat_fix_rr<-mvrnorm(ndraws,mu=pv_f_rr_robinia,Sigma=vcov(fit5_rec_d_c_t_f)) 
vmat_non_rr<-mvrnorm(ndraws,mu=pv_n_rr_robinia,Sigma=vcov(fit2_rec_d_c_t_n))

ys_fix_rr<-array(dim=c(ndraws,length(xs)))
ys_non_rr<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix_rr[i,]<-fn4_rec_d_t(vmat_fix_rr[i,1],vmat_fix_rr[i,2],vmat_fix_rr[i,3],vmat_fix_rr[i,4],mD_fix,xs,mD,1)*100
	ys_non_rr[i,]<-fn1_rec_d_t(vmat_non_rr[i,1],vmat_non_rr[i,2],vmat_non_rr[i,3],vmat_non_rr[i,4],mD_non,xs,mD,1)*100
}

civec_fix_rr<-array(dim=c(2,length(xs)))
civec_non_rr<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix_rr[,j]<-quantile(ys_fix_rr[,j],c(lowb,upb),na.rm=TRUE)
	civec_non_rr[,j]<-quantile(ys_non_rr[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_rr_robinia<-civec_fix_rr
civec_non_rr_robinia<-civec_non_rr

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

#setwd("")
pdf("Fig10_RRfit_wCI_Robinia_150928.pdf",width=4.5,height=4.5)

x<-c(0,180)
y<-c(0,100)
plot(y~x,type="n",xlim=c(0,180),y=c(0,100),xlab="Stand Age (Years)",ylab="Recruitment Rate (% Recruited/year)",main=expression(italic(Robinia~pseudoacacia)))
points(navg_rr_robinia$meanRR*100~navg_rr_robinia$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",col="blue",main="U.S. Trees",pch=2,cex=0.8)
points(fxavg_rr_robinia$meanRR*100~fxavg_rr_robinia$stdage,col="red",pch=1,cex=0.8)

# fixers
	# Ricker
curve(fn4_rec_d_t(pv_f_rr_robinia[1],pv_f_rr_robinia[2],pv_f_rr_robinia[3],pv_f_rr_robinia[4],mD_fix,x,mD,1)*100,from=min(agevec_rr_robinia),to=max(agevec_rr_robinia),lwd=3,lty=1,add=TRUE,col="red")

# non-fixers
	# Saturating
curve(fn1_rec_d_t(pv_n_rr_robinia[1],pv_n_rr_robinia[2],pv_n_rr_robinia[3],pv_n_rr_robinia[4],mD_non,x,mD,1)*100,from=min(agevec_rr_robinia),to=max(agevec_rr_robinia),lwd=3,col="blue",lty=1,add=TRUE,xlab="Stand Age",ylab="Estimated Relative Recruitment Rate (%)")

legend("topright",c(expression(italic(Robinia~pseudoacacia)),"Non-fixers"),lty=1,bty="n",cex=1.2,lwd=2,col=c("red","blue"),pch=c(1,2))

points(civec_fix_rr_robinia[1,]~xs_rr_robinia,lwd=2,typ="l",lty=2,col="red")
points(civec_fix_rr_robinia[2,]~xs_rr_robinia,lwd=2,typ="l",lty=2,col="red")

points(civec_non_rr_robinia[1,]~xs_rr_robinia,lwd=2,typ="l",lty=2,col="blue") # This CI is really wierd. Need to troubleshoot.
points(civec_non_rr_robinia[2,]~xs_rr_robinia,lwd=2,typ="l",lty=2,col="blue")

dev.off()





#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################
#################################






####################################################################################################
#############################         Crown Class Growth Analysis         ##########################
####################################################################################################

# setwd("")
t<-read.csv("tree_small.csv",header=T)
t_original<-t
tcn_cclcd<-t_original[,c("tcn","cclcd")]

# setwd("")
write.csv(tcn_cclcd,file="FIA_tcn_cclcd_150715.csv")
rm(t)

# setwd("")
gdata<-read.csv("FIAGRSTAT_final_data_for_fit_150707.csv",header=T)

gdata<-gdata[,c(-1)]
gdata_original<-gdata

gdata_cclcd<-merge(gdata,tcn_cclcd,by.x="tcn",by.y="tcn",all.x=TRUE,all.y=FALSE)
# setwd("")
write.csv(gdata_cclcd,file="FIA_GR_CCLCD_data_for_fit_150715.csv")

######################################################
##############    Mean Calculations   ################
######################################################

# cclcd==1,2,3 canopy
gdata<-gdata_cclcd[gdata_cclcd$cclcd==1 | gdata_cclcd$cclcd==2 | gdata_cclcd$cclcd==3,]

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

u_stdage_canopy<-u_stdage
u_meanLGR_fix_canopy<-u_meanLGR_fix
u_meanLGR_non_canopy<-u_meanLGR_non

# cclcd==4,5 understory
gdata<-gdata_cclcd[gdata_cclcd$cclcd==4 | gdata_cclcd$cclcd==5,]

u_stdage <- sort(unique(gdata$stdage))
u_meanLGR_fix <- rep(NA,length(u_stdage))
u_meanLGR_non <- rep(NA,length(u_stdage))
for(j in 1:length(u_stdage)){
	print(j)
	u_meanLGR_fix[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==1,]$LGRp)))
	u_meanLGR_non[j] <- exp(mean(log(gdata[gdata$stdage==u_stdage[j] & gdata$FIX==0,]$LGRp)))
}

u_stdage_understory<-u_stdage
u_meanLGR_fix_understory<-u_meanLGR_fix
u_meanLGR_non_understory<-u_meanLGR_non

############################################################################################################
################################ Other Other Other Other  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# setwd("")
pdf("Fig7_GRmean_CrownClass_150715.pdf",width=4.5,height=4.5)

plot(u_meanLGR_fix_canopy*100~u_stdage_canopy,xlab="Stand Age",ylab="Growth Rate (% diamter growth/year)",main="Canopy and Understory Trees",pch=19,cex=0.8,col="red",xlim=c(0,250),ylim=c(0,6),cex.lab=1.5,cex.main=2)
points(u_meanLGR_non_canopy*100~u_stdage_canopy,pch=17,cex=0.8,col="blue")
points(u_meanLGR_fix_understory*100~u_stdage_understory,pch=19,cex=0.8,col="green")
points(u_meanLGR_non_understory*100~u_stdage_understory,pch=17,cex=0.8,col="yellow")
legend("topright",legend=c("Canopy Fixers","Canopy Non-fixers","Understory Fixers","Understory Non-fixers"),col=c("red","blue","green","yellow"),pch=c(19,17,19,17),cex=1,bty="n")

dev.off()

######################################################
##################    STAT Analysis   ################
######################################################

MDL<-0.05
lMDL<-function(D,T){
	log((log(D+MDL)-log(D))/T)
}

### Canopy

gdata<-gdata_cclcd[gdata_cclcd$cclcd==1 | gdata_cclcd$cclcd==2 | gdata_cclcd$cclcd==3,]
gdata_canopy<-gdata_cclcd[gdata_cclcd$cclcd==1 | gdata_cclcd$cclcd==2 | gdata_cclcd$cclcd==3,]

################ Saturating

# fixers

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)

	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1.2,-10,-4,30,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))

u_stdage<-u_stdage_canopy
u_meanLGR_fix<-u_meanLGR_fix_canopy

plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

# non-fixers

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5.6,-3.8,50,-1)

	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5.6,-3.8,50,-1)

u_stdage<-u_stdage_canopy
u_meanLGR_non<-u_meanLGR_non_canopy

plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################ Ricker

# fixers

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

# non-fixers

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)

	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
###############    Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_4_f)))
# [1] 0
# [1] 36.01922

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_4_n)))
# [1] 0
# [1] 211.9182

pv_f_gr_canopy<-coef(fit_logGR_1_f)
pv_n_gr_canopy<-coef(fit_logGR_1_n)

######################################################
###############     CI Computation    ################
######################################################

library(MASS)

agevec_gr_canopy<-unique(gdata_canopy$stdage)

ndraws<-1000
xs<-seq(min(agevec_gr_canopy),max(agevec_gr_canopy), 1)
xs_gr_canopy<-seq(min(agevec_gr_canopy),max(agevec_gr_canopy), 1)

vmat_fix<-mvrnorm(ndraws, mu=pv_f_gr_canopy[2:5],Sigma=vcov(fit_logGR_1_f)[c(2:5),c(2:5)])
vmat_non<-mvrnorm(ndraws, mu=pv_n_gr_canopy[2:5],Sigma=vcov(fit_logGR_1_n)[c(2:5),c(2:5)])

ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))


for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_1_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_1_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_gr_canopy<-civec_fix
civec_non_gr_canopy<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################
# setwd("")
pdf("Fig11_GRfit_wCI_Canopy_Updated_150928.pdf",width=4.5,height=4.5)

plot(u_stdage_canopy,u_meanLGR_fix_canopy*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (% diameter growth/year)", main="Canopy",cex=1,pch=1,ylim=c(0,4))
points(u_stdage_canopy,u_meanLGR_non_canopy*100,col="blue",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f_gr_canopy[2],pv_f_gr_canopy[3],pv_f_gr_canopy[4],pv_f_gr_canopy[5],exp(mean(log(gdata$dbh)))
,x))*100,from=min(agevec_gr_canopy),to=max(agevec_gr_canopy),lwd=4,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
curve(exp(gMLF_1_n(pv_n_gr_canopy[2],pv_n_gr_canopy[3],pv_n_gr_canopy[4],pv_n_gr_canopy[5],exp(mean(log(gdata$dbh)))
,x))*100,from=min(agevec_gr_canopy),to=max(agevec_gr_canopy),lwd=4,lty=1,col="blue",add=TRUE)

points(civec_fix_gr_canopy[1,]~xs_gr_canopy,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_gr_canopy[2,]~xs_gr_canopy,lwd=2,lty=2,typ="l",col="red")

points(civec_non_gr_canopy[1,]~xs_gr_canopy,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_gr_canopy[2,]~xs_gr_canopy,lwd=2,lty=2,typ="l",col="blue")

legend("topright", bty="n",legend=c("Canopy N fixers","Canopy Non-fixers "),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2),cex=1.2)

dev.off()

### Understory

gdata<-gdata_cclcd[gdata_cclcd$cclcd==4 | gdata_cclcd$cclcd==5,]
gdata_understory<-gdata_cclcd[gdata_cclcd$cclcd==4 | gdata_cclcd$cclcd==5,]

################ Saturating

# fixers

	# Functions
gMLF_1_f<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_f<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_f(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_f(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1.2,-10,-4,69,-0.2)

	# MLE
fit_logGR_1_f<-mle2(normNNLF_1_f,start=list(sd=1,a=pv_f[2],b=pv_f[3],c=pv_f[4],d=pv_f[5]), 
	data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_1_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1.2,-10,-4,69,-0.2)
agevec_noOG<-sort(unique(gdata$stdage))

u_stdage<-u_stdage_understory
u_meanLGR_fix<-u_meanLGR_fix_understory

plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Saturating Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_1_f)
curve(exp(gMLF_1_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

# non-fixers

	# Functions
gMLF_1_n<-function(a,b,c,d,D,A){
	a<-exp(a)
	b<-exp(b)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_1_n<-function(sd,a,b,c,d,G,D,A,T){
	up<-gMLF_1_n(a,b,c,d,D[G>0],A[G>0])#positive
	un<-gMLF_1_n(a,b,c,d,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5.6,-4.5,45,-0.59)

	# MLE
fit_logGR_1_n<-mle2(normNNLF_1_n,start=list(sd=1,a=pv_n[2],b=pv_n[3],c=pv_n[4],d=pv_n[5]), 
	data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_1_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5.6,-4.5,45,-0.59)

u_stdage<-u_stdage_understory
u_meanLGR_non<-u_meanLGR_non_understory

plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Saturating Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_1_n)
curve(exp(gMLF_1_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

################ Ricker

# fixers

	# Functions
gMLF_4_f<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_f<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_f(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_f(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)

	# MLE
fit_logGR_4_f<-mle2(normNNLF_4_f,start=list(sd=1,a=pv_f[2],c=pv_f[3],d=pv_f[4],k=pv_f[5]), data=list(G=gdata[gdata$FIX==1,]$LGR,D=gdata[gdata$FIX==1,]$dbh,A=gdata[gdata$FIX==1,]$stdage,T=gdata[gdata$FIX==1,]$t))
summary(fit_logGR_4_f)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_f<-c(1,-6.6,-2.3,-0.2,0.04)
agevec_noOG<-sort(unique(gdata$stdage))
plot(u_stdage,u_meanLGR_fix*100,col="red",xlab="Stand Age",ylab="Growth Rate", main="N Fixer Ricker Curve Fit",cex=1,pch=1,ylim=c(0,4))
curve(exp(gMLF_4_f(pv_f[2],pv_f[3],pv_f[4],pv_f[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")

pv_f_fit<-coef(fit_logGR_4_f)
curve(exp(gMLF_4_f(pv_f_fit[2],pv_f_fit[3],pv_f_fit[4],pv_f_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=4,lty=2,col="green",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("red","green"),lwd=2)

# non-fixers

	# Functions
gMLF_4_n<-function(a,c,d,k,D,A){
	a<-exp(a)
	c<-exp(c)
	Dhat<-exp(mean(log(D))) 
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

	# Negative loglikelihood
normNNLF_4_n<-function(sd,a,c,d,k,G,D,A,T){
	up<-gMLF_4_n(a,c,d,k,D[G>0],A[G>0])#positive
	un<-gMLF_4_n(a,c,d,k,D[G<=0],A[G<=0]) #negative
	m<-lMDL(D[G<=0],T[G<=0])
	-(sum(dnorm(log(G[G>0]),mean=up,sd=sd,log=TRUE))+sum(pnorm(m,mean=un,sd=sd,log=TRUE)))
} 

	# Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)

	# MLE
fit_logGR_4_n<-mle2(normNNLF_4_n,start=list(sd=1,a=pv_n[2],c=pv_n[3],d=pv_n[4],k=pv_n[5]), data=list(G=gdata[gdata$FIX==0,]$LGR,D=gdata[gdata$FIX==0,]$dbh,A=gdata[gdata$FIX==0,]$stdage,T=gdata[gdata$FIX==0,]$t))
summary(fit_logGR_4_n)

	# Plot Fits with Initial Values and MLE Parameter Values
pv_n<-c(1,-5,-2.3,-0.6,0.05)
plot(u_stdage,u_meanLGR_non*100,col="blue",xlab="Stand Age",ylab="Growth Rate", main="Non-fixer Ricker Curve Fit",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_4_n(pv_n[2],pv_n[3],pv_n[4],pv_n[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=1,col="blue",add=TRUE)

pv_n_fit<-coef(fit_logGR_4_n)
curve(exp(gMLF_4_n(pv_n_fit[2],pv_n_fit[3],pv_n_fit[4],pv_n_fit[5],exp(mean(log(gdata$dbh))),x))*100,from=min(agevec_noOG),to=max(agevec_noOG),lwd=6,lty=2,col="green",add=TRUE)
legend("topright",legend=c("Prior","Fit"),lty=c(1,2),col=c("blue","green"),lwd=2)

######################################################
###############    Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

	# For N Fixers
deltaAICs(c(AIC(fit_logGR_1_f),AIC(fit_logGR_4_f)))
# [1] 0
# [1] 0.1300251

	# For Non Fixers
deltaAICs(c(AIC(fit_logGR_1_n),AIC(fit_logGR_4_n)))
# [1] 0
# [1] 24.42865

pv_f_gr_understory<-coef(fit_logGR_1_f)
pv_n_gr_understory<-coef(fit_logGR_1_n)

######################################################
###############     CI Computation    ################
######################################################

library(MASS)

agevec_gr_understory<-unique(gdata_understory$stdage)

ndraws<-1000
xs<-seq(min(agevec_gr_understory),max(agevec_gr_understory), 1)
xs_gr_understory<-seq(min(agevec_gr_understory),max(agevec_gr_understory), 1)

vmat_fix<-mvrnorm(ndraws, mu=pv_f_gr_understory[2:5],Sigma=vcov(fit_logGR_1_f)[c(2:5),c(2:5)])
vmat_non<-mvrnorm(ndraws, mu=pv_n_gr_understory[2:5],Sigma=vcov(fit_logGR_1_n)[c(2:5),c(2:5)])

ys_fix<-array(dim=c(ndraws, length(xs)))
ys_non<-array(dim=c(ndraws, length(xs)))


for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-exp(gMLF_1_f(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(gdata$dbh))),xs))*100
	ys_non[i,]<-exp(gMLF_1_n(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(gdata$dbh))),xs))*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_gr_understory<-civec_fix
civec_non_gr_understory<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# setwd("")

pdf("Fig12_GRfit_wCI_Understory_Updated_150928.pdf",width=4.5,height=4.5)

plot(u_stdage_understory,u_meanLGR_fix_understory*100,col="red",xlab="Stand Age (Years)",ylab="Growth Rate (% diameter growth/year)", main="Understory",cex=1,pch=1,ylim=c(0,4))
points(u_stdage_understory,u_meanLGR_non_understory*100,col="blue",cex=1,pch=2,ylim=c(0,4))
curve(exp(gMLF_1_f(pv_f_gr_understory[2],pv_f_gr_understory[3],pv_f_gr_understory[4],pv_f_gr_understory[5],exp(mean(log(gdata$dbh)))
,x))*100,from=min(agevec_gr_understory),to=max(agevec_gr_understory),lwd=3,lty=1,col="red",add=TRUE,xlab="Stand Age",ylab="Estimated Relative Growth Rate",main="Average-size Tree")
curve(exp(gMLF_1_n(pv_n_gr_understory[2],pv_n_gr_understory[3],pv_n_gr_understory[4],pv_n_gr_understory[5],exp(mean(log(gdata$dbh)))
,x))*100,from=min(agevec_gr_understory),to=max(agevec_gr_understory),lwd=3,lty=1,col="blue",add=TRUE)

points(civec_fix_gr_understory[1,]~xs_gr_understory,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_gr_understory[2,]~xs_gr_understory,lwd=2,lty=2,typ="l",col="red")

points(civec_non_gr_understory[1,]~xs_gr_understory,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_gr_understory[2,]~xs_gr_understory,lwd=2,lty=2,typ="l",col="blue")

legend("topright", bty="n",legend=c("Understory N-fixers","Understory Non-fixers "),col=c("red","blue"),lty=1,lwd=2,pch=c(1,2),cex=1.2)

dev.off()

####################################################################################################
##########################         Crown Class Mortality Analysis         ##########################
####################################################################################################

# setwd("")

mdata<-read.csv("FIAMRSTAT_morality_data_stat_final_150709.csv",header=T)
tcn_cclcd<-read.csv("FIA_tcn_cclcd_150715.csv",header=T)
mdata<-mdata[,c(-1)]
tcn_cclcd<-tcn_cclcd[,c(-1)]

mdata_cclcd<-merge(mdata,tcn_cclcd,by.x="tcn0",by.y="tcn",all.x=FALSE,all.y=FALSE)

######################################################
#############     Mean Calculations    ###############
######################################################

### Canopy

mdata_canopy<-mdata_cclcd[mdata_cclcd$cclcd==1 | mdata_cclcd$cclcd==2 | mdata_cclcd$cclcd==3,]

##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1

#Now: Fixers

mdata<-mdata_canopy
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_canopy_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_canopy
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_canopy_v1<-favg
navg_mr<-favg

##################  Mean Calculation  ################ Version 2: plot with S=0, omit

#Now: Fixers

mdata<-mdata_canopy
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_canopy_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_canopy
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_canopy_v2<-favg
navg_mr<-favg

######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)
mdata_plotInfo

# setwd("")

write.csv(mdata_plotInfo,file="FIA_mdata_canopy_plotInfo_S0_150728.csv")

### Understory

mdata_understory<-mdata_cclcd[mdata_cclcd$cclcd==4 | mdata_cclcd$cclcd==5,]

##################  Mean Calculation  ################ Version 1: plot with S=0, assigned mr=1

#Now: Fixers

mdata<-mdata_understory
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_fix_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_understory_v1<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_understory
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS[pAgeNS$S==0,]$MR<-1
# hist(pAgeNS$MR)
# plot(pAgeNS$MR~pAgeNS$stdage,xlim=c(0,280))

pAgeNS_non_wS0<-pAgeNS

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_understory_v1<-favg
navg_mr<-favg

##################  Mean Calculation  ################ Version 2: plot with S=0, omit

#Now: Fixers

mdata<-mdata_understory
mdata<-mdata[mdata$FIX==1,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS[pAgeNS$S!=0,]$MR<-log(pAgeNS[pAgeNS$S!=0,]$N/pAgeNS[pAgeNS$S!=0,]$S)/pAgeNS[pAgeNS$S!=0,]$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

fxavg_mr_understory_v2<-favg
fxavg_mr<-favg

# Now: For non fixers
mdata<-mdata_understory
mdata<-mdata[mdata$FIX==0,]

plotN<-tapply(mdata$N,mdata$pcn,sum)
plotS<-tapply(mdata$S,mdata$pcn,sum)

pN<-t(plotN)
pN<-t(pN)
plN<-data.frame(pcn=as.numeric(rownames(pN)),N=pN)
pS<-t(plotS)
pS<-t(pS)
plS<-data.frame(pcn=as.numeric(rownames(pS)),S=pS)

pNS<-merge(plN,plS,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
age_mt<-data.frame(pcn=mdata$pcn,stdage=mdata$stdage,mt=mdata$mt)
age_mt<-subset(age_mt,!duplicated(age_mt$pcn))
pAgeNS<-merge(age_mt,pNS,by.x="pcn",by.y="pcn",all.x=FALSE,all.y=FALSE)

pAgeNS.original<-pAgeNS
pAgeNS$MR<-NA
pAgeNS$MR<-log(pAgeNS$N/pAgeNS$S)/pAgeNS$mt
pAgeNS<-pAgeNS[pAgeNS$S!=0,]

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

avg<-tapply(pAgeNS$MR,pAgeNS$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanMR=tavg)

navg_mr_understory_v2<-favg
navg_mr<-favg

######################################################
###############    Check Plots with S=0   ############
######################################################

mdata_plotInfo<-data.frame(FIX=c(1,0),NoPlotS0=NA,NoPlotN0=NA,percPlotS0=NA,NoTreeinPlotS0=NA,percTreeinPlotS0=NA)
mdata_plotInfo$NoPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])
mdata_plotInfo$NoPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])
mdata_plotInfo$NoPlotN0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$NoPlotN0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_fix_wS0$N==0,])
mdata_plotInfo$percPlotS0[1]<-nrow(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,])/nrow(pAgeNS_fix_wS0)
mdata_plotInfo$percPlotS0[2]<-nrow(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,])/nrow(pAgeNS_non_wS0)
mdata_plotInfo$NoTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)
mdata_plotInfo$NoTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)
mdata_plotInfo$percTreeinPlotS0[1]<-sum(pAgeNS_fix_wS0[pAgeNS_fix_wS0$S==0,]$N)/sum(pAgeNS_fix_wS0$N)
mdata_plotInfo$percTreeinPlotS0[2]<-sum(pAgeNS_non_wS0[pAgeNS_non_wS0$S==0,]$N)/sum(pAgeNS_non_wS0$N)
mdata_plotInfo

# setwd("")

write.csv(mdata_plotInfo,file="FIA_mdata_understory_plotInfo_S0_150728.csv")

######################################################
###############     STAT Analysis     ################
######################################################

library(bbmle)

### Canopy

mdata<-mdata_canopy

############### Saturating

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	}

# fixers

# MLE Fit
x<-0:300
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

fxavg_mr<-fxavg_mr_canopy_v2

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

navg_mr<-navg_mr_canopy_v2

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

############### Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
###############   Model Comparison    ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
# [1] 0.6372045
# [1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
# [1] 13.30091
# [1] 0

pv_f_mr_canopy<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_canopy<-coef(fit_ricker_surv_mi_d_t_n)

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

agevec_mr_canopy<-unique(mdata_canopy$stdage)

ndraws<-1000
xs_mr_canopy<-seq(min(agevec_mr_canopy),max(agevec_mr_canopy),1)
xs<-seq(min(agevec_mr_canopy),max(agevec_mr_canopy),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_canopy,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_canopy,Sigma=vcov(fit_ricker_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata_canopy$dbh0))),xs,exp(mean(log(mdata_canopy$dbh0))),1)*100
	ys_non[i,]<-fn_ricker_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata_canopy$dbh0))),xs,exp(mean(log(mdata_canopy$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_canopy<-civec_fix
civec_non_mr_canopy<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# setwd("")

pdf("Fig13_MRfit_wCI_Canopy_version1_Updated_150731.pdf",width=4.5,height=4.5)

plot(fxavg_mr_canopy_v1$meanMR*100~fxavg_mr_canopy_v1$stdage,xlim=c(min(agevec_mr_canopy),max(agevec_mr_canopy)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Canopy",ylim=c(0,100))
points(navg_mr_canopy_v1$meanMR*100~navg_mr_canopy_v1$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_canopy[1],pv_f_mr_canopy[2],pv_f_mr_canopy[3],pv_f_mr_canopy[4],exp(mean(log(mdata_canopy$dbh0))),x,exp(mean(log(mdata_canopy$dbh0))),1)*100,from=min(agevec_mr_canopy),to=max(agevec_mr_canopy),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_canopy[1],pv_n_mr_canopy[2],pv_n_mr_canopy[3],pv_n_mr_canopy[4],exp(mean(log(mdata_canopy$dbh0))),x,exp(mean(log(mdata_canopy$dbh0))),1)*100,from=min(agevec_mr_canopy),to=max(agevec_mr_canopy),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_canopy[1,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_canopy[2,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_canopy[1,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_canopy[2,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("Canopy Fixers","Canopy Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

# setwd("")

pdf("Fig14_MRfit_wCI_Canopy_version2_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_canopy_v2$meanMR*100~fxavg_mr_canopy_v2$stdage,xlim=c(min(agevec_mr_canopy),max(agevec_mr_canopy)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Canopy",ylim=c(0,15))
points(navg_mr_canopy_v2$meanMR*100~navg_mr_canopy_v2$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_canopy[1],pv_f_mr_canopy[2],pv_f_mr_canopy[3],pv_f_mr_canopy[4],exp(mean(log(mdata_canopy$dbh0))),x,exp(mean(log(mdata_canopy$dbh0))),1)*100,from=min(agevec_mr_canopy),to=max(agevec_mr_canopy),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_ricker_mort_d_t(pv_n_mr_canopy[1],pv_n_mr_canopy[2],pv_n_mr_canopy[3],pv_n_mr_canopy[4],exp(mean(log(mdata_canopy$dbh0))),x,exp(mean(log(mdata_canopy$dbh0))),1)*100,from=min(agevec_mr_canopy),to=max(agevec_mr_canopy),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_canopy[1,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_canopy[2,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_canopy[1,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_canopy[2,]~xs_mr_canopy,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("Canopy N fixers","Canopy Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

######################################################
###############     STAT Analysis     ################
######################################################

library(bbmle)

### Understory

mdata<-mdata_understory

############### Saturating

fn_sat_surv_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m*t)
	} # saturating with diameter effect

fn_sat_mort_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	}

# fixers

# MLE Fit
x<-0:300
pv<-c(-4.5,-4.2,  -3.2,-4,  10,12,  -0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_sat_surv_mort_d_t_f <- function(a,b,c,d,N,S,D,A,F,t){
	uf <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==1], A[F==1], exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_sat_surv_mi_d_t_f <- mle2(NLL_sat_surv_mort_d_t_f,start=list(a=pv_f[1],b=pv_f[2],c=pv_f[3],d=pv_f[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-4.5,-4.2,  -3.2,-4,  10,12,  -0.1,-0.1) 
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_sat_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_sat_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_sat_surv_mi_d_t_f)
mort_f_new<-fn_sat_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

fxavg_mr<-fxavg_mr_understory_v2

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_sat_surv_mort_d_t_n <- function(a,b,c,d,N,S,D,A,F,t){
	un <- fn_sat_surv_mort_d_t(a,b,c,d, D[F==0], A[F==0], exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_sat_surv_mi_d_t_n <- mle2(NLL_sat_surv_mort_d_t_n,start=list(a=pv_n[1],b=pv_n[2],c=pv_n[3],d=pv_n[4]),
	data=list(N=mdata$N,S=mdata$S,D=mdata$dbh0,A=mdata$stdage,F=mdata$FIX,t=mdata$mt))
summary(fit_sat_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-4.5,-4.2,-3.2,-2.9,10,10,-0.1,-0.1) 
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_sat_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_sat_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_sat_surv_mi_d_t_n)
mort_n_new<-fn_sat_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

navg_mr<-navg_mr_understory_v2

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Saturating",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

############### Ricker

fn_ricker_surv_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m*t) # This is the survival probability
}

fn_ricker_mort_d_t<-function(a,c,d,k,D,A,Dhat,t){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

# fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-5.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])

NLL_ricker_surv_mort_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	uf<-fn_ricker_surv_mort_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
	-(sum(dbinom(S[F==1],size=N[F==1],prob=uf,log=TRUE)))
}

fit_ricker_surv_mi_d_t_f<-mle2(NLL_ricker_surv_mort_d_t_f,start=list(af=pv_f[1],cf=pv_f[2],df=pv_f[3],kf=pv_f[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_f)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_f<-c(pv[2],pv[4],pv[6],pv[8])
gamma_f<-fn_ricker_surv_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_f<-fn_ricker_mort_d_t(pv_f[1],pv_f[2],pv_f[3],pv_f[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_f_new<-coef(fit_ricker_surv_mi_d_t_f)
mort_f_new<-fn_ricker_mort_d_t(pv_f_new[1],pv_f_new[2],pv_f_new[3],pv_f_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(fxavg_mr$meanMR~fxavg_mr$stdage,xlim=c(0,300),col="red",cex=0.8,pch=1,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate N Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_f~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_f_new~x,col="red",lwd=3,typ="l")
legend("topright",legend=c("Initial","Fit"),col=c("dark green","red"),lty=c(2,1),lwd=2,bty="n")

# non fixers

# MLE Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])

NLL_ricker_surv_mort_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	un<-fn_ricker_surv_mort_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-(sum(dbinom(S[F==0],size=N[F==0],prob=un,log=TRUE)))
}

fit_ricker_surv_mi_d_t_n<-mle2(NLL_ricker_surv_mort_d_t_n,start=list(a0=pv_n[1],c0=pv_n[2],d0=pv_n[3],k0=pv_n[4]),data=list(N=mdata$N,S=mdata$S,A=mdata$stdage,F=mdata$FIX,D=mdata$dbh0,t=mdata$mt))
summary(fit_ricker_surv_mi_d_t_n)

# Plot out the Fit
pv<-c(-5.6,-3.6,-2.9,-3.0,-0.1,-0.1,0.08,0.05)
pv_n<-c(pv[1],pv[3],pv[5],pv[7])
gamma_n<-fn_ricker_surv_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)
mort_n<-fn_ricker_mort_d_t(pv_n[1],pv_n[2],pv_n[3],pv_n[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

pv_n_new<-coef(fit_ricker_surv_mi_d_t_n)
mort_n_new<-fn_ricker_mort_d_t(pv_n_new[1],pv_n_new[2],pv_n_new[3],pv_n_new[4],exp(mean(log(mdata$dbh0))),x,exp(mean(log(mdata$dbh0))),1)

plot(navg_mr$meanMR~navg_mr$stdage,xlim=c(0,300),col="blue",cex=0.8,pch=2,xlab="Stand Age",ylab="Relative Mortality Rate",main="Mortality Rate Non Fixers - Ricker",ylim=c(0,0.2))
agevec_noOG<-unique(mdata$stdage)
points(mort_n~x,col="dark green",typ="l",lwd=3,lty=2)
points(mort_n_new~x,col="blue",typ="l",lwd=3)
legend("topright",legend=c("Initial","Fit"),col=c("dark green","blue"),lty=c(2,1),lwd=3,bty="n")

######################################################
###############   Model Comparison    ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_f),AIC(fit_ricker_surv_mi_d_t_f)))
# [1] 27.322
# [1] 0

deltaAICs(c(AIC(fit_sat_surv_mi_d_t_n),AIC(fit_ricker_surv_mi_d_t_n)))
# [1] 0
# [1] 3.346047

pv_f_mr_understory<-coef(fit_ricker_surv_mi_d_t_f)
pv_n_mr_understory<-coef(fit_sat_surv_mi_d_t_n)

######################################################
############    CI Computation   #####################
######################################################

library(MASS)

agevec_mr_understory<-unique(mdata_understory$stdage)

ndraws<-1000
xs_mr_understory<-seq(min(agevec_mr_understory),max(agevec_mr_understory),1)
xs<-seq(min(agevec_mr_understory),max(agevec_mr_understory),1)
vmat_fix<-mvrnorm(ndraws,mu=pv_f_mr_understory,Sigma=vcov(fit_ricker_surv_mi_d_t_f))
vmat_non<-mvrnorm(ndraws,mu=pv_n_mr_understory,Sigma=vcov(fit_sat_surv_mi_d_t_n))

ys_fix<-array(dim=c(ndraws,length(xs)))
ys_non<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix[i,]<-fn_ricker_mort_d_t(vmat_fix[i,1],vmat_fix[i,2],vmat_fix[i,3],vmat_fix[i,4],exp(mean(log(mdata_understory$dbh0))),xs,exp(mean(log(mdata_understory$dbh0))),1)*100
	ys_non[i,]<-fn_sat_mort_d_t(vmat_non[i,1],vmat_non[i,2],vmat_non[i,3],vmat_non[i,4],exp(mean(log(mdata_understory$dbh0))),xs,exp(mean(log(mdata_understory$dbh0))),1)*100
}

civec_fix<-array(dim=c(2,length(xs)))
civec_non<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix[,j]<-quantile(ys_fix[,j],c(lowb,upb),na.rm=TRUE)
	civec_non[,j]<-quantile(ys_non[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_mr_understory<-civec_fix
civec_non_mr_understory<-civec_non

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################
# setwd("")

pdf("Fig13_MRfit_wCI_Understory_version1_Updated_150731.pdf",width=4.5,height=4.5)

plot(fxavg_mr_understory_v1$meanMR*100~fxavg_mr_understory_v1$stdage,xlim=c(min(agevec_mr_understory),max(agevec_mr_understory)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Understory",ylim=c(0,100))
points(navg_mr_understory_v1$meanMR*100~navg_mr_understory_v1$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_understory[1],pv_f_mr_understory[2],pv_f_mr_understory[3],pv_f_mr_understory[4],exp(mean(log(mdata_understory$dbh0))),x,exp(mean(log(mdata_understory$dbh0))),1)*100,from=min(agevec_mr_understory),to=max(agevec_mr_understory),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_sat_mort_d_t(pv_n_mr_understory[1],pv_n_mr_understory[2],pv_n_mr_understory[3],pv_n_mr_understory[4],exp(mean(log(mdata_understory$dbh0))),x,exp(mean(log(mdata_understory$dbh0))),1)*100,from=min(agevec_mr_understory),to=max(agevec_mr_understory),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_understory[1,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_understory[2,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_understory[1,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_understory[2,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("Understory Fixers","Understory Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

# setwd("")

pdf("Fig14_MRfit_wCI_Understory_version2_Updated_150928.pdf",width=4.5,height=4.5)

plot(fxavg_mr_understory_v2$meanMR*100~fxavg_mr_understory_v2$stdage,xlim=c(min(agevec_mr_understory),max(agevec_mr_understory)),col="red",cex=1,pch=1,xlab="Stand Age (Years)",ylab="Mortality Rate (%dead/year)",main="Understory",ylim=c(0,25))
points(navg_mr_understory_v2$meanMR*100~navg_mr_understory_v2$stdage,xlim=c(0,2250),col="blue",cex=1,pch=2)

curve(fn_ricker_mort_d_t(pv_f_mr_understory[1],pv_f_mr_understory[2],pv_f_mr_understory[3],pv_f_mr_understory[4],exp(mean(log(mdata_understory$dbh0))),x,exp(mean(log(mdata_understory$dbh0))),1)*100,from=min(agevec_mr_understory),to=max(agevec_mr_understory),lwd=3,col="red",lty=1,add=TRUE)
curve(fn_sat_mort_d_t(pv_n_mr_understory[1],pv_n_mr_understory[2],pv_n_mr_understory[3],pv_n_mr_understory[4],exp(mean(log(mdata_understory$dbh0))),x,exp(mean(log(mdata_understory$dbh0))),1)*100,from=min(agevec_mr_understory),to=max(agevec_mr_understory),lwd=3,col="blue",lty=1,add=TRUE)

points(civec_fix_mr_understory[1,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="red")
points(civec_fix_mr_understory[2,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="red")
points(civec_non_mr_understory[1,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="blue")
points(civec_non_mr_understory[2,]~xs_mr_understory,lwd=2,lty=2,typ="l",col="blue")

legend('topright',c("Understory N fixers","Understory Non-fixers"),col=c("red","blue"),lty=1,bty="n",cex=1,lwd=2,pch=c(1,2))

dev.off()

####################################################################################################
########################         Crown Class Recruitment Analysis         ##########################
####################################################################################################

# setwd("")

rdata<-read.csv("FIARRSTAT_recruitment_data_stat_final_150715.csv",header=T)
rdata<-rdata[,c(-1)]

rdata<-rdata[rdata$stdage>=0&!is.na(rdata$stdage),]
rdata_fp_pcn<-unique(rdata[rdata$FIX==1,]$pcn) #fixer present
rdata_fp<-data.frame(pcn=rdata_fp_pcn,fixer_present=1)
rdata_fp1<-merge(rdata_fp,rdata,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=TRUE)
rdata_fp2<-rdata_fp1[!is.na(rdata_fp1$fixer_present),]

rdata.original<-rdata
rdata<-rdata_fp2

# setwd("")

pts<-read.csv("plot_ts.csv",header=T)
p<-read.csv("plot.csv",header=T)
p1<-data.frame(pcn0=pts$pcn1,pcnT=pts$pcn2)
p2<-data.frame(pcn0=pts$pcn2,pcnT=pts$pcn3)
p3<-data.frame(pcn0=pts$pcn3,pcnT=pts$pcn4)
pcn<-rbind(p1,p2,p3)
pcn<-pcn[!is.na(pcn$pcn0)&!is.na(pcn$pcnT),]

pcn0<-data.frame(pcn0=pcn$pcn0)
pcnT<-data.frame(pcnT=pcn$pcnT)

rdata_original<-rdata

# setwd("")

tcn_cclcd<-read.csv("FIA_tcn_cclcd_150715.csv",header=T)
tcn_cclcd<-tcn_cclcd[,c(-1)]
rdata_cclcd<-merge(rdata,tcn_cclcd,by.x="tcn",by.y="tcn",all.x=FALSE,all.y=FALSE)

# setwd("")

write.csv(rdata_cclcd,file="FIA_RecruitmentDataforFit_wCCLCD_150722.csv")

######################################################
#############     Mean Calculations    ###############
######################################################

### Canopy

rdata<-rdata_cclcd[rdata_cclcd$cclcd==1 | rdata_cclcd$cclcd==2 | rdata_cclcd$cclcd==3,]

rdata0<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcn0,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=FALSE)
colnames(rdata0)<-c("pcn0","spcd_0","stdage_0","tcn_0","R_0","dbh_0","N_0","FIX_0")
rdataT<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcnT,by.x="pcn",by.y="pcnT",all.x=FALSE,all.y=FALSE)
colnames(rdataT)<-c("pcnT","spcd_T","stdage_T","tcn_T","R_T","dbh_T","N_T","FIX_T")

sum.na<-function(x){
	x<-x[!is.na(x)]
	sum(x)
}

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

rdata.original<-rdata
# For pcn0: time 0

rdata<-rdata0

nonN<-tapply(rdata[rdata$FIX_0==0,]$N_0,rdata[rdata$FIX_0==0,]$pcn0,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_0==0,]$dbh_0,rdata[rdata$FIX_0==0,]$pcn0,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na) # for time0, calculate mean diameter for all trees

fixN<-tapply(rdata[rdata$FIX_0==1,]$N_0,rdata[rdata$FIX_0==1,]$pcn0,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_0==1,]$dbh_0,rdata[rdata$FIX_0==1,]$pcn0,mean.na)
fixmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na)

n1<-data.frame(pcn0=as.numeric(rownames(nonN)),N=nonN)
n4<-data.frame(pcn0=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcn0=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcn0=as.numeric(rownames(fixN)),N=fixN)
f4<-data.frame(pcn0=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcn0=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)

fix_0<-mf2
fix_0$FIX<-1
non_0<-mn2
non_0$FIX<-0
final_0<-rbind(fix_0,non_0)
final_0<-final_0[!is.na(final_0$N)&!is.na(final_0$pcn0)&!is.na(final_0$mDBH)&!is.na(final_0$mDBHbt)&!is.na(final_0$FIX),]

# For pcnT: time T

rdata<-rdataT

nonR<-tapply(rdata[rdata$FIX_T==0,]$R_T,rdata[rdata$FIX_T==0,]$pcnT,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_T==0,]$dbh_T,rdata[rdata$FIX_T==0,]$pcnT,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na) # for time0, calculate mean diameter for all trees

fixR<-tapply(rdata[rdata$FIX_T==1,]$R_T,rdata[rdata$FIX_T==1,]$pcnT,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_T==1,]$dbh_T,rdata[rdata$FIX_T==1,]$pcnT,mean.na)
fixmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na)

n1<-data.frame(pcnT=as.numeric(rownames(nonR)),R=nonR)
n4<-data.frame(pcnT=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcnT=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcnT=as.numeric(rownames(fixR)),R=fixR)
f4<-data.frame(pcnT=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcnT=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)

fix_T<-mf2
fix_T$FIX<-1
non_T<-mn2
non_T$FIX<-0
final_T<-rbind(fix_T,non_T)
final_T<-final_T[!is.na(final_T$R)&!is.na(final_T$pcnT)&!is.na(final_T$mDBH)&!is.na(final_T$mDBHbt)&!is.na(final_T$FIX),]

pcn1<-merge(pcn,final_0,by.x="pcn0",by.y="pcn0",all.x=FALSE,all.y=FALSE)
pcn2<-merge(pcn1,final_T[,c(1,2,5)],by=c("pcnT","FIX"),all.x=FALSE,all.y=FALSE)

final<-pcn2

fp<-p[,c(1,2,5,6,8)]
MRdat<-merge(fp,final,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

mt<-data.frame(pcn=p$pcn,mt=p$meastime)
MRdat_mt0<-merge(MRdat,mt,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mt0)[12]<-"mt0"
MRdat_mtT<-merge(MRdat_mt0,mt,by.x="pcnT",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mtT)[13]<-"mtT"

MRdat_mtT$mt<-MRdat_mtT$mtT-MRdat_mtT$mt0
MRdat<-MRdat_mtT
MRdat$RR<-NA
MRdat<-MRdat[MRdat$N>0,]
MRdat$RR<-MRdat$R/(MRdat$N*MRdat$mt)

MRdat<-MRdat[MRdat$stdage>=0,]
MRdat_canopy<-MRdat

############ Canopy
MRdat<-MRdat_canopy

fixer_rrall<-MRdat[MRdat$FIX==1,]
nonfixer_rrall<-MRdat[MRdat$FIX==0,]

avg<-tapply(MRdat$RR,MRdat$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
allavg_rr<-favg
plot(allavg_rr$meanRR~allavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. All Trees")

avg<-tapply(MRdat[MRdat$FIX==1,]$RR,MRdat[MRdat$FIX==1,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
fxavg_rr<-favg
plot(fxavg_rr$meanRR~fxavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Fixers")

avg<-tapply(MRdat[MRdat$FIX==0,]$RR,MRdat[MRdat$FIX==0,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
navg_rr<-favg
plot(navg_rr$meanRR~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Non Fixers")

plot(navg_rr$meanRR*100~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",lwd=1,col="blue",main="U.S. Trees",pch=2,ylim=c(0,100))
points(fxavg_rr$meanRR*100~fxavg_rr$stdage,lwd=1,col="red",pch=1)
legend("topright",c("N Fixers","Non Fixers"),col=c("red","blue"),pch=c(1,2),bty="n")

par(mfrow=c(1,2))
hist(navg_rr$meanRR,main="Non Fixer Recruitment Rate",xlab="Recruitment Rate ",breaks=seq(0,0.8,0.01))
hist(fxavg_rr$meanRR,main="N Fixer Recruitment Rate",xlab="Recruitment Rate ",breaks=seq(0,0.9,0.01))

navg_rr_canopy<-navg_rr
fxavg_rr_canopy<-fxavg_rr

### Understory

rdata<-rdata_cclcd[rdata_cclcd$cclcd==4 | rdata_cclcd$cclcd==5,]

rdata0<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcn0,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=FALSE)
colnames(rdata0)<-c("pcn0","spcd_0","stdage_0","tcn_0","R_0","dbh_0","N_0","FIX_0")
rdataT<-merge(rdata[,c("pcn","spcd0","stdage","tcn","R","dbh0","N","FIX")],pcnT,by.x="pcn",by.y="pcnT",all.x=FALSE,all.y=FALSE)
colnames(rdataT)<-c("pcnT","spcd_T","stdage_T","tcn_T","R_T","dbh_T","N_T","FIX_T")

sum.na<-function(x){
	x<-x[!is.na(x)]
	sum(x)
}

mean.na<-function(x){
	x<-x[!is.na(x)]
	mean(x)
}

rdata.original<-rdata
# For pcn0: time 0

rdata<-rdata0

nonN<-tapply(rdata[rdata$FIX_0==0,]$N_0,rdata[rdata$FIX_0==0,]$pcn0,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_0==0,]$dbh_0,rdata[rdata$FIX_0==0,]$pcn0,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na) # for time0, calculate mean diameter for all trees

fixN<-tapply(rdata[rdata$FIX_0==1,]$N_0,rdata[rdata$FIX_0==1,]$pcn0,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_0==1,]$dbh_0,rdata[rdata$FIX_0==1,]$pcn0,mean.na)
fixmDBHbt<-tapply(rdata$dbh_0,rdata$pcn0,mean.na)

n1<-data.frame(pcn0=as.numeric(rownames(nonN)),N=nonN)
n4<-data.frame(pcn0=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcn0=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcn0=as.numeric(rownames(fixN)),N=fixN)
f4<-data.frame(pcn0=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcn0=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcn0",by.y="pcn0",all.x=TRUE,all.y=TRUE)

fix_0<-mf2
fix_0$FIX<-1
non_0<-mn2
non_0$FIX<-0
final_0<-rbind(fix_0,non_0)
final_0<-final_0[!is.na(final_0$N)&!is.na(final_0$pcn0)&!is.na(final_0$mDBH)&!is.na(final_0$mDBHbt)&!is.na(final_0$FIX),]

# For pcnT: time T

rdata<-rdataT

nonR<-tapply(rdata[rdata$FIX_T==0,]$R_T,rdata[rdata$FIX_T==0,]$pcnT,sum.na) # for time0, calculate toal number alive
nonmDBH<-tapply(rdata[rdata$FIX_T==0,]$dbh_T,rdata[rdata$FIX_T==0,]$pcnT,mean.na) # for time0, calculate mean diameter for trees in the same functional group
nonmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na) # for time0, calculate mean diameter for all trees

fixR<-tapply(rdata[rdata$FIX_T==1,]$R_T,rdata[rdata$FIX_T==1,]$pcnT,sum.na)
fixmDBH<-tapply(rdata[rdata$FIX_T==1,]$dbh_T,rdata[rdata$FIX_T==1,]$pcnT,mean.na)
fixmDBHbt<-tapply(rdata$dbh_T,rdata$pcnT,mean.na)

n1<-data.frame(pcnT=as.numeric(rownames(nonR)),R=nonR)
n4<-data.frame(pcnT=as.numeric(rownames(nonmDBH)),mDBH=nonmDBH)
n5<-data.frame(pcnT=as.numeric(rownames(nonmDBHbt)),mDBHbt=nonmDBHbt)

f1<-data.frame(pcnT=as.numeric(rownames(fixR)),R=fixR)
f4<-data.frame(pcnT=as.numeric(rownames(fixmDBH)),mDBH=fixmDBH)
f5<-data.frame(pcnT=as.numeric(rownames(fixmDBHbt)),mDBHbt=fixmDBHbt)

mf1<-merge(f1,f4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mf2<-merge(mf1,f5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn1<-merge(n1,n4,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)
mn2<-merge(mn1,n5,by.x="pcnT",by.y="pcnT",all.x=TRUE,all.y=TRUE)

fix_T<-mf2
fix_T$FIX<-1
non_T<-mn2
non_T$FIX<-0
final_T<-rbind(fix_T,non_T)
final_T<-final_T[!is.na(final_T$R)&!is.na(final_T$pcnT)&!is.na(final_T$mDBH)&!is.na(final_T$mDBHbt)&!is.na(final_T$FIX),]

pcn1<-merge(pcn,final_0,by.x="pcn0",by.y="pcn0",all.x=FALSE,all.y=FALSE)
pcn2<-merge(pcn1,final_T[,c(1,2,5)],by=c("pcnT","FIX"),all.x=FALSE,all.y=FALSE)

final<-pcn2

fp<-p[,c(1,2,5,6,8)]
MRdat<-merge(fp,final,by.x="pcn",by.y="pcn0",all.x=FALSE,all.y=TRUE)

mt<-data.frame(pcn=p$pcn,mt=p$meastime)
MRdat_mt0<-merge(MRdat,mt,by.x="pcn",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mt0)[12]<-"mt0"
MRdat_mtT<-merge(MRdat_mt0,mt,by.x="pcnT",by.y="pcn",all.x=TRUE,all.y=FALSE)
colnames(MRdat_mtT)[13]<-"mtT"

MRdat_mtT$mt<-MRdat_mtT$mtT-MRdat_mtT$mt0
MRdat<-MRdat_mtT
MRdat$RR<-NA
MRdat<-MRdat[MRdat$N>0,]
MRdat$RR<-MRdat$R/(MRdat$N*MRdat$mt)

MRdat<-MRdat[MRdat$stdage>=0,]
MRdat_understory<-MRdat

MRdat<-MRdat_understory

fixer_rrall<-MRdat[MRdat$FIX==1,]
nonfixer_rrall<-MRdat[MRdat$FIX==0,]

avg<-tapply(MRdat$RR,MRdat$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
allavg_rr<-favg
plot(allavg_rr$meanRR~allavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. All Trees")

avg<-tapply(MRdat[MRdat$FIX==1,]$RR,MRdat[MRdat$FIX==1,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
fxavg_rr<-favg
plot(fxavg_rr$meanRR~fxavg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Fixers")

avg<-tapply(MRdat[MRdat$FIX==0,]$RR,MRdat[MRdat$FIX==0,]$stdage,mean.na)
tavg<-t(avg)
tavg<-t(tavg)
favg<-data.frame(stdage=as.numeric(rownames(tavg)),meanRR=tavg)
navg_rr<-favg
plot(navg_rr$meanRR~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",main="U.S. Non Fixers")

plot(navg_rr$meanRR*100~navg_rr$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",lwd=1,col="blue",main="U.S. Trees",pch=2,ylim=c(0,100))
points(fxavg_rr$meanRR*100~fxavg_rr$stdage,lwd=1,col="red",pch=1)
legend("topright",c("N Fixers","Non Fixers"),col=c("red","blue"),pch=c(1,2),bty="n")

par(mfrow=c(1,2))
hist(navg_rr$meanRR,main="Non Fixer Recruitment Rate",xlab="Recruitment Rate ")
hist(fxavg_rr$meanRR,main="N Fixer Recruitment Rate",xlab="Recruitment Rate ")

navg_rr_understory<-navg_rr
fxavg_rr_understory<-fxavg_rr

######################################################
###############     STAT Analysis     ################
######################################################
library(bbmle)

### Canopy

MRdat<-MRdat_canopy
navg_rr<-navg_rr_canopy
fxavg_rr<-fxavg_rr_canopy

################ Saturating

fn1_rec_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
	y <- y*t
} # Saturating with diameter

# fixers

# MLE
pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_f<-pv

NLL2_rec_d_t_f <- function(af,bf,cf,df,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	yf <- fn1_rec_d_t(af, bf, cf, df, D[F==1], A[F==1], Dhat, t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit2_rec_d_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_f)

fit2_rec_d_c_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_f)

# non-fixers

pv<-c(-20,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_n<-pv

NLL2_rec_d_t_n <- function(a0,b0,c0,d0,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	y0 <- fn1_rec_d_t(a0, b0, c0, d0, D[F==0], A[F==0], Dhat, t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit2_rec_d_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_n)

fit2_rec_d_c_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_n)

################ Ricker

fn4_rec_d_t<-function(a,c,d,k,D,A,Dhat,t){
		a<-exp(a)
	c<-exp(c)
y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	y<-y*t
	} #Right-skewed with diameter effect

# fixers

# MLE

pv<-c(-4,-8.6,-3.7,-1.6,-0.6,-0.2,0.2,0.18)
pv_f<-pv

NLL5_rec_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	yf<-fn4_rec_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit5_rec_d_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_f)

fit5_rec_d_c_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_f)

# non-fixers

# MLE
pv<-c(-6,-8.6,3,-1.6,-0.6,-0.2,0.2,0.18)

NLL5_rec_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	y0<-fn4_rec_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit5_rec_d_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_n)

fit5_rec_d_c_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_n)

######################################################
################   Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

# fixer
deltaAICs(c(AIC(fit2_rec_d_t_f),AIC(fit2_rec_d_c_t_f),AIC(fit5_rec_d_t_f),AIC(fit5_rec_d_c_t_f)))

# [1] 1.477254
# [1] 0
# [1] 1.614316
# [1] 0.5074987

# non-fixers
deltaAICs(c(AIC(fit2_rec_d_t_n),AIC(fit2_rec_d_c_t_n),AIC(fit5_rec_d_t_n),AIC(fit5_rec_d_c_t_n)))
# [1] 3.834806
# [1] 0
# [1] 105.2227
# [1] 91.0209

pv_f_rr_canopy<-coef(fit2_rec_d_c_t_f)
pv_n_rr_canopy<-coef(fit2_rec_d_c_t_n)

######################################################
#################    CI Computation   ################
######################################################

library(MASS)

agevec_rr_canopy<-unique(MRdat_canopy$stdage)
mD<-mean(MRdat_canopy$mDBH)
mD_fix<-mean(MRdat_canopy[MRdat_canopy$FIX==1,]$mDBH)
mD_non<-mean(MRdat_canopy[MRdat_canopy$FIX==0,]$mDBH)

ndraws<-1000
xs_rr_canopy<-seq(min(agevec_rr_canopy),max(agevec_rr_canopy),1)
xs<-seq(min(agevec_rr_canopy),max(agevec_rr_canopy),1)

vmat_fix_rr<-mvrnorm(ndraws,mu=pv_f_rr_canopy,Sigma=vcov(fit2_rec_d_c_t_f)) 
vmat_non_rr<-mvrnorm(ndraws,mu=pv_n_rr_canopy,Sigma=vcov(fit2_rec_d_c_t_n))

ys_fix_rr<-array(dim=c(ndraws,length(xs)))
ys_non_rr<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix_rr[i,]<-fn1_rec_d_t(vmat_fix_rr[i,1],vmat_fix_rr[i,2],vmat_fix_rr[i,3],vmat_fix_rr[i,4],mD_fix,xs,mD,1)*100
	ys_non_rr[i,]<-fn1_rec_d_t(vmat_non_rr[i,1],vmat_non_rr[i,2],vmat_non_rr[i,3],vmat_non_rr[i,4],mD_non,xs,mD,1)*100
}

civec_fix_rr<-array(dim=c(2,length(xs)))
civec_non_rr<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix_rr[,j]<-quantile(ys_fix_rr[,j],c(lowb,upb),na.rm=TRUE)
	civec_non_rr[,j]<-quantile(ys_non_rr[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_rr_canopy<-civec_fix_rr
civec_non_rr_canopy<-civec_non_rr

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# setwd("")

pdf("Fig16_RRfit_wCI_Canopy_150930.pdf",width=4.5,height=4.5)

x<-c(0,180)
y<-c(0,60)
plot(y~x,type="n",xlim=c(0,180),y=c(0,60),xlab="Stand Age (Years)",ylab="Recruitment Rate (% Recruited/year)",main="Canopy")
points(navg_rr_canopy$meanRR*100~navg_rr_canopy$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",col="blue",main="U.S. Trees",pch=2,cex=0.8)
points(fxavg_rr_canopy$meanRR*100~fxavg_rr_canopy$stdage,col="red",pch=1,cex=0.8)

# fixers
	# Saturating
curve(fn1_rec_d_t(pv_f_rr_canopy[1],pv_f_rr_canopy[2],pv_f_rr_canopy[3],pv_f_rr_canopy[4],mD_fix,x,mD,1)*100,from=min(agevec_rr_canopy),to=max(agevec_rr_canopy),lwd=3,lty=1,add=TRUE,col="red")

# non-fixers
	# Saturating
curve(fn1_rec_d_t(pv_n_rr_canopy[1],pv_n_rr_canopy[2],pv_n_rr_canopy[3],pv_n_rr_canopy[4],mD_non,x,mD,1)*100,from=min(agevec_rr_canopy),to=max(agevec_rr_canopy),lwd=3,col="blue",lty=1,add=TRUE,xlab="Stand Age",ylab="Estimated Relative Recruitment Rate (%)")

legend("topright",c("Canopy N fixers","Canopy Non-fixers"),lty=1,bty="n",cex=1.2,lwd=2,col=c("red","blue"),pch=c(1,2))

points(civec_fix_rr_canopy[1,]~xs_rr_canopy,lwd=2,typ="l",lty=2,col="red")
points(civec_fix_rr_canopy[2,]~xs_rr_canopy,lwd=2,typ="l",lty=2,col="red")

points(civec_non_rr_canopy[1,]~xs_rr_canopy,lwd=2,typ="l",lty=2,col="blue") # This CI is really wierd. Need to troubleshoot.
points(civec_non_rr_canopy[2,]~xs_rr_canopy,lwd=2,typ="l",lty=2,col="blue")

dev.off()

### Understory

MRdat<-MRdat_understory
navg_rr<-navg_rr_understory
fxavg_rr<-fxavg_rr_understory

################ Saturating

fn1_rec_d_t <- function(a, b, c, d, D, A, Dhat, t){
	a<-exp(a)
	b<-exp(b)
y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
	y <- y*t
} # Saturating with diameter

# fixers

# MLE
pv<-c(-3,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_f<-pv

NLL2_rec_d_t_f <- function(af,bf,cf,df,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	yf <- fn1_rec_d_t(af, bf, cf, df, D[F==1], A[F==1], Dhat, t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit2_rec_d_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_f)

fit2_rec_d_c_t_f <- mle2(NLL2_rec_d_t_f,start=list(af=pv_f[2],bf=pv_f[4],cf=pv_f[6],df=pv_f[8]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_f)

# non-fixers

pv<-c(-30,-6.5,-1.2,-0.4,10,8,-0.59,-0.2)
pv_n<-pv

NLL2_rec_d_t_n <- function(a0,b0,c0,d0,N,R,D,A,F,t){
	Dhat <- mean(MRdat$mDBH)
	y0 <- fn1_rec_d_t(a0, b0, c0, d0, D[F==0], A[F==0], Dhat, t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit2_rec_d_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBH,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_t_n)

fit2_rec_d_c_t_n <- mle2(NLL2_rec_d_t_n,start=list(a0=pv_n[1],b0=pv_n[3],c0=pv_n[5],d0=pv_n[7]),
	data=list(N=MRdat$N,R=MRdat$R,D=MRdat$mDBHbt,A=MRdat$stdage,F=MRdat$FIX,t=MRdat$mt))
summary(fit2_rec_d_c_t_n)

################ Ricker

fn4_rec_d_t<-function(a,c,d,k,D,A,Dhat,t){
		a<-exp(a)
	c<-exp(c)
y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	y<-y*t
	} #Right-skewed with diameter effect

# fixers

# MLE

pv<-c(-4,-8.6,-3.7,-4.6,-0.6,-0.2,0.2,0.18)
pv_f<-pv

NLL5_rec_d_t_f<-function(af,cf,df,kf,N,S,D,A,F,t){
	yf<-fn4_rec_d_t(af,cf,df,kf,D[F==1],A[F==1],exp(mean(log(D))),t[F==1])
- sum(dpois(R[F==1],lambda=yf*N[F==1],log=TRUE))
}

fit5_rec_d_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_f)

fit5_rec_d_c_t_f<-mle2(NLL5_rec_d_t_f,start=list(af=pv[2],cf=pv[4],df=pv[6],kf=pv[8]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_f)

# non-fixers

# MLE
pv<-c(-6,-8.6,2,-1.6,-0.6,-0.2,0.2,0.18)

NLL5_rec_d_t_n<-function(a0,c0,d0,k0,N,S,D,A,F,t){
	y0<-fn4_rec_d_t(a0,c0,d0,k0,D[F==0],A[F==0],exp(mean(log(D))),t[F==0])
	-sum(dpois(R[F==0],lambda=y0*N[F==0],log=TRUE))
}

fit5_rec_d_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBH,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_t_n)

fit5_rec_d_c_t_n<-mle2(NLL5_rec_d_t_n,start=list(a0=pv[1],c0=pv[3],d0=pv[5],k0=pv[7]),data=list(N=MRdat$N,R=MRdat$R,A=MRdat$stdage,F=MRdat$FIX,D=MRdat$mDBHbt,Dhat=mean(MRdat$mDBH),t=MRdat$mt))
summary(fit5_rec_d_c_t_n)

######################################################
################   Model Comparison   ################
######################################################

deltaAICs <- function(AICvec){
	for(i in 1:length(AICvec)){
		print(AICvec[i]-min(AICvec))
	}
}

# fixer
deltaAICs(c(AIC(fit2_rec_d_t_f),AIC(fit2_rec_d_c_t_f),AIC(fit5_rec_d_t_f),AIC(fit5_rec_d_c_t_f)))

# [1] 0.3482467
# [1] 0.1898743
# [1] 0.1898556
# [1] 0

# non-fixers
deltaAICs(c(AIC(fit2_rec_d_t_n),AIC(fit2_rec_d_c_t_n),AIC(fit5_rec_d_t_n),AIC(fit5_rec_d_c_t_n)))
# [1] 0
# [1] 1.301578
# [1] 7.910167
# [1] 8.817104

pv_f_rr_understory<-coef(fit5_rec_d_c_t_f)
pv_n_rr_understory<-coef(fit2_rec_d_t_n)

######################################################
#################    CI Computation   ################
######################################################

library(MASS)

agevec_rr_understory<-unique(MRdat_understory$stdage)
mD<-mean(MRdat_understory$mDBH)
mD_fix<-mean(MRdat_understory[MRdat_understory$FIX==1,]$mDBH)
mD_non<-mean(MRdat_understory[MRdat_understory$FIX==0,]$mDBH)

ndraws<-10000
xs_rr_understory<-seq(min(agevec_rr_understory),max(agevec_rr_understory),1)
xs<-seq(min(agevec_rr_understory),max(agevec_rr_understory),1)

vmat_fix_rr<-mvrnorm(ndraws,mu=pv_f_rr_understory,Sigma=vcov(fit5_rec_d_c_t_f)) 
vmat_non_rr<-mvrnorm(ndraws,mu=pv_n_rr_understory,Sigma=vcov(fit2_rec_d_t_n))

ys_fix_rr<-array(dim=c(ndraws,length(xs)))
ys_non_rr<-array(dim=c(ndraws,length(xs)))

for(i in 1:ndraws){
	print(i/ndraws*100)
	ys_fix_rr[i,]<-fn4_rec_d_t(vmat_fix_rr[i,1],vmat_fix_rr[i,2],vmat_fix_rr[i,3],vmat_fix_rr[i,4],mD_fix,xs,mD,1)*100
	ys_non_rr[i,]<-fn1_rec_d_t(vmat_non_rr[i,1],vmat_non_rr[i,2],vmat_non_rr[i,3],vmat_non_rr[i,4],mD_non,xs,mD,1)*100
}

civec_fix_rr<-array(dim=c(2,length(xs)))
civec_non_rr<-array(dim=c(2,length(xs)))

lowb<-0.025
upb<-0.975
for(j in 1:length(xs)){
	civec_fix_rr[,j]<-quantile(ys_fix_rr[,j],c(lowb,upb),na.rm=TRUE)
	civec_non_rr[,j]<-quantile(ys_non_rr[,j],c(lowb,upb),na.rm=TRUE)
}

civec_fix_rr_understory<-civec_fix_rr
civec_non_rr_understory<-civec_non_rr

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

# setwd("")

pdf("Fig17_RRfit_wCI_Understory_150930.pdf",width=4.5,height=4.5)

x<-c(0,180)
y<-c(0,100)
plot(y~x,type="n",xlim=c(0,180),y=c(0,100),xlab="Stand Age (Years)",ylab="Recruitment Rate (% Recruited/year)",main="Understory")
points(navg_rr_understory$meanRR*100~navg_rr_understory$stdage,xlab="Stand Age",ylab="Mean Recruitment at Each Stand Age",col="blue",main="U.S. Trees",pch=2,cex=0.8)
points(fxavg_rr_understory$meanRR*100~fxavg_rr_understory$stdage,col="red",pch=1,cex=0.8)

# fixers
	# Ricker
curve(fn4_rec_d_t(pv_f_rr_understory[1],pv_f_rr_understory[2],pv_f_rr_understory[3],pv_f_rr_understory[4],mD_fix,x,mD,1)*100,from=min(agevec_rr_understory),to=max(agevec_rr_understory),lwd=3,lty=1,add=TRUE,col="red")

# non-fixers
	# Sat
curve(fn1_rec_d_t(pv_n_rr_understory[1],pv_n_rr_understory[2],pv_n_rr_understory[3],pv_n_rr_understory[4],mD_non,x,mD,1)*100,from=min(agevec_rr_understory),to=max(agevec_rr_understory),lwd=3,col="blue",lty=1,add=TRUE,xlab="Stand Age",ylab="Estimated Relative Recruitment Rate (%)")

legend("topright",c("Understory N fixers","Understory Non-fixers"),lty=1,bty="n",cex=1.2,lwd=2,col=c("red","blue"),pch=c(1,2))

points(civec_fix_rr_understory[1,]~xs_rr_understory,lwd=2,typ="l",lty=2,col="red")
points(civec_fix_rr_understory[2,]~xs_rr_understory,lwd=2,typ="l",lty=2,col="red")

points(civec_non_rr_understory[1,]~xs_rr_understory,lwd=2,typ="l",lty=2,col="blue") 
points(civec_non_rr_understory[2,]~xs_rr_understory,lwd=2,typ="l",lty=2,col="blue")

dev.off()





################################
################################
################################
################################
################################
################################
################################
################################
################################








####################################################################################################
################################                    IBM               ##############################
####################################################################################################

### To increase computational efficiency, when # of alive trees exceeds 10000, we reduce population size to 30% for each functional group.

### 10 runs for each IBM scenario

#########################################################
##############     Load in IBM Data    ################## # All
#########################################################

# setwd("")

gdata<-read.csv("FIAGRSTAT_final_data_for_fit_150707.csv",header=T)
gdata<-gdata[,c(-1)]

#########################################################
########   Load in Parameters/Functions    ############## # All 
#########################################################

### Growth
# fixer fit_logGR_4_f
# non-fixer fit_logGR_1_n

pv_f_gr_all<-c(2.43991324,-7.47700936,-2.25434542,-0.55436896,0.02089067) # sd,a,c,d,k
pv_n_gr_all<-c(1.422643623, -5.636759485, -4.285533715, 48.22197703, -0.491180703) # sd,a,b,c,d

gMLF_4_f_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

gMLF_1_n_Dhat<-function(a,b,c,d,D,A,Dhat){
	a<-exp(a)
	b<-exp(b)
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

### Mortality
# fixer fit_ricker_surv_mi_d_t_f
# non-fixer fit_ricker_surv_mi_d_t_n

pv_f_mr_all<-c(-3.6843464,-3.3737355,-0.3540852,0.1018018) #af,cf,df,kf
pv_n_mr_all<-c(-4.16235239,-3.04883890,-0.62008068,0.02406838) # a0,c0,d0,k0

fn_ricker_surv_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m) # This is the survival probability
}

fn_ricker_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

### Recruitment
# fixer fit5_rec_d_c_t_f
# non-fixer fit2_rec_d_c_t_n

pv_f_rr_all<-c(-4.2749295,-1.4506284,1.4402179,0.7805302) #af,cf,df,kf
pv_n_rr_all<-c(-3.4848887,0.1021638,6.0737946,0.3888109) #a0,b0,c0,d0

fn4_rec_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	} #Right-skewed with diameter effect # got rid of measurement time, now it's annual rate

fn1_rec_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
} # Saturating with diameter # got rid of measurement time, now it's annual rate

#########################################################
##############           Runs          ################## # Gr Mr Rr Diff
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-1
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=startage:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")

write.csv(looptrack_runs,file="FIA_IBM_GrMrRrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Gr Diff
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 0
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-1
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=startage:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")

write.csv(looptrack_runs,file="FIA_IBM_GrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Mr Diff
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-1
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=startage:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")

write.csv(looptrack_runs,file="FIA_IBM_MrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Rr Diff
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 0
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-1
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=startage:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_RrDiff_Updated_150729.csv")

#########################################################
##############     Load in IBM Data    ################## # Robinia
#########################################################

# setwd("")
gdata<-read.csv("FIA_GrowthData_RobiniaPseudoacacia_150709.csv",header=T)
gdata<-gdata[,c(-1)]

#########################################################
########   Load in Parameters/Functions    ############## # Robinia
#########################################################

### Growth
# fixer fit_logGR_4_f
# non-fixer fit_logGR_4_n

pv_f_gr_all<-c(2.483432016, -7.803591349, -2.285318111, -0.444301579, 0.019968956) # sd,a,c,d,k
pv_n_gr_all<-c(1.410749312, -5.230214624, -2.572961848, -0.475666825, 0.024059197) # sd,a,c,d,k

gMLF_4_f_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

gMLF_4_n_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

### Mortality
# fixer fit_ricker_surv_mi_d_t_f
# non-fixer fit_sat_surv_mi_d_t_n

pv_f_mr_all<-c(-2.700711528, -3.23250364, -0.531328682, 0.011390397) #af,cf,df,kf
pv_n_mr_all<-c(-4.118112193, -2.284634107, 2.961686703, -0.713159055) # a,b,c,d

fn_ricker_surv_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m) # This is the survival probability
}

fn_ricker_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

fn_sat_surv_mort_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
	gamma<-exp(-m)
} # saturating with diameter effect

fn_sat_mort_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	m <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d 
}

### Recruitment
# fixer fit5_rec_d_c_t_f
# non-fixer fit2_rec_d_c_t_n

pv_f_rr_all<-c(-4.267900034, -1.44177184, 1.376702378, 0.809134683) #af,cf,df,kf
pv_n_rr_all<-c(-3.47517344, 1.185911335, 4.68522797, 0.293531281) #a0,b0,c0,d0

fn4_rec_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	} #Right-skewed with diameter effect # got rid of measurement time, now it's annual rate

fn1_rec_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
} # Saturating with diameter # got rid of measurement time, now it's annual rate

#########################################################
##############           Runs          ################## # Gr Mr Rr Diff # Robinia
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-300

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Robinia_GrMrRrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Gr Diff # Robinia
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 0
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Robinia_GrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Mr Diff # Robinia
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Robinia_MrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Rr Diff # Robinia
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 0
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=10,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 28 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 135 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.2),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_4_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_sat_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Robinia_RrDiff_Updated_150729.csv")

#########################################################
##############     Load in IBM Data    ################## # Cercocarpus
#########################################################

# setwd("")
gdata<-read.csv("FIA_GrowthData_CercocarpusLedifolius_150709.csv",header=T)
gdata<-gdata[,c(-1)]

#########################################################
########   Load in Parameters/Functions    ############## # Cercocarpus 
#########################################################

### Growth
# fixer fit_logGR_2_f
# non-fixer fit_logGR_1_n

pv_f_gr_all<-c(2.434070539, -7.791209655, -3.352334228, 15.12812733, -0.805774517, -0.043412716) # sd,a,b,c,d,k
pv_n_gr_all<-c(1.346496925, -5.569811309, -26.74682926, 27.11697916, -0.928119366) # sd,a,b,c,d

gMLF_2_f_Dhat<-function(a,b,c,d,k,D,A,Dhat){
	a<-exp(a)
	b<-exp(b)
	y<-log(a+(b-a)/(1+exp(-k*(A-c))))+d*(log(D)-log(Dhat))
}

gMLF_1_n_Dhat<-function(a,b,c,d,D,A,Dhat){
	a<-exp(a)
	b<-exp(b)
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

### Mortality
# fixer fit_ricker_surv_mi_d_t_f
# non-fixer fit_ricker_surv_mi_d_t_n

pv_f_mr_all<-c(-7.75403971, -10.07565381, 0.285573233, 3.899748963) #af,cf,df,kf
pv_n_mr_all<-c(-3.377629851, -3.787105095, -0.298348607, -0.030641217) # a0,c0,d0,k0

fn_ricker_surv_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m) # This is the survival probability
}

fn_ricker_mort_d_t_Dhat <-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

### Recruitment
# fixer fit5_rec_d_c_t_f
# non-fixer fit2_rec_d_c_t_n

pv_f_rr_all<-c(-4.2749295,-1.4506284,1.4402179,0.7805302) #af,cf,df,kf
pv_n_rr_all<-c(-3.4848887,0.1021638,6.0737946,0.3888109) #a0,b0,c0,d0

fn4_rec_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	} #Right-skewed with diameter effect # got rid of measurement time, now it's annual rate

fn1_rec_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
} # Saturating with diameter # got rid of measurement time, now it's annual rate

#########################################################
##############           Runs          ################## # Gr Mr Rr Diff # Cercocarpus 
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-300

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_2_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],gvec_fix[6],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Cercocarpus_GrMrRrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Gr Diff # Cercocarpus 
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 0
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-300

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_2_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],gvec_fix[6],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Cercocarpus_GrDiff_150930.csv")

#########################################################
##############           Runs          ################## # Mr Diff # Cercocarpus 
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-300

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_2_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],gvec_fix[6],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Cercocarpus_MrDiff_150930.csv")

#########################################################
##############           Runs          ################## # Gr Mr Diff # Cercocarpus 
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-300

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_2_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],gvec_fix[6],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_1_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Cercocarpus_GrMrDiff_150930.csv")

#########################################################
##############     Load in IBM Data    ################## # Alnus
#########################################################

# setwd("")
gdata<-read.csv("FIA_GrowthData_AlnusRubra_150709.csv",header=T)
gdata<-gdata[,c(-1)]

#########################################################
########   Load in Parameters/Functions    ############## # Alnus 
#########################################################

### Growth
# fixer fit_logGR_1_f
# non-fixer fit_logGR_4_n

pv_f_gr_all<-c(1.610766808, -6.489136908, -2.847599361, 14.40147934, -0.350026919) # sd,a,b,c,d
pv_n_gr_all<-c(1.952419097, -9.264741841, -3.99129076, -0.334999839, 0.01339732) # sd,a,c,d,k

gMLF_1_f_Dhat<-function(a,b,c,d,D,A,Dhat){
	a<-exp(a)
	b<-exp(b)
	y<-log(a+(b-a)*exp(-A/c))+d*(log(D)-log(Dhat)) 
}

gMLF_4_n_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-log(a+c*(k-a)*A*exp(-c*A))+d*(log(D)-log(Dhat)) 
}

### Mortality
# fixer fit_ricker_surv_mi_d_t_f
# non-fixer fit_ricker_surv_mi_d_t_n

pv_f_mr_all<-c(-3.662369213, -8.105139266, -1.033022901, 0.808692778) #af,cf,df,kf
pv_n_mr_all<-c(-2.85773833, -3.971103717, -0.653392455, -0.049975799) # a0,c0,d0,k0

fn_ricker_surv_mort_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	gamma<-exp(-m) # This is the survival probability
}

fn_ricker_mort_d_t_Dhat <-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	m<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
}

### Recruitment
# fixer fit5_rec_d_c_t_f
# non-fixer fit2_rec_d_c_t_n

pv_f_rr_all<-c(-4.2749295,-1.4506284,1.4402179,0.7805302) #af,cf,df,kf
pv_n_rr_all<-c(-3.4848887,0.1021638,6.0737946,0.3888109) #a0,b0,c0,d0

fn4_rec_d_t_Dhat<-function(a,c,d,k,D,A,Dhat){
	a<-exp(a)
	c<-exp(c)
	y<-(a+c*(k-a)*A*exp(-c*A))*(D/Dhat)^d
	} #Right-skewed with diameter effect # got rid of measurement time, now it's annual rate

fn1_rec_d_t_Dhat <- function(a, b, c, d, D, A, Dhat){
	a<-exp(a)
	b<-exp(b)
	y <- (a + (b-a)*exp(-A/c))*(D/Dhat)^d
} # Saturating with diameter # got rid of measurement time, now it's annual rate

#########################################################
##############           Runs          ################## # Gr Mr Rr Diff # Alnus
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 1

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_1_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Alnus_GrMrRrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Gr Diff # Alnus
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 0
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_1_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")

write.csv(looptrack_runs,file="FIA_IBM_Alnus_GrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Mr Diff # Alnus
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 0
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_1_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Alnus_MrDiff_Updated_150729.csv")

#########################################################
##############           Runs          ################## # Gr Mr Diff # Alnus
#########################################################

startage<-1
looptrack_runs<-data.frame(age=startage:250,run1=NA,run2=NA,run3=NA,run4=NA,run5=NA,run6=NA,run7=NA,run8=NA,run9=NA,run10=NA)

for(run_no in 1:4){

g_fixeffect<- 1
m_fixeffect <- 1
r_fixeffect <- 0

# Get initial condition
start_gdata<-gdata[gdata$stdage<=90,] # 163 data points: stick to this for now.

nruns<-1
start_fix<-start_gdata[start_gdata$FIX==1,]
start_non<-start_gdata[start_gdata$FIX==0,]
nflong<-start_fix[!is.na(start_fix$dbh),]$dbh # 26 N fixers
n0long<-start_non[!is.na(start_non$dbh),]$dbh # 175 non fixers

startage<-min(gdata$stdage)
endage<-250
ages<-seq(startage,endage)
Dhat_f<-exp(mean(log(start_fix$dbh),na.rm=TRUE))
Dhat_0<-exp(mean(log(start_non$dbh),na.rm=TRUE))

nfstart<-c(nflong)
n0start<-c(n0long)
rf<-100

start_pba_offset<-1 # start with N fixers being half of what they are. W.L. This is starting with n fixers being 10% larger than what they are.
if(start_pba_offset<1){
	nfstart<-sample(nflong,start_pba_offset*length(nflong)) # slightly more fixers than nonfixers
	n0start<-c(n0long,sample(n0long,length(nfstart)))
} # take a subset of nfs, then replicate the same number of n0s
if(start_pba_offset>1){
	n0start<-sample(n0long,(0.9)*length(n0long)) # W.L.: take random sample from nonfixer(start), with size slightly smaller than the true state.
	nfstart<-c(nflong,sample(nflong,length(n0long)-length(n0start))) # W.L.: keep the fixers data(start), and add 0.1*length(n0long)
} # take a subset of nfs, then replicate the same number of nos. W.L.: shrink nonfixers size by 10% and increase fixer size by 10% of fixer size. The total size is the same as the original starting size.

# create array to store simulated data: 3 dimensions 
nf<-array(NA,dim=c(length(nfstart)*rf,length(ages),nruns))
n0<-array(NA,dim=c(length(n0start)*rf,length(ages),nruns))

ba_f<-array(NA,dim=c(nruns,ncol(nf[,,1]))) # ncol(nf[,,1]) is the number of ages
ba_0<-array(NA,dim=c(nruns,ncol(n0[,,1])))

looptrack<-data.frame(j=1:250,nf_tot=NA,n0_tot=NA,nf_alive=NA,n0_alive=NA,BA_fixer=NA,BA_nonfixer=NA,perc_nf=NA)

k<-1

print(paste("Run number",k))
j<-1
nf[1:length(nfstart),j,k]<-nfstart # Fills in starting age diameter for each run
n0[1:length(n0start),j,k]<-n0start 
	
nf_tot<-nrow(nf[!is.na(nf[,1,k]),,k]) # count how many rows for start age is not NA: number of alive at age class 1 for each run
n0_tot<-nrow(n0[!is.na(n0[,1,k]),,k])
nf_alive<-nf_tot
n0_alive<-n0_tot

BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
perc_nf<-BA_fix/(BA_fix+BA_non)

plot(startage,perc_nf,xlim=c(1,250),ylim=c(0,0.25),xlab="Stand Age",ylab="Proportional Basal Area (%)", main="All Fixers")

	for(t in startage:(endage-1)){
	
		print(paste("Age",t+1))
		print(paste("total N fixers",nf_tot))
		print(paste("total Non fixers",n0_tot))
		print(paste("alive N fixers",nf_alive))
		print(paste("alive Non fixers",n0_alive))
		print(paste("Basal Area N fixers",BA_fix))
		print(paste("Basal Area Non fixers",BA_non))
		print(paste("Percent BA N fixers",BA_fix/(BA_fix+BA_non)))

		looptrack[t,]$nf_tot<-nf_tot
		looptrack[t,]$n0_tot<-n0_tot
		looptrack[t,]$nf_alive<-nf_alive
		looptrack[t,]$n0_alive<-n0_alive
		
		looptrack[t,]$BA_fixer<-BA_fix
		looptrack[t,]$BA_nonfixer<-BA_non
		looptrack[t,]$perc_nf<-BA_fix/(BA_fix+BA_non)
		
		points(t,looptrack[t,]$perc_nf)

		j<-j+1

		if(nf_alive>5000|n0_alive>5000) {
			nfstart_reduced<-sample(nf[,j-1,k][!is.na(nf[,j-1,k])],0.1*length(nf[,j-1,k][!is.na(nf[,j-1,k])]),replace=FALSE)
			n0start_reduced<-sample(n0[,j-1,k][!is.na(n0[,j-1,k])],0.1*length(n0[,j-1,k][!is.na(n0[,j-1,k])]),replace=FALSE)
			nf[,j-1,k]<-NA
			n0[,j-1,k]<-NA
			
			if(length(nfstart_reduced)==0|length(n0start_reduced)==0) break
			
			nf[,j-1,k][1:length(nfstart_reduced)]<-nfstart_reduced
			n0[,j-1,k][1:length(n0start_reduced)]<-n0start_reduced
			}

		# calculate what new size would be if it survives

		gvec_fix<-pv_f_gr_all
		gvec_non<-pv_n_gr_all
		gf<-exp(gMLF_1_f_Dhat(gvec_fix[2],gvec_fix[3],gvec_fix[4],gvec_fix[5],nf[,j-1,k],t+1,Dhat_f))
		if(g_fixeffect==0){
			gf<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],nf[,j-1,k],t+1,Dhat_f))
		}
		g0<-exp(gMLF_4_n_Dhat(gvec_non[2],gvec_non[3],gvec_non[4],gvec_non[5],n0[,j-1,k],t+1,Dhat_0))
		nf[,j,k]<-nf[,j-1,k]*exp(gf) # Next stand age, update diameter
		n0[,j,k]<-n0[,j-1,k]*exp(g0)	
				
		# Find out if each individual survives

		mvec_fix<-pv_f_mr_all
		mvec_non<-pv_n_mr_all
		mf<-fn_ricker_mort_d_t_Dhat(mvec_fix[1],mvec_fix[2],mvec_fix[3],mvec_fix[4],nf[,j-1,k],t+1,Dhat_f) #vector of length n: 0 if dead, 1 if alive
		if(m_fixeffect==0){
			mf<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],nf[,j-1,k],t+1,Dhat_f) # same for non-fixer
			}	
		m0<-fn_ricker_mort_d_t_Dhat(mvec_non[1],mvec_non[2],mvec_non[3],mvec_non[4],n0[,j-1,k],t+1,Dhat_0)
		
		surv_f<-exp(-mf)
		surv_0<-exp(-m0)
		
		surv_count_f<-NULL
		surv_count_0<-NULL

		for(step1 in 1:nf_tot){
				surv_count_f[step1]<-rbinom(1,1,surv_f[step1])	
		}
		for(step2 in 1:n0_tot){
			surv_count_0[step2]<-rbinom(1,1,surv_0[step2])					
		}# if 0 is generated, it means that this tree is dead. 		
		
		nf[which(surv_count_f==0),j,k]<-0
		n0[which(surv_count_0==0),j,k]<-0

		# Now add # of new recruits # Set as the same with Robinia (change it to U.S. trend)

		rvec_fix<-pv_f_rr_all
		rvec_non<-pv_n_rr_all
		rf<-fn4_rec_d_t_Dhat(rvec_fix[1],rvec_fix[2],rvec_fix[3],rvec_fix[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f) 
		if(r_fixeffect==0){
			rf<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(nf[,j-1,k],na.rm=T),t+1,Dhat_f)
		}
		r0<-fn1_rec_d_t_Dhat(rvec_non[1],rvec_non[2],rvec_non[3],rvec_non[4],mean(n0[,j-1,k],na.rm=T),t+1,Dhat_0)
		
		new_fixers<-round(rf*nf_alive)
		new_nonfixers<-round(r0*n0_alive)
		startsize<-12.99546 # calculated this from rdata: FIARRSTAT_recruitment_data_stat_final_140809.csv R==1, mean(dbh0)  -->startsize<-mean(exp(log(rec[!is.na(rec$dbh0),]$dbh0)))
		nf[(nf_tot+1):(nf_tot+new_fixers),j,k]<- startsize
		n0[(n0_tot+1):(n0_tot+new_nonfixers),j,k]<-startsize

		# Update alive counts and total counts
				
		nf[,j,k][nf[,j,k]==0]<-NA
		n0[,j,k][n0[,j,k]==0]<-NA

		nf_tot<-nf_tot+new_fixers
		n0_tot<-n0_tot+new_nonfixers

		nf_alive<-length(nf[,j,k][!is.na(nf[,j,k])&nf[,j,k]!=0])
		n0_alive<-length(n0[,j,k][!is.na(n0[,j,k])&n0[,j,k]!=0])
		if(nf_tot>nrow(nf[,,k])) print("ERROR! NOT ENOUGH FIXER ROWS.")
		if(n0_tot>nrow(n0[,,k])) print("ERROR! NOT ENOUGH NON-FIXER ROWS.")
				
		BA_fix<-sum((nf[,j,k]/2)^2*pi,na.rm=TRUE)
		BA_non<-sum((n0[,j,k]/2)^2*pi,na.rm=TRUE)
		
		if(nf_alive==0|n0_alive==0) break
	}	

looptrack_runs[,run_no+1]<-looptrack$perc_nf

}

# setwd("")
write.csv(looptrack_runs,file="FIA_IBM_Alnus_GrMrDiff_Updated_150729.csv")

#######################################################
#################       Plotting     ##################
#######################################################

# setwd("")
all_grmrrr<-read.csv("FIA_IBM_All_GrMrRrDiff_Updated_150729.csv",header=T)
all_gr<-read.csv("FIA_IBM_All_GrDiff_Updated_150729.csv",header=T)
all_mr<-read.csv("FIA_IBM_All_MrDiff_Updated_150729.csv",header=T)
all_rr<-read.csv("FIA_IBM_All_RrDiff_Updated_150729.csv",header=T)

robinia_grmrrr<-read.csv("FIA_IBM_Robinia_GrMrRrDiff_Updated_150729.csv",header=T)
robinia_gr<-read.csv("FIA_IBM_Robinia_GrDiff_Updated_150729.csv",header=T)
robinia_mr<-read.csv("FIA_IBM_Robinia_MrDiff_Updated_150729.csv",header=T)
robinia_rr<-read.csv("FIA_IBM_Robinia_RrDiff_Updated_150729.csv",header=T)

cerc_grmrrr<-read.csv("FIA_IBM_Cercocarpus_GrMrRrDiff_Updated_150729.csv",header=T)
cerc_gr<-read.csv("FIA_IBM_Cercocarpus_GrDiff_Updated_150729.csv",header=T)
cerc_mr<-read.csv("FIA_IBM_Cercocarpus_MrDiff_Updated_150729.csv",header=T)
cerc_grmr<-read.csv("FIA_IBM_Cercocarpus_GrMrDiff_Updated_150729.csv",header=T)

alnus_grmrrr<-read.csv("FIA_IBM_Alnus_GrMrRrDiff_Updated_150729.csv",header=T)
alnus_gr<-read.csv("FIA_IBM_Alnus_GrDiff_Updated_150729.csv",header=T)
alnus_mr<-read.csv("FIA_IBM_Alnus_MrDiff_Updated_150729.csv",header=T)
alnus_grmr<-read.csv("FIA_IBM_Alnus_GrMrDiff_Updated_150729.csv",header=T)

############################################################################################################
################################ MS MS MS MS MS MS MS MS  PLOT PLOT PLOT PLOT  #############################
################################ PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT PLOT  #############################
############################################################################################################

### All N Fixers

# setwd("")
pdf("Fig18_IBM_All_150930.pdf",width=4.5,height=4.5)

plot(all_grmrrr$run1*100~all_grmrrr$age,xlab="Stand Age (Years)",xlim=c(0,250),ylim=c(0,20),ylab="N Fixer Proportional Basal Area (%)",main="All N Fixers",typ="l",col="black")
points(all_grmrrr$run2*100~all_grmrrr$age,typ="l",col="black")
points(all_grmrrr$run3*100~all_grmrrr$age,typ="l",col="black")
points(all_grmrrr$run4*100~all_grmrrr$age,typ="l",col="black")

points(all_gr$run1*100~all_gr$age,typ="l",col="blue")
points(all_gr$run2*100~all_gr$age,typ="l",col="blue")
points(all_gr$run3*100~all_gr$age,typ="l",col="blue")
points(all_gr$run4*100~all_gr$age,typ="l",col="blue")

points(all_mr$run1*100~all_mr$age,typ="l",col="red")
points(all_mr$run2*100~all_mr$age,typ="l",col="red")
points(all_mr$run3*100~all_mr$age,typ="l",col="red")
points(all_mr$run4*100~all_mr$age,typ="l",col="red")

points(all_rr$run1*100~all_rr$age,typ="l",col="dark green")
points(all_rr$run2*100~all_rr$age,typ="l",col="dark green")
points(all_rr$run3*100~all_rr$age,typ="l",col="dark green")
points(all_rr$run4*100~all_rr$age,typ="l",col="dark green")

legend("topright",col=c("black","blue","red","dark green"),legend=c("All Different","Growth Different","Mortality Different","Recruitment Different"),lty=1,lwd=2,bty="n")

dev.off()

### Robinia

# setwd("")
pdf("Fig19_IBM_Robinia_150930.pdf",width=4.5,height=4.5)

plot(robinia_grmrrr$run1*100~robinia_grmrrr$age,xlab="Stand Age (Years)",xlim=c(0,250),ylim=c(0,20),ylab="N Fixer Proportional Basal Area (%)",main=expression(italic(Robinia~pseudoacacia)),typ="l",col="black")
points(robinia_grmrrr$run2*100~robinia_grmrrr$age,typ="l",col="black")
points(robinia_grmrrr$run3*100~robinia_grmrrr$age,typ="l",col="black")
points(robinia_grmrrr$run4*100~robinia_grmrrr$age,typ="l",col="black")

points(robinia_gr$run1*100~robinia_gr$age,typ="l",col="blue")
points(robinia_gr$run2*100~robinia_gr$age,typ="l",col="blue")
points(robinia_gr$run3*100~robinia_gr$age,typ="l",col="blue")
points(robinia_gr$run4*100~robinia_gr$age,typ="l",col="blue")

points(robinia_mr$run1*100~robinia_mr$age,typ="l",col="red")
points(robinia_mr$run2*100~robinia_mr$age,typ="l",col="red")
points(robinia_mr$run3*100~robinia_mr$age,typ="l",col="red")
points(robinia_mr$run4*100~robinia_mr$age,typ="l",col="red")

points(robinia_rr$run1*100~robinia_rr$age,typ="l",col="dark green")
points(robinia_rr$run2*100~robinia_rr$age,typ="l",col="dark green")
points(robinia_rr$run3*100~robinia_rr$age,typ="l",col="dark green")
points(robinia_rr$run4*100~robinia_rr$age,typ="l",col="dark green")

legend("topright",col=c("black","blue","red","dark green"),legend=c("All Different","Growth Different","Mortality Different","Recruitment Different"),lty=1,lwd=2,bty="n")

dev.off()

### Cercocarpus

# setwd("")
pdf("Fig20_IBM_Cercocarpus_150930.pdf",width=4.5,height=4.5)

plot(cerc_grmrrr$run1*100~cerc_grmrrr$age,xlab="Stand Age (Years)",xlim=c(0,250),ylim=c(0,20),ylab="N Fixer Proportional Basal Area (%)",main=expression(italic(Cercocarpus~ledifolius)),typ="l",col="black")
points(cerc_grmrrr$run2*100~cerc_grmrrr$age,typ="l",col="black")
points(cerc_grmrrr$run3*100~cerc_grmrrr$age,typ="l",col="black")
points(cerc_grmrrr$run4*100~cerc_grmrrr$age,typ="l",col="black")

points(cerc_gr$run1*100~cerc_gr$age,typ="l",col="blue")
points(cerc_gr$run2*100~cerc_gr$age,typ="l",col="blue")
points(cerc_gr$run3*100~cerc_gr$age,typ="l",col="blue")
points(cerc_gr$run4*100~cerc_gr$age,typ="l",col="blue")

points(cerc_mr$run1*100~cerc_mr$age,typ="l",col="red")
points(cerc_mr$run2*100~cerc_mr$age,typ="l",col="red")
points(cerc_mr$run3*100~cerc_mr$age,typ="l",col="red")
points(cerc_mr$run4*100~cerc_mr$age,typ="l",col="red")

points(cerc_grmr$run1*100~cerc_grmr$age,typ="l",col="orange")
points(cerc_grmr$run2*100~cerc_grmr$age,typ="l",col="orange")
points(cerc_grmr$run3*100~cerc_grmr$age,typ="l",col="orange")
points(cerc_grmr$run4*100~cerc_grmr$age,typ="l",col="orange")

legend("topright",col=c("black","blue","red","orange"),legend=c("All Different","Growth Different","Mortality Different","Growth Mortality Different"),lty=1,lwd=2,bty="n")

dev.off()

### Alnus

# setwd("")
pdf("Fig21_IBM_Alnus_150930.pdf",width=4.5,height=4.5)

plot(alnus_grmrrr$run1*100~alnus_grmrrr$age,xlab="Stand Age (Years)",xlim=c(0,250),ylim=c(0,30),ylab="N Fixer Proportional Basal Area (%)",main=expression(italic(Alnus~rubra)),typ="l",col="black")
points(alnus_grmrrr$run2*100~alnus_grmrrr$age,typ="l",col="black")
points(alnus_grmrrr$run3*100~alnus_grmrrr$age,typ="l",col="black")
points(alnus_grmrrr$run4*100~alnus_grmrrr$age,typ="l",col="black")

points(alnus_gr$run1*100~alnus_gr$age,typ="l",col="blue")
points(alnus_gr$run2*100~alnus_gr$age,typ="l",col="blue")
points(alnus_gr$run3*100~alnus_gr$age,typ="l",col="blue")
points(alnus_gr$run4*100~alnus_gr$age,typ="l",col="blue")

points(alnus_mr$run1*100~alnus_mr$age,typ="l",col="red")
points(alnus_mr$run2*100~alnus_mr$age,typ="l",col="red")
points(alnus_mr$run3*100~alnus_mr$age,typ="l",col="red")
points(alnus_mr$run4*100~alnus_mr$age,typ="l",col="red")

points(alnus_grmr$run1*100~alnus_grmr$age,typ="l",col="orange")
points(alnus_grmr$run2*100~alnus_grmr$age,typ="l",col="orange")
points(alnus_grmr$run3*100~alnus_grmr$age,typ="l",col="orange")
points(alnus_grmr$run4*100~alnus_grmr$age,typ="l",col="orange")

legend("topright",col=c("black","blue","red","orange"),legend=c("All Different","Growth Different","Mortality Different","Growth Mortality Different"),lty=1,lwd=2,bty="n")

dev.off()