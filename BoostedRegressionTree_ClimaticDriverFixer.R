### Climatic Driver
### Author: WL
### Boosted Regression Tree Analysis

library(gbm) # 2008 paper, package
library(dismo) # 2016 tutorial, package

###################################
###########  Functions  ###########
###################################

mean.na<-function(data){
  data<-data[!is.na(data)]
  mean(data)
}

matrix_fun<-function(data,var_toExtract){
  data_new<-data[,c("flat","flon",var_toExtract)]
  data_cast<-cast(data_new,flon~flat,mean.na,value=var_toExtract)
  data_matrix<-as.matrix(data_cast[,c(-1)])
  return(data_matrix)
}

matrix_plot<-function(data,main_title,legend_text,zlim){
  image.plot(x_flon+0.5,y_flat+0.5,data,nlevel=20,zlim=zlim,main=main_title,xlab=expression(Longitude~degree),ylab=expression(Latitude~degree~N),cex.lab=1.5,cex.main=2,cex.axis=1.5,legend.args=list(text=legend_text,cex=1),ylim=c(14,80),xlim=c(-170,-60))
  map('world',regions='usa',add=TRUE)
  map('world',regions='mexico',add=TRUE)
  map('world',regions='puerto rico',add=TRUE)
}

min_BA <- 0.005
BA_minadded<-function(min_BA,BA_matrix,BA_matrix_min){
  for(i in 1:nrow(BA_matrix)){
    for(j in 1:ncol(BA_matrix)){
      if(!is.na(BA_matrix[i,j])) BA_matrix_min[i,j]<-max(BA_matrix[i,j],min_BA)
    }
  }
  return(BA_matrix_min*100)
}

ticks<-c(0.5,5,10,50,100)
matrix_plot_log <- function(data, main_title, legend_text, zlim, min_BA){
  data_minadded <- BA_minadded(min_BA, data, data)
  image.plot(x_flon+0.5, y_flat+0.5, log(data_minadded),nlevel=20, zlim=c(log(min(zlim)), log(max(zlim))), axis.args=list(at=log(ticks),labels=ticks), main=main_title, xlab=expression(Longitude~degree),ylab=expression(Latitude~degree~N),cex.lab=1.5,cex.main=2,cex.axis=1.5,legend.args=list(text=legend_text,cex=1),ylim=c(14,80),xlim=c(-170,-60))
  map('world',regions='usa',add=TRUE)
  map('world',regions='mexico',add=TRUE)
  map('world',regions='puerto rico',add=TRUE)
}

meangrid_cal <- function(var_toCal, gridID, unigrid, var_name){
  BA_perc <- tapply(var_toCal,gridID,mean.na)
  BA_perc <- as.data.frame(BA_perc)
  BA_perc$gridID <- rownames(BA_perc)
  BA_perc_latlon <- merge(BA_perc, unigrid, by="gridID")
  colnames(BA_perc_latlon)[2] <- var_name
  return(BA_perc_latlon)
}

poseffect <- function(datamat, clmmat){
  newmat <- matrix(NA, nrow=90, ncol=48)
  for(i in 1:nrow(datamat)){
    for(j in 1:ncol(datamat)){
      if(!is.na(clmmat[i,j]) & clmmat[i,j] > 0)  newmat[i,j] <- datamat[i,j]
    }
  }
  return(newmat)
}

negeffect <- function(datamat, clmmat){
  newmat <- matrix(NA, nrow=90, ncol=48)
  for(i in 1:nrow(datamat)){
    for(j in 1:ncol(datamat)){
      if(!is.na(clmmat[i,j]) & clmmat[i,j] < 0)  newmat[i,j] <- datamat[i,j]
    }
  }
  return(newmat)
}

toNA_fun <- function(matrix, argument){
  newmatrix <- matrix
  if(argument == "less than 0"){
    for(i in 1:nrow(matrix)){
      for(j in 1:ncol(matrix)){
        if(matrix[i,j] <0 & !is.na(matrix[i,j]) ) newmatrix[i, j] <- NA
      }
    }
  }
  if(argument == "greater than 0"){
    for(i in 1:nrow(matrix)){
      for(j in 1:ncol(matrix)){
        if(matrix[i,j] >0 & !is.na(matrix[i,j])) newmatrix[i, j] <- NA
      }
    }
  }
  return(newmatrix)
}

### Load in data
#setwd()
dat<-read.csv("CRU_All_BA_perc_clm10min_FinalwClimate_actvsriz_150312.csv")[,c(-1)]
dat$BA_fix_pos <- floor(dat$BA_fix)
dat$BA_all_pos <- floor(dat$BA_all)

set.seed(1)
train_ind <- sample(seq_len(nrow(dat)), size=round(nrow(dat)/2))
dat_train <- dat[train_ind,]
dat_test <- dat[-train_ind,]

#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
yr2080clm <- clm_tot[,c("map_yr2080","mat_yr2080","flat","flon")]
colnames(yr2080clm)[1:2] <- c("MAP", "MAT")
yr2080clm_MAPonly <- clm_tot[,c("map_yr2080","mat_mean","flat","flon")]
colnames(yr2080clm_MAPonly)[1:2] <- c("MAP", "MAT")
yr2080clm_MATonly <- clm_tot[,c("map_mean","mat_yr2080","flat","flon")]
colnames(yr2080clm_MATonly)[1:2] <- c("MAP", "MAT")

### Climate variable plots
mapmat <- matrix_fun(dat, "MAT")
matrix_plot(mapmat, "MAP","mm",c(-5,30))

########################
######################## # Only 2 variables: MAT and MAP
########################

### All fixers

prfix.tc5.lr01 <- gbm.step(data=dat_train,
                           gbm.x=c(3,7),
                           gbm.y=22, # 22 perc_fix_all, 23 riz, 24 act
                           family="gaussian",
                           tree.complexity=5,
                           learning.rate=0.01,
                           bag.fraction=0.5)
summary(prfix.tc5.lr01)
gbm.plot(prfix.tc5.lr01, write.title = FALSE, n.plots=2, rug=FALSE)
find.int <- gbm.interactions(prfix.tc5.lr01)
find.int$interactions
find.int$rank.list
gbm.perspec(prfix.tc5.lr01, 1, 2)
gbm.plot.fits(prfix.tc5.lr01)

MATpredall <- plot(prfix.tc5.lr01, 'MAT', return.grid=TRUE)
MAPpredall <- plot(prfix.tc5.lr01, 'MAP', return.grid=TRUE)

preds <- predict.gbm(prfix.tc5.lr01,
                     dat_test,
                     n.trees=prfix.tc5.lr01$gbm.call$best.trees,
                     type="response")

dat_preds <- cbind(dat_test$perc_fix_all, preds)
dat$pred_perc_fix <- predict.gbm(prfix.tc5.lr01,
                                 dat,
                                 n.trees=prfix.tc5.lr01$gbm.call$best.trees,
                                 type="response")
yr2080clm$pred_perc_fix <- predict.gbm(prfix.tc5.lr01,
                                       yr2080clm,
                                       n.trees=prfix.tc5.lr01$gbm.call$best.trees,
                                       type="response") # make predictions based on BRT results prfix.tc5.lr01 using MAT and MAP data in yr2080clm
yr2080clm_MAPonly$pred_perc_fix <- predict.gbm(prfix.tc5.lr01,
                                       yr2080clm_MAPonly,
                                       n.trees=prfix.tc5.lr01$gbm.call$best.trees,
                                       type="response") # make predictions based on BRT results prfix.tc5.lr01 using MAT and MAP data in yr2080clm
yr2080clm_MATonly$pred_perc_fix <- predict.gbm(prfix.tc5.lr01,
                                       yr2080clm_MATonly,
                                       n.trees=prfix.tc5.lr01$gbm.call$best.trees,
                                       type="response") # make predictions based on BRT results prfix.tc5.lr01 using MAT and MAP data in yr2080clm


yr2080clm_all <- yr2080clm
yr2080clm_all_MAPonly <- yr2080clm_MAPonly
yr2080clm_all_MATonly <- yr2080clm_MATonly
datall <- dat

# calculate R-square value
1-sum((dat_preds[,1]-dat_preds[,2])^2)/sum((dat_preds[,1]-mean(dat_preds[,1]))^2) # 0.3954

# calculate observed means and plot against predicted means
datall$MATclass <- floor(datall$MAT)
datall$MAPclass <- floor(datall$MAP/10)*10

mean.na <- function(x){
  x <- x[!is.na(x)]
  mean(x)
}

meancal_clm <- function(var, class){
#  var <- eval(parse(text=vartext))
#  class <- eval(parse(text=classtext))
  cal <- tapply(var, class, mean.na)
  cal <- as.data.frame(cal)
  table <- data.frame(var=cal, class=as.numeric(rownames(cal)))
  return(table)
}

MATobsall <- meancal_clm(datall$perc_fix_all, datall$MATclass)
colnames(MATobsall) <- c("perc_fix_all_mean","MATclass")
MAPobsall <- meancal_clm(datall$perc_fix_all, datall$MAPclass)
colnames(MAPobsall) <- c("perc_fix_all_mean","MAPclass")

### Riz fixers
prriz.tc5.lr01 <- gbm.step(data=dat_train,
                           gbm.x=c(3,7),
                           gbm.y=23, # 22 perc_fix_all, 23 riz, 24 act
                           family="gaussian",
                           tree.complexity=5,
                           learning.rate=0.01,
                           bag.fraction=0.5)
summary(prriz.tc5.lr01)
gbm.plot(prriz.tc5.lr01, n.plot=2, write.title=FALSE)
gbm.plot.fits(prriz.tc5.lr01)

MATpredriz <- plot(prriz.tc5.lr01, 'MAT', return.grid=TRUE)
MAPpredriz <- plot(prriz.tc5.lr01, 'MAP', return.grid=TRUE)

preds_riz <- predict.gbm(prriz.tc5.lr01,
                         dat_test,
                         n.trees=prriz.tc5.lr01$gbm.call$best.trees,
                         type="response")

dat_preds_riz <- cbind(dat_test$perc_riz_all, preds_riz)
dat$pred_perc_riz <- predict.gbm(prriz.tc5.lr01,
                                 dat,
                                 n.trees=prriz.tc5.lr01$gbm.call$best.trees,
                                 type="response")
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
yr2080clm <- clm_tot[,c("map_yr2080","mat_yr2080","flat","flon")]
colnames(yr2080clm)[1:2] <- c("MAP", "MAT")
yr2080clm_MAPonly <- clm_tot[,c("map_yr2080","mat_mean","flat","flon")]
colnames(yr2080clm_MAPonly)[1:2] <- c("MAP", "MAT")
yr2080clm_MATonly <- clm_tot[,c("map_mean","mat_yr2080","flat","flon")]
colnames(yr2080clm_MATonly)[1:2] <- c("MAP", "MAT")

yr2080clm$pred_perc_riz <- predict.gbm(prriz.tc5.lr01,
                                       yr2080clm,
                                       n.trees=prriz.tc5.lr01$gbm.call$best.trees,
                                       type="response")
yr2080clm_MAPonly$pred_perc_riz <- predict.gbm(prriz.tc5.lr01,
                                       yr2080clm_MAPonly,
                                       n.trees=prriz.tc5.lr01$gbm.call$best.trees,
                                       type="response")
yr2080clm_MATonly$pred_perc_riz <- predict.gbm(prriz.tc5.lr01,
                                       yr2080clm_MATonly,
                                       n.trees=prriz.tc5.lr01$gbm.call$best.trees,
                                       type="response")

yr2080clm_riz <- yr2080clm
yr2080clm_riz_MAPonly <- yr2080clm_MAPonly
yr2080clm_riz_MATonly <- yr2080clm_MATonly
datriz <- dat

# calculate observed means and plot against predicted means
datriz$MATclass <- floor(datriz$MAT)
datriz$MAPclass <- floor(datriz$MAP/10)*10

# calculate R-square value
1-sum((dat_preds_riz[,1]-dat_preds_riz[,2])^2)/sum((dat_preds_riz[,1]-mean(dat_preds_riz[,1]))^2) # 0.4221

MATobsriz <- meancal_clm(datriz$perc_riz_all, datriz$MATclass)
colnames(MATobsriz) <- c("perc_riz_all_mean","MATclass")
MAPobsriz <- meancal_clm(datriz$perc_riz_all, datriz$MAPclass)
colnames(MAPobsriz) <- c("perc_riz_all_mean","MAPclass")

### Act fixers
# Re-load dat data.frame first
pract.tc5.lr01 <- gbm.step(data=dat_train,
                           gbm.x=c(3,7),
                           gbm.y=24, # 22 perc_fix_all, 23 riz, 24 act
                           family="gaussian",
                           tree.complexity=5,
                           learning.rate=0.01,
                           bag.fraction=0.5)
summary(pract.tc5.lr01)
gbm.plot(pract.tc5.lr01, n.plot=2, write.title=FALSE)
gbm.plot.fits(pract.tc5.lr01)

MATpredact <- plot(pract.tc5.lr01, 'MAT', return.grid=TRUE)
MAPpredact <- plot(pract.tc5.lr01, 'MAP', return.grid=TRUE)

preds_act <- predict.gbm(pract.tc5.lr01,
                         dat_test,
                         n.trees=pract.tc5.lr01$gbm.call$best.trees,
                         type="response")

dat_preds_act <- cbind(dat_test$perc_act_all, preds_act)

dat$pred_perc_act <- predict.gbm(pract.tc5.lr01,
                                 dat,
                                 n.trees=pract.tc5.lr01$gbm.call$best.trees,
                                 type="response")

clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
yr2080clm <- clm_tot[,c("map_yr2080","mat_yr2080","flat","flon")]
colnames(yr2080clm)[1:2] <- c("MAP", "MAT")
yr2080clm_MAPonly <- clm_tot[,c("map_yr2080","mat_mean","flat","flon")]
colnames(yr2080clm_MAPonly)[1:2] <- c("MAP", "MAT")
yr2080clm_MATonly <- clm_tot[,c("map_mean","mat_yr2080","flat","flon")]
colnames(yr2080clm_MATonly)[1:2] <- c("MAP", "MAT")

yr2080clm$pred_perc_act <- predict.gbm(pract.tc5.lr01,
                                       yr2080clm,
                                       n.trees=pract.tc5.lr01$gbm.call$best.trees,
                                       type="response")
yr2080clm_MAPonly$pred_perc_act <- predict.gbm(pract.tc5.lr01,
                                       yr2080clm_MAPonly,
                                       n.trees=pract.tc5.lr01$gbm.call$best.trees,
                                       type="response")
yr2080clm_MATonly$pred_perc_act <- predict.gbm(pract.tc5.lr01,
                                       yr2080clm_MATonly,
                                       n.trees=pract.tc5.lr01$gbm.call$best.trees,
                                       type="response")

yr2080clm_act <- yr2080clm
yr2080clm_act_MAPonly <- yr2080clm_MAPonly
yr2080clm_act_MATonly <- yr2080clm_MATonly
datact <- dat

# calculate observed means and plot against predicted means
datact$MATclass <- floor(datact$MAT)
datact$MAPclass <- floor(datact$MAP/10)*10

# calculate R-square value
1-sum((dat_preds_act[,1]-dat_preds_act[,2])^2)/sum((dat_preds_act[,1]-mean(dat_preds_act[,1]))^2) # 0.05

MATobsact <- meancal_clm(datact$perc_act_all, datact$MATclass)
colnames(MATobsact) <- c("perc_act_all_mean","MATclass")
MAPobsact <- meancal_clm(datact$perc_act_all, datact$MAPclass)
colnames(MAPobsact) <- c("perc_act_all_mean","MAPclass")

#########
######### # Now, plot the three together.
#########

#setwd()
pdf("CRU_160807_predvsobs_mean_threetogether.pdf",width=5.33,height=8)

par(mfrow=c(3,2))
plot(MATpredall$y*100 ~ MATpredall$MAT, xlim=c(-5,30), pch=2, ylim=c(0, 45),col="red", cex=2, xlab=expression(Mean~Annual~Temperature~(~degree~C)), ylab="Proportional Basal Area (%)", main="All N Fixers") # predictions
points(MATobsall$perc_fix_all_mean*100 ~ MATobsall$MATclass, pch=1, col="blue", cex=2) # observations
legend("topleft", bty="n", col=c("red","blue"), pch=c(2,1), legend=c("Predicted Means", "Observed Means"))
plot(MAPpredall$y*100 ~ MAPpredall$MAP, xlim=c(0,410), pch=2,  ylim=c(0, 20), col="red", cex=2, xlab=expression(Mean~Annual~Precipitation~(mm~month^-1)), ylab="Proportional Basal Area (%)") # predictions
points(MAPobsall$perc_fix_all_mean*100 ~ MAPobsall$MAPclass, pch=1, col="blue", cex=2) # observations

plot(MATpredriz$y*100 ~ MATpredriz$MAT, xlim=c(-5,30), pch=2, ylim=c(0,45), col="red", cex=2, xlab=expression(Mean~Annual~Temperature~(~degree~C)), ylab="Proportional Basal Area (%)", main="Rhizobial N Fixers") # predictions
points(MATobsriz$perc_riz_all_mean*100 ~ MATobsriz$MATclass, pch=1, col="blue", cex=2) # observations
# legend("topleft", bty="n", col=c("red","blue"), pch=c(2,1), legend=c("Predicted Means", "Observed Means"))
plot(MAPpredriz$y*100 ~ MAPpredriz$MAP, xlim=c(0,410), pch=2, ylim=c(0,20), col="red", cex=2, xlab=expression(Mean~Annual~Precipitation~(mm~month^-1)), ylab="Proportional Basal Area (%)") # predictions
points(MAPobsriz$perc_riz_all_mean*100 ~ MAPobsriz$MAPclass, pch=1, col="blue", cex=2) # observations

plot(MATpredact$y*100 ~ MATpredact$MAT, xlim=c(-5,30), pch=2, ylim=c(0,2),col="red", cex=2, xlab=expression(Mean~Annual~Temperature~(~degree~C)), ylab="Proportional Basal Area (%)", main="Actinorhizal N Fixers") # predictions
points(MATobsact$perc_act_all_mean*100 ~ MATobsact$MATclass, pch=1, col="blue", cex=2) # observations
# legend("topleft", bty="n", col=c("red","blue"), pch=c(2,1), legend=c("Predicted Means", "Observed Means"))
plot(MAPpredact$y*100 ~ MAPpredact$MAP, xlim=c(0,410), pch=2, ylim=c(0,5),col="red", cex=2, xlab=expression(Mean~Annual~Precipitation~(mm~month^-1)), ylab="Proportional Basal Area (%)") # predictions
points(MAPobsact$perc_act_all_mean*100 ~ MAPobsact$MAPclass, pch=1, col="blue", cex=2) # observations

dev.off()

#########
######### # Plotting all
#########

library(rjags)
library(R2jags)
library(fields) # for image.plot()
library(reshape) # for cast()
library(maps)
library(boot)
data(mexicoMapEnv)


dat <- datall[,!(names(datall) %in% c("gridID"))]
ftdat <- yr2080clm_all

ftdat_map <- yr2080clm_all_MAPonly
ftdat_mat <- yr2080clm_all_MATonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]

### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)

### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat", "flon"), all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat", "flon"), all=TRUE)

# all fixers
perc_fix <- meangrid_cal(dat$perc_fix_all, dat$gridID, unigrid, "perc_fix")
pred_fix <- meangrid_cal(dat$pred_perc_fix, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_fix, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_fix, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_fix, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"

clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_all <- obs_ba
crt_ba_all <- crt_ba
ft_ba_all <- ft_ba
crt_obs_ba_all <- crt_obs_ba
ft_crt_ba_all <- ft_crt_ba

#setwd()
pdf("CRU_160821_MATMAPOnly_vs_Current.pdf", width=7, height=3.5)
par(mfrow=c(1,2))
plot(x_2080mat ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="yr2080 MAT + Current MAP", xlim=c(-0.5, 1), ylim=c(-0.5, 1))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
plot(x_2080map ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="Current MAT + yr2080 MAP", xlim=c(-0.5, 1), ylim=c(-0.5, 1))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
dev.off()

map_ba <- matrix_fun(clm_tot, "map_mean")
mat_ba <- matrix_fun(clm_tot, "mat_mean")
map_ba_2080 <- matrix_fun(clm_tot, "map_yr2080")
mat_ba_2080 <- matrix_fun(clm_tot, "mat_yr2080")

matrix_plot_log(crt_obs_ba_all, "Over-prediction","BA%",c(0.1,100),0.001) # here, remember to not multiply by 100. The multiplication messes up with the answer.
# matrix_plot(crt_obs_ba_all*100, "Over-prediction", "BA%", c(0,100))

#setwd()
pdf("CRU_160907_ClimateVariables.pdf",height=7,width=7)
par(mfrow=c(2,2))
# matrix_plot(map_ba,"Current Precipitation",expression(mm~month^-1),c(0,351))
matrix_plot(mat_ba,"Current Temperature",expression(~degree~C),c(-4,31))
matrix_plot_log(map_ba/100, "Current Precipitation", expression(mm~month^-1),c(4,421))
# matrix_plot(map_ba_2080,"yr2080 Precipitation",expression(mm~month^-1),c(0,421))
matrix_plot(mat_ba_2080,"yr2080 Temperature",expression(~degree~C),c(-4,31))
matrix_plot_log(map_ba_2080/100, "yr2080 Precipitation",expression(mm~month^-1),c(4,421))
dev.off()

change_mat <- mat_ba_2080 - mat_ba
change_map <- map_ba_2080 - map_ba

#setwd()
pdf("CRU_160907_ClimateVariables_change_SI.pdf",height=7,width=7)
par(mfrow=c(2,2))
matrix_plot(change_mat, "Increase in Temperature", expression(~degree~C), c(0, 11))
matrix_plot(-change_mat, "Decrease in Temperature", expression(~degree~C), c(0, 4))
matrix_plot_log(change_map/100, "Increase in Precipitation", expression(mm~month^-1), c(0.5,290))
matrix_plot_log(-change_map/100, "Decrease in Precipitation", expression(mm~month^-1), c(0.5,80))
# matrix_plot(change_map, "Increase in Precipitation", expression(mm~month^-1), c(0,290))
# matrix_plot(-change_map, "Decrease in Precipitation", expression(mm~month^-1), c(0,80))
dev.off()

# effects into categories
all_warming <- poseffect(x_2080mat, mat_2080_crt)
all_cooling <- negeffect(x_2080mat, mat_2080_crt)
all_wetting <- poseffect(x_2080map, map_2080_crt)
all_drying <- negeffect(x_2080map, map_2080_crt)

par(mfrow=c(2,2))
matrix_plot(all_warming*100,"All - Warming Effect","BA%",c(0,100))
matrix_plot(all_cooling*100,"Cooling Effect","BA%",c(0,2))
matrix_plot(all_wetting*100, "Wetting Effect", "BA%", c(0,100))
matrix_plot(all_drying*100, "Drying Effect", "BA%", c(-5,100))

#########
######### # Plotting riz
#########

dat <- datriz[,!(names(datriz) %in% c("gridID"))]
ftdat <- yr2080clm_riz
ftdat_mat <- yr2080clm_riz_MATonly
ftdat_map <- yr2080clm_riz_MAPonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]

### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)

### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat","flon"),all=TRUE)

# rhizobial fixers
perc_fix <- meangrid_cal(dat$perc_riz_all, dat$gridID, unigrid, "perc_fix") # here, extract rhizobial %
pred_fix <- meangrid_cal(dat$pred_perc_riz, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_riz, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_riz, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_riz, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"

clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_riz <- obs_ba
crt_ba_riz <- crt_ba
ft_ba_riz <- ft_ba
crt_obs_ba_riz <- crt_obs_ba
ft_crt_ba_riz <- ft_crt_ba

#setwd()
pdf("CRU_160821_MATMAPOnly_vs_Current_Rhizobial.pdf", width=7, height=3.5)
par(mfrow=c(1,2))
plot(x_2080mat ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="yr2080 MAT + Current MAP", xlim=c(-0.5, 1), ylim=c(-0.5, 1))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
plot(x_2080map ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="Current MAT + yr2080 MAP", xlim=c(-0.5, 1), ylim=c(-0.5, 1))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
dev.off()

# effects into categories
riz_warming <- poseffect(ft_ba_mat, mat_2080_crt)
riz_cooling <- negeffect(ft_ba_mat, mat_2080_crt)
riz_wetting <- poseffect(ft_ba_map, map_2080_crt)
riz_drying <- negeffect(ft_ba_map, map_2080_crt)

par(mfrow=c(2,2))
matrix_plot(riz_warming*100,"Riz - Warming Effect","BA%",c(0,100))
matrix_plot(riz_cooling*100,"Cooling Effect","BA%",c(0,2))
matrix_plot(riz_wetting*100, "Wetting Effect", "BA%", c(0,100))
matrix_plot(riz_drying*100, "Drying Effect", "BA%", c(-5,100))


#########
######### # Plotting act
#########

dat <- datact[,!(names(datact) %in% c("gridID"))]
ftdat <- yr2080clm_act
ftdat_map <- yr2080clm_act_MAPonly
ftdat_mat <- yr2080clm_act_MATonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]

### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)

### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat","flon"),all=TRUE)

# actinorhizal fixers
perc_fix <- meangrid_cal(dat$perc_act_all, dat$gridID, unigrid, "perc_fix") # here, extract actinorhizal %
pred_fix <- meangrid_cal(dat$pred_perc_act, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_act, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_act, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_act, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"

clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_act <- obs_ba
crt_ba_act <- crt_ba
ft_ba_act <- ft_ba
crt_obs_ba_act <- crt_obs_ba
ft_crt_ba_act <- ft_crt_ba

#setwd()
pdf("CRU_160821_MATMAPOnly_vs_Current_Actinorhizal.pdf", width=7, height=3.5)
par(mfrow=c(1,2))
plot(x_2080mat ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="yr2080 MAT + Current MAP", xlim=c(-0.1, 0.15), ylim=c(-0.1, 0.15))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
plot(x_2080map ~ x_2080both, xlab="yr2080 MAT + MAP", ylab="Current MAT + yr2080 MAP", xlim=c(-0.1, 0.15), ylim=c(-0.1, 0.15))
abline(0, 1, col="red", lwd=2)
abline(0,0, col="black", lty=2)
abline(v=0, lty=2)
dev.off()

# effects into categories
act_warming <- poseffect(ft_ba_mat, mat_2080_crt)
act_cooling <- negeffect(ft_ba_mat, mat_2080_crt)
act_wetting <- poseffect(ft_ba_map, map_2080_crt)
act_drying <- negeffect(ft_ba_map, map_2080_crt)

par(mfrow=c(2,2))
matrix_plot(act_warming*100,"Act - Warming Effect","BA%",c(0,11))
matrix_plot(act_cooling*100,"Cooling Effect","BA%",c(0,1))
matrix_plot(act_wetting*100, "Wetting Effect", "BA%", c(0,16))
matrix_plot(act_drying*100, "Drying Effect", "BA%", c(0,10))

#########
######### # Plot observations vs. current predictions
#########

#setwd()
pdf("CRU_160907_Comparisons_obs_crt_yr2080.pdf",height=10.5,width=14)

par(mfrow=c(3,4))

# all
matrix_plot_log(obs_ba_all, "Observation","BA%",c(0.5,100))
matrix_plot_log(crt_ba_all, "Current Prediction","BA%",c(0.5,100))
# matrix_plot_log(ft_ba, "yr2080 Prediction","BA%",c(0.1,100))
range(crt_obs_ba_all*100,na.rm=T)
matrix_plot_log(toNA_fun(crt_obs_ba_all, "less than 0"), "Over-prediction","BA%",c(0.5,100)) # here, remember to not multiply by 100. The multiplication messes up with the answer.
# matrix_plot(crt_obs_ba_all*100, "Over-prediction", "BA%", c(0,100))
matrix_plot_log(-toNA_fun(crt_obs_ba_all, "greater than 0"), "Under-prediction", "BA%", c(0.5,100))
# matrix_plot(-crt_obs_ba_all*100, "Under-prediction", "BA%", c(0,100))
# matrix_plot(ft_crt_ba*100, "Future Increase", "BA%", c(0,100))
# matrix_plot(-ft_crt_ba*100, "Future Decrease", "BA%", c(0,100))

# riz
matrix_plot_log(obs_ba_riz, " ","BA%",c(0.5,100))
matrix_plot_log(crt_ba_riz, " ","BA%",c(0.5,100))
# matrix_plot_log(ft_ba, " ","BA%",c(0.1,100))
matrix_plot_log(toNA_fun(crt_obs_ba_riz,"less than 0"), " ", "BA%", c(0.5,100))
matrix_plot_log(-toNA_fun(crt_obs_ba_riz,"greater than 0"), " ", "BA%", c(0.5,100))
# matrix_plot(crt_obs_ba_riz*100, " ", "BA%", c(0,100))
# matrix_plot(-crt_obs_ba_riz*100, " ", "BA%", c(0,100))
# matrix_plot(ft_crt_ba*100, "Future Increase", "BA%", c(0,100))
# matrix_plot(-ft_crt_ba*100, "Future Decrease", "BA%", c(0,100))

# act
matrix_plot_log(obs_ba_act, " ","BA%",c(0.5,100))
matrix_plot_log(crt_ba_act, " ","BA%",c(0.5,100))
# matrix_plot_log(ft_ba, " ","BA%",c(0.1,100))
matrix_plot_log(toNA_fun(crt_obs_ba_act,"less than 0"), " ", "BA%", c(0.5,100))
matrix_plot_log(-toNA_fun(crt_obs_ba_act, "greater than 0"), " ", "BA%", c(0.5,100))
# matrix_plot(crt_obs_ba_act*100, " ", "BA%", c(0,10))
# matrix_plot(-crt_obs_ba_act*100, " ", "BA%", c(0,40))
# matrix_plot(ft_crt_ba*100, "Future Increase", "BA%", c(0,11))
# matrix_plot(-ft_crt_ba*100, "Future Decrease", "BA%", c(0,11))

dev.off()

########
######## # Plotting future and inc. vs. dec.
########

#setwd()
pdf("CRU_160907_Comparisons_ft_inc_dec_MATonly.pdf",height=10.5,width=10.5)

par(mfrow=c(3,3))

# all
# matrix_plot_log(obs_ba_all, "Observation","BA%",c(0.1,100))
# matrix_plot_log(crt_ba_all, "Current Prediction","BA%",c(0.1,100))
matrix_plot_log(ft_ba_all, "yr2080 Prediction","BA%",c(0.5,100))
# matrix_plot(crt_obs_ba_all*100, "Over-prediction", "BA%", c(0,100))
# matrix_plot(-crt_obs_ba_all*100, "Under-prediction", "BA%", c(0,100))
matrix_plot_log(toNA_fun(ft_crt_ba_all, "less than 0"), "Future Increase", "BA%", c(0.5,100))
matrix_plot_log(-toNA_fun(ft_crt_ba_all,"greater than 0"), "Future Decrease", "BA%", c(0.5,100))
# matrix_plot(ft_crt_ba_all*100, "Future Increase", "BA%", c(0,100))
# matrix_plot(-ft_crt_ba_all*100, "Future Decrease", "BA%", c(0,100))

# riz
# matrix_plot_log(obs_ba_riz, " ","BA%",c(0.1,100))
# matrix_plot_log(crt_ba_riz, " ","BA%",c(0.1,100))
matrix_plot_log(ft_ba_riz, " ","BA%",c(0.5,100))
# matrix_plot(crt_obs_ba_riz*100, " ", "BA%", c(0,100))
# matrix_plot(-crt_obs_ba_riz*100, " ", "BA%", c(0,100))
matrix_plot_log(toNA_fun(ft_crt_ba_riz,"less than 0"), " ", "BA%", c(0.5,100))
matrix_plot_log(-toNA_fun(ft_crt_ba_riz,"greater than 0"), " ", "BA%", c(0.5,100))
# matrix_plot(ft_crt_ba_riz*100, " ", "BA%", c(0,100))
# matrix_plot(-ft_crt_ba_riz*100, " ", "BA%", c(0,100))

# act
# matrix_plot_log(obs_ba_act, " ","BA%",c(0.1,100))
# matrix_plot_log(crt_ba_act, " ","BA%",c(0.1,100))
matrix_plot_log(ft_ba_act, " ","BA%",c(0.5,100))
# matrix_plot(crt_obs_ba_act*100, " ", "BA%", c(0,10))
# matrix_plot(-crt_obs_ba_act*100, " ", "BA%", c(0,40))
matrix_plot_log(toNA_fun(ft_crt_ba_act,"less than 0"), " ", "BA%", c(0.5,11))
matrix_plot_log(-toNA_fun(ft_crt_ba_act,"greater than 0"), " ", "BA%", c(0.5,11))
# matrix_plot(ft_crt_ba_act*100, " ", "BA%", c(0,11))
# matrix_plot(-ft_crt_ba_act*100, " ", "BA%", c(0,11))

dev.off()

### Difference in Temperature and Precipitation

# map_ba <- matrix_fun(clm_tot, "map_mean")
# mat_ba <- matrix_fun(clm_tot, "mat_mean")
# map_ba_2080 <- matrix_fun(clm_tot, "map_yr2080")
# mat_ba_2080 <- matrix_fun(clm_tot, "mat_yr2080")

#setwd()
pdf("CRU_160807_Change_ft_inc_dec.pdf",height=7,width=7)

par(mfrow=c(2,2))
mat_2080_crt <- mat_ba_2080 - mat_ba
map_2080_crt <- map_ba_2080 - map_ba
matrix_plot(mat_2080_crt, "Increase in Temperature", expression(~degree~C), c(0, 11))
matrix_plot(-mat_2080_crt, "Decrease in Temperature", expression(~degree~C), c(0, 4))
matrix_plot(map_2080_crt, "Increase in Precipitation", "mm/month", c(0, 290))
matrix_plot(-map_2080_crt, "Decrease in Precipitation", "mm/month", c(0, 80))

dev.off()

############################################################
############################################################ # Plot: effects into categories
############################################################

#########
######### # Plotting all
#########

dat <- datall[,!(names(datall) %in% c("gridID"))]
ftdat <- yr2080clm_all
ftdat_map <- yr2080clm_all_MAPonly
ftdat_mat <- yr2080clm_all_MATonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)
### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat", "flon"), all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat", "flon"), all=TRUE)

# all fixers
perc_fix <- meangrid_cal(dat$perc_fix_all, dat$gridID, unigrid, "perc_fix")
pred_fix <- meangrid_cal(dat$pred_perc_fix, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_fix, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_fix, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_fix, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"
clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_all <- obs_ba
crt_ba_all <- crt_ba
ft_ba_all <- ft_ba
crt_obs_ba_all <- crt_obs_ba
ft_crt_ba_all <- ft_crt_ba

# effects into categories
all_warming <- poseffect(x_2080mat, mat_2080_crt)
all_cooling <- negeffect(x_2080mat, mat_2080_crt)
all_wetting <- poseffect(x_2080map, map_2080_crt)
all_drying <- negeffect(x_2080map, map_2080_crt)
all_comb <- x_2080both - (x_2080mat + x_2080map)

all_warming_inc <- poseffect(all_warming, ft_crt_ba_all)
all_cooling_inc <- poseffect(all_cooling, ft_crt_ba_all)
all_drying_inc <- poseffect(all_drying, ft_crt_ba_all)
all_wetting_inc <- poseffect(all_wetting, ft_crt_ba_all)
all_comb_inc <- poseffect(all_comb, ft_crt_ba_all)

all_warming_dec <- negeffect(all_warming, ft_crt_ba_all)
all_cooling_dec <- negeffect(all_cooling, ft_crt_ba_all)
all_drying_dec <- negeffect(all_drying, ft_crt_ba_all)
all_wetting_dec <- negeffect(all_wetting, ft_crt_ba_all)
all_comb_dec <- negeffect(all_comb, ft_crt_ba_all)

#setwd()
pdf("CRU_160902_all_yr2080_effects.pdf",width=7,height=7)
par(mfrow=c(2,2))
matrix_plot(all_warming*100,"All - Warming Effect","BA%",c(-80,80))
matrix_plot(all_cooling*100,"Cooling Effect","BA%",c(-0.1,0.1))
matrix_plot(all_wetting*100, "Wetting Effect", "BA%", c(-40,40))
matrix_plot(all_drying*100, "Drying Effect", "BA%", c(-40,40))
dev.off()

#setwd()
pdf("CRU_160902_all_inc_yr2080_effects_categories.pdf",width=15,height=9)
par(mfrow=c(3,5))
# increase
matrix_plot_log(ft_crt_ba_all, "Future Increase", "BA%", c(0.5,100))
plot.new()
plot.new()
plot.new()
plot.new()
###############################
############################### 
###############################
ratio <- 3.5/4
width <- 6
height <- width*ratio

#setwd()
pdf("2014CRU_SI_Fig3A_Oct2016.pdf",width=width, height=height)
matrix_plot_log(all_warming_inc, "Warming Effect", "BA%", c(0.5,80),0.005)
matrix_plot_log(all_wetting_inc, "Wetting Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(all_drying_inc, "Drying Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(all_comb_inc, "Combined Effect", "BA%", c(0.5, 45),0.005)
matrix_plot_log(-all_warming_inc, " ", "BA%", c(0.5, 30),0.005)
matrix_plot_log(-all_wetting_inc, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-all_drying_inc, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-all_comb_inc, " ", "BA%", c(0.5, 40),0.005)
dev.off()

#setwd()
pdf("2014CRU_SI_Fig3B_Oct2016.pdf",width=width, height=height)
matrix_plot_log(all_warming_dec, "Warming Effect", "BA%", c(0.5,35),0.005)
matrix_plot_log(all_cooling_dec, "Cooling Effect", "BA%", c(0.5,0.6),0.005)
matrix_plot_log(all_wetting_dec, "Wetting Effect", "BA%", c(0.5, 16),0.005)
matrix_plot_log(all_drying_dec, "Drying Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(all_comb_dec, "Combined Effect", "BA%", c(0.5, 35),0.005)
matrix_plot_log(-all_warming_dec, " ", "BA%", c(0.5, 55),0.005)
matrix_plot_log(-all_cooling_dec, " ", "BA%", c(0.5, 0.6),0.005)
matrix_plot_log(-all_wetting_dec, " ", "BA%", c(0.5, 30),0.005)
matrix_plot_log(-all_drying_dec, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-all_comb_dec, " ", "BA%", c(0.5, 30),0.005)
dev.off()

#########
######### # Plotting riz
#########

dat <- datriz[,!(names(datriz) %in% c("gridID"))]
ftdat <- yr2080clm_riz
ftdat_mat <- yr2080clm_riz_MATonly
ftdat_map <- yr2080clm_riz_MAPonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)
### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat", "flon"), all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat", "flon"), all=TRUE)

# rhizobial fixers
perc_fix <- meangrid_cal(dat$perc_riz_all, dat$gridID, unigrid, "perc_fix") # here, extract rhizobial %
pred_fix <- meangrid_cal(dat$pred_perc_riz, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_riz, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_riz, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_riz, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"
clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_all <- obs_ba
crt_ba_all <- crt_ba
ft_ba_all <- ft_ba
crt_obs_ba_all <- crt_obs_ba
ft_crt_ba_all <- ft_crt_ba

# effects into categories
riz_warming <- poseffect(x_2080mat, mat_2080_crt)
riz_cooling <- negeffect(x_2080mat, mat_2080_crt)
riz_wetting <- poseffect(x_2080map, map_2080_crt)
riz_drying <- negeffect(x_2080map, map_2080_crt)
riz_comb <-x_2080both - (x_2080mat + x_2080map)

riz_warming_inc <- poseffect(riz_warming, ft_crt_ba_all)
riz_cooling_inc <- poseffect(riz_cooling, ft_crt_ba_all)
riz_drying_inc <- poseffect(riz_drying, ft_crt_ba_all)
riz_wetting_inc <- poseffect(riz_wetting, ft_crt_ba_all)
riz_comb_inc <- poseffect(riz_comb, ft_crt_ba_all)

riz_warming_dec <- negeffect(riz_warming, ft_crt_ba_all)
riz_cooling_dec <- negeffect(riz_cooling, ft_crt_ba_all)
riz_drying_dec <- negeffect(riz_drying, ft_crt_ba_all)
riz_wetting_dec <- negeffect(riz_wetting, ft_crt_ba_all)
riz_comb_dec <- negeffect(riz_comb, ft_crt_ba_all)

#setwd()
pdf("CRU_160902_riz_yr2080_effects.pdf",width=7,height=7)
par(mfrow=c(2,2))
matrix_plot(riz_warming*100,"Riz - Warming Effect","BA%",c(-80,80))
matrix_plot(riz_cooling*100,"Cooling Effect","BA%",c(-0.1,0.1))
matrix_plot(riz_wetting*100, "Wetting Effect", "BA%", c(-40,40))
matrix_plot(riz_drying*100, "Drying Effect", "BA%", c(-40,40))
dev.off()

#setwd()
pdf("2014CRU_SI_Fig3C_Oct2016.pdf",width=width, height=height)
matrix_plot_log(riz_warming_inc, "Warming Effect", "BA%", c(0.5,80),0.005)
matrix_plot_log(riz_cooling_inc, "Cooling Effect", "BA%", c(0.5,0.6),0.005)
matrix_plot_log(riz_wetting_inc, "Wetting Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(riz_drying_inc, "Drying Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(riz_comb_inc, "Combined Effect", "BA%", c(0.5, 50),0.005)
matrix_plot_log(-riz_warming_inc, " ", "BA%", c(0.5, 30),0.005)
matrix_plot_log(-riz_cooling_inc, " ", "BA%", c(0.5, 30),0.005)
matrix_plot_log(-riz_wetting_inc, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-riz_drying_inc, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-riz_comb_inc, " ", "BA%", c(0.5, 41),0.005)
dev.off()

#setwd()
pdf("2014CRU_SI_Fig3D_Oct2016.pdf",width=width, height=height)
matrix_plot_log(riz_warming_dec, "Warming Effect", "BA%", c(0.5,35),0.005)
plot.new()
matrix_plot_log(riz_wetting_dec, "Wetting Effect", "BA%", c(0.5, 16),0.005)
matrix_plot_log(riz_drying_dec, "Drying Effect", "BA%", c(0.5, 40),0.005)
matrix_plot_log(riz_comb_dec, "Combined Effect", "BA%", c(0.5, 35),0.005)
matrix_plot_log(-riz_warming_dec, " ", "BA%", c(0.5, 55),0.005)
plot.new()
matrix_plot_log(-riz_wetting_dec, " ", "BA%", c(0.5, 30),0.005)
matrix_plot_log(-riz_drying_dec, " ", "BA%", c(0.5, 40),0.005)
matrix_plot_log(-riz_comb_dec, " ", "BA%", c(0.5, 35),0.005)
dev.off()

#########
######### # Plotting act
#########

dat <- datact[,!(names(datact) %in% c("gridID"))]
ftdat <- yr2080clm_act
ftdat_map <- yr2080clm_act_MAPonly
ftdat_mat <- yr2080clm_act_MATonly

### Load in climate data
#setwd()
clm_tot<-read.csv("CRU_climate_current_mapCorrected_yr2080_160126.csv",header=T)[,c(-1)]
### Create matrix for plotting
x_flon<-seq(-154,-65,1)
y_flat<-seq(14,61,1)
### Create a grid cell matrix
xy_flatlon<-data.frame(flat=rep(y_flat,each=length(x_flon)),flon=rep(x_flon,length(y_flat)))
clm_tot<-merge(xy_flatlon,clm_tot,by=c("flat","flon"),all=TRUE)

### Calculate per grid means
unigrid <- unique(dat[,c("flat","flon")])
unigrid$gridID <- 1:nrow(unigrid)
dat <- merge(dat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat <- merge(ftdat, unigrid, by=c("flat","flon"),all=TRUE)
ftdat_mat <- merge(ftdat_mat, unigrid, by=c("flat", "flon"), all=TRUE)
ftdat_map <- merge(ftdat_map, unigrid, by=c("flat", "flon"), all=TRUE)

# actinorhizal fixers
perc_fix <- meangrid_cal(dat$perc_act_all, dat$gridID, unigrid, "perc_fix") # here, extract actinorhizal %
pred_fix <- meangrid_cal(dat$pred_perc_act, dat$gridID, unigrid, "pred_fix")
pred_fix_ft <- meangrid_cal(ftdat$pred_perc_act, ftdat$gridID, unigrid, "pred_fix_ft")
pred_fix_ft_mat <- meangrid_cal(ftdat_mat$pred_perc_act, ftdat_mat$gridID, unigrid, "pred_fix_ft_mat")
colnames(pred_fix_ft_mat)[2] <- "pred_fix_ft_mat"
pred_fix_ft_map <- meangrid_cal(ftdat_map$pred_perc_act, ftdat_map$gridID, unigrid, "pred_fix_ft_map")
colnames(pred_fix_ft_map)[2] <- "pred_fix_ft_map"

clm_tot <- merge(clm_tot, perc_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_mat[,-1],by=c("flat","flon"),all=TRUE)
clm_tot <- merge(clm_tot, pred_fix_ft_map[,-1],by=c("flat","flon"),all=TRUE)
obs_ba <- matrix_fun(clm_tot,"perc_fix")
crt_ba <- matrix_fun(clm_tot,"pred_fix")
ft_ba <- matrix_fun(clm_tot, "pred_fix_ft")
ft_ba_mat <- matrix_fun(clm_tot, "pred_fix_ft_mat")
ft_ba_map <- matrix_fun(clm_tot, "pred_fix_ft_map")
crt_obs_ba <- crt_ba - obs_ba
ft_crt_ba <- ft_ba - crt_ba
x_2080mat <- ft_ba_mat - crt_ba
x_2080map <- ft_ba_map - crt_ba
x_2080both <- ft_ba - crt_ba

obs_ba_all <- obs_ba
crt_ba_all <- crt_ba
ft_ba_all <- ft_ba
crt_obs_ba_all <- crt_obs_ba
ft_crt_ba_all <- ft_crt_ba

# effects into categories
act_warming <- poseffect(x_2080mat, mat_2080_crt)
act_cooling <- negeffect(x_2080mat, mat_2080_crt)
act_wetting <- poseffect(x_2080map, map_2080_crt)
act_drying <- negeffect(x_2080map, map_2080_crt)
act_comb <- x_2080both - (x_2080mat + x_2080map)

act_warming_inc <- poseffect(act_warming, ft_crt_ba_all)
act_cooling_inc <- poseffect(act_cooling, ft_crt_ba_all)
act_drying_inc <- poseffect(act_drying, ft_crt_ba_all)
act_wetting_inc <- poseffect(act_wetting, ft_crt_ba_all)
act_comb_inc <- poseffect(act_comb, ft_crt_ba_all)

act_warming_dec <- negeffect(act_warming, ft_crt_ba_all)
act_cooling_dec <- negeffect(act_cooling, ft_crt_ba_all)
act_drying_dec <- negeffect(act_drying, ft_crt_ba_all)
act_wetting_dec <- negeffect(act_wetting, ft_crt_ba_all)
act_comb_dec <- negeffect(act_comb, ft_crt_ba_all)

#setwd()
pdf("CRU_160902_act_yr2080_effects.pdf",width=7,height=7)
par(mfrow=c(2,2))
matrix_plot(act_warming*100,"Act - Warming Effect","BA%",c(-10,10))
matrix_plot(act_cooling*100,"Cooling Effect","BA%",c(-0.05,0.05))
matrix_plot(act_wetting*100, "Wetting Effect", "BA%", c(-15,15))
matrix_plot(act_drying*100, "Drying Effect", "BA%", c(-5,5))
dev.off()

#setwd()
pdf("2014CRU_SI_Fig3E_Oct2016.pdf",width=width, height=height)
matrix_plot_log(act_warming_inc, "Warming Effect", "BA%", c(0.5,11),0.005)
plot.new()
# matrix_plot_log(act_cooling_inc, "Cooling Effect", "BA%", c(0.5,0.6))
matrix_plot_log(act_wetting_inc, "Wetting Effect", "BA%", c(0.5, 15),0.005)
matrix_plot_log(act_drying_inc, "Drying Effect", "BA%", c(0.5, 0.6),0.005)
matrix_plot_log(act_comb_inc, "Combined Effect", "BA%", c(0.5, 10),0.005)
matrix_plot_log(-act_warming_inc, " ", "BA%", c(0.5, 0.6),0.005)
plot.new()
# matrix_plot_log(-act_cooling_inc, " ", "BA%", c(0.5, 30))
matrix_plot_log(-act_wetting_inc, " ", "BA%", c(0.5, 0.6),0.005)
matrix_plot_log(-act_drying_inc, " ", "BA%", c(0.5, 0.6),0.005)
matrix_plot_log(-act_comb_inc, " ", "BA%", c(0.5, 12),0.005)
dev.off()

#setwd()
pdf("2014CRU_SI_Fig3F_Oct2016.pdf",width=width, height=height)
matrix_plot_log(act_warming_dec, "Warming Effect", "BA%", c(0.5,0.6),0.005)
matrix_plot_log(act_cooling_dec, "Cooling Effect", "BA%", c(0.5,0.6),0.005)
matrix_plot_log(act_wetting_dec, "Wetting Effect", "BA%", c(0.5, 5),0.005)
matrix_plot_log(act_drying_dec, "Drying Effect", "BA%", c(0.5, 5),0.005)
matrix_plot_log(act_comb_dec, "Combined Effect", "BA%", c(0.5, 5),0.005)
matrix_plot_log(-act_warming_dec, " ", "BA%", c(0.5, 10),0.005)
matrix_plot_log(-act_cooling_dec, " ", "BA%", c(0.5, 0.6),0.005)
matrix_plot_log(-act_wetting_dec, " ", "BA%", c(0.5, 5),0.005)
matrix_plot_log(-act_drying_dec, " ", "BA%", c(0.5, 5),0.005)
matrix_plot_log(-act_comb_dec, " ", "BA%", c(0.5, 6),0.005)
dev.off()