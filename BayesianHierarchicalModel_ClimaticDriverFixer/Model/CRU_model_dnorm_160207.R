### Climatic Driver of Symbiotic Nitrogen-fixing Tree Success in North America
### Model: dnorm
### Author: WL
 
###############
############### 1=MAT, 2=MAP
###############

model{
	
	# priors
	a1~dunif(-100,100)
	a2~dunif(-100,100)
	int~dunif(0,1) # intercept
	sigma~dunif(0,1)

	MAT.bar <- mean(MAT)
	MAT.sdev <- sd(MAT)
	MAP.bar <- mean(MAP)
	MAP.sdev <- sd(MAP)

	# Likelihood for theta
	for(i in 1:n){
		# m <- int+a1*(MAT[i]-MAT.bar)/MAT.sdev+a2*(MAP[i]-MAP.bar)/MAP.sdev
		y[i]~dnorm(int+a1*(MAT[i]-MAT.bar)/MAT.sdev+a2*(MAP[i]-MAP.bar)/MAP.sdev
,sigma)
	}

	# Prediction
	for (k in 1: length(MAT_x)){
	y_pred_MAT[k] <- int+a1*(MAT_x[k]-MAT.bar)/MAT.sdev
	}

	# Prediction - MAP
	for (k in 1: length(MAP_x)){
	y_pred_MAP[k] <- int+a2*(MAP_x[k]-MAP.bar)/MAP.sdev
	}

}