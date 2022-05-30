### Climatic Driver of Symbiotic Nitrogen-fixing Tree Success in North America
### Model: dnorm
### Author: WL
 
###############
############### 1=MAT, 2=MAP
###############

model{
	
	# priors
	a1~dunif(-100,100)
	b1~dunif(-100,100)
	a2~dunif(-100,100)
	b2~dunif(-100,100)
	c2~dunif(-100,100)
	int~dunif(0,1) # intercept
	sigma~dunif(0,1)

	MAT.bar <- mean(MAT)
	MAT.sdev <- sd(MAT)
	MAP.bar <- mean(MAP)
	MAP.sdev <- sd(MAP)

	# Likelihood for theta
	for(i in 1:n){
		# m <- int+a1*(MAT[i]-MAT.bar)/MAT.sdev+a2*(MAP[i]-MAP.bar)/MAP.sdev
		# mat_m <- exp(a1+ b1*(MAT[i]-MAT.bar)/MAT.sdev) /(1+ exp(a1+b1*(MAT[i]-MAT.bar)/MAT.sdev) ) # Logistic shape function
		y[i]~dnorm(int+exp(a1+b1*(MAT[i]-MAT.bar)/MAT.sdev)/(1+exp(a1+b1*(MAT[i]-MAT.bar)/MAT.sdev))+exp(a2+b2*(MAP[i]-MAP.bar)/MAP.sdev)/(1+exp(a2+b2*(MAP[i]-MAP.bar)/MAP.sdev)),sigma)
	}

	# Prediction
	for (k in 1: length(MAT_x)){
	y_pred_MAT[k] <- int+exp(a1+b1*(MAT_x[k]-MAT.bar)/MAT.sdev)/(1+exp(a1+b1*(MAT_x[k]-MAT.bar)/MAT.sdev))
	}

	# Prediction - MAP
	for (k in 1: length(MAP_x)){
	y_pred_MAP[k] <- int+exp(a2+b2*(MAP_x[k]-MAP.bar)/MAP.sdev)/(1+exp(a2+b2*(MAP_x[k]-MAP.bar)/MAP.sdev))
	}

}