#integration-analytical derivation
mean.ia<-function(G1){ 
	eg.val<-eigen(G1)$values
	k <- nrow(G1)
	E <- mean(eg.val)
	H<-1/mean(1/eg.val)
	I <- var(eg.val)/(mean(eg.val)^2)
	I2<-var(1/eg.val)/(mean(1/eg.val)^2)
	mean.a<-(H/E)*(1+2*(I+I2-1+H/E+2*I*I2/(k+2))/(k+2))
	mean.i<-1-mean.a
	return(mean.i)
	}