mean.ca<-function(G1){
	#conditional evolvability,
	#Analytical mean cE
	
	k <- nrow(G1)
	eg.val<-eigen(G1)$values
	H<-1/mean(1/eg.val)
	I <- var(1/eg.val)/(mean(1/eg.val)^2)
	mean.ce.a.1<-H*(1+(2*I/(k+2)))
	return(mean.ce.a.1)
	}