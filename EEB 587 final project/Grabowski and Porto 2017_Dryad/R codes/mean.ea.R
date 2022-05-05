mean.ea<-function(G1){
	#evolvability, 
	#Analytical mean E
	
	eg.val<-eigen(G1)$values
	mean.e.a.1<-sum(eg.val)/length(eg.val)
	return(mean.e.a.1)
	}