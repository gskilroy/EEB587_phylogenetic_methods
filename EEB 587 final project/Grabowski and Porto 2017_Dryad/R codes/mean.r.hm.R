mean.r.hm<-function(G1){
	#mean respondibility - mean of 1000 iterations of random betas
	it<-1000
	saved.r<-1:it*0
	for (i in 1:it){
		r.beta<-rnorm(dim(G1)[2],0,1)
		r.beta<-r.beta/sqrt(sum(r.beta^2))
		saved.r[i]<-sqrt(t(r.beta)%*%(G1%*%G1)%*%r.beta)
		}
	mean.r<-sum(saved.r)/(length(saved.r))
	return(mean.r)
	}