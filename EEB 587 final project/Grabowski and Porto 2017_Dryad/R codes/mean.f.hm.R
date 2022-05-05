mean.f.hm<-function(G1){
#flexibility =mean of 	#1000 iterations of random betas
	it<-1000
	saved.f.1<-1:it*0
	for (i in 1:it){
		r.beta<-rnorm(dim(G1)[2],0,1)
		r.beta<-r.beta/sqrt(sum(r.beta^2))
		evol.response<-G1%*%r.beta
		evol.response=evol.response/sqrt(sum(evol.response^2))
		#print(r.beta)
		#print(G1)
		saved.f.1[i]<-sum(r.beta*evol.response)
		}
	return(sum(saved.f.1)/(length(saved.f.1)))
	}