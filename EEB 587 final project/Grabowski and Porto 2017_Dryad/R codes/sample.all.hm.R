q r2<-function(cv){
	cor.m<-cov2cor(cv)
	int<-mean(cor.m[lower.tri(cor.m)==T]^2)
	return(int)
	}

var.ev<-function(cv){
	cor.m<-cov2cor(cv)
	#int<-sum((eigen(cor.m)$values-1)^2)/length(eigen(cor.m))
	int<-var(eigen(cor.m)$values)/(length(eigen(cor.m)$values)-1)
	return(int)	
	}

 sample.all.hm<-function(cv1,num.traits,anal.type,ind){
 	pop.size<-10000
	reps<-100
	ev1.dist<-NULL
	#Create vector of means
	means<-rep(0,num.traits)
	
	
	fun<-match.fun(anal.type)
	store.hm<-matrix(NA,reps,ind)
	store.bias.hm<-matrix(NA,1,ind)
	store.imp.hm<-matrix(NA,1,ind)	
	store.inac.hm<-matrix(NA,1,ind)



	#Calculate parameter P based on n=pop.size
	pop.dist<-mvrnorm(pop.size,means,cv1)
	pop.cv<-cov(pop.dist,use="pairwise.complete.obs")
	parameter.hm<-fun(pop.cv)
		
	#print("OK")
	for (n in 5:ind){
		#print(paste("Individual #",n))
		for (p in 1:reps){
			sampled.data<-sample(1:pop.size,n,replace=FALSE)
			sampled.data<-pop.dist[sampled.data,]
			sampled.cv<-cov(sampled.data,use="pairwise.complete.obs")

			
			if (is.positive.definite(sampled.cv)==FALSE){
				next()
				}
							
			store.hm[p,n]<-fun(sampled.cv)	

			}
		
			
			
		}

	out<-list(parameter.hm,store.hm)
	return(out)

	}
