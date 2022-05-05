CalcR2<-function(cv){
	cor.m<-cov2cor(cv)
	int<-mean(cor.m[lower.tri(cor.m)==T]^2)
	return(int)
	}
