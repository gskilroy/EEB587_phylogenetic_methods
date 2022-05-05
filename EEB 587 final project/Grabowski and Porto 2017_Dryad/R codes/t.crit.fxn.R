t.crit.fxn<-function(n,p){
	t.crit<-1:p*NA
	#gives the t.crit values for n-p individuals	
	for (i in n:p){
		t.crit[i]<-abs(qt(.025,i-1))
		}
	return(t.crit)
	}	