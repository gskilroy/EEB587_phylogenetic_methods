###SIMULATION
#Example of the simulation starter 
#.r =respondability/.r2= coefficient of determination/.f=flexibility/.e=evolvability/.i=integration/.c=conditional evolvability

#Source appropriate files
source("/path/sample.all.hm.R") #sampling function
source("/path/mean.ca.R") #conditional evolvability, analytical solution
source("/path/mean.ea.R") # evolvability, analytical solution
source("/path/mean.f.hm.R") #flexibility, random betas
source("/path/mean.ia.R") #integration, analytical solution
source("/path/mean.r.hm.R") #respondability, random betas
source("/path/t.crit.fxn.R")
source("/path/CalcR2.R") #coefficient of determination

#Required packages
library(clusterGeneration)
library(truncnorm)
library(matrixcalc)
library(foreach)
library(doParallel)
library(moments)


#Start simulation

ind<-150 #number of individuals in the simulated populations
ntraits<-10 #size of the covariance matrix
#number of iterations
iters<-1 #change the number of iterations

G2=list()
repeat{#This is the code to generate the random cov matrix...i used the repeat function to ensure the correct r2
  G1=genPositiveDefMat(ntraits,covMethod=c("c-vine"), eta=1,rangeVar=c(0.5,0.6))
  D=r2(G1$Sigma)
  E=eigen(G1$Sigma)
  G=log(E$values)
  H=skewness(G)
  if (D>=0.17 && D<=0.171 && H>=-0.5){break}
  # change the interval if you want to change the r2 
}


t.crit<-t.crit.fxn(ntraits,ind)

#setup parallel backend to use 8 processors
cl<-makeCluster(6)
registerDoParallel(cl)
out.Final2<-foreach(icount(iters),factor_mult=c(1:100)/20, .export=c('mean.ca','mean.ea','mean.f.hm','mean.ia','mean.r.hm','r2','sample.all.hm','t.crit.fxn','var.ev'),.packages=c('clusterGeneration','truncnorm','moments','matrixcalc')) %dopar% {
  print(iters)
  e.values<-eigen(G1$Sigma)$value
  e.vectors<-eigen(G1$Sigma)$vector
  sum.e.values<-sum(e.values)
  #e.values[1]<-e.values[1]*factor_mult #this is just in case you want to modify magnitude of integration (it's a scaling factor)
  new.sum.e.values<-sum(e.values)
  quotient.ev<-sum.e.values/new.sum.e.values
  e.values<-quotient.ev*e.values
  e.values<-diag(e.values)
  G2<-((e.vectors)%*%e.values%*%t(e.vectors))
  r2.G2<-r2(G2)
  ind.list=c(1:150)
  
  #This part uses another R code to run the simulations
  out.c.G<-sample.all.hm(G2,ntraits,"mean.ca",ind)
  out.e.G<-sample.all.hm(G2,ntraits,"mean.ea",ind)
  out.i.G<-sample.all.hm(G2,ntraits,"mean.ia",ind)
  out.f.G<-sample.all.hm(G2,ntraits,"mean.f.hm",ind)
  out.r2.G<-sample.all.hm(G2,ntraits,"r2",ind)
  out.vc.G<-sample.all.hm(G2,ntraits,"var.ev",ind)
  out.r.G<-sample.all.hm(G2,ntraits,"mean.r.hm",ind)
  
  
  #The remaining of the code is purely getting the output of the simulations and organizing it in a data frame.
  #.r =respondability/.r2= coefficient of determination/.f=flexibility/.e=evolvability/.i=integration/.c=conditional evolvability
  
  output.r.means<-as.numeric(apply(out.r.G[[2]],2,mean)) 
  output.r.sd<-as.numeric(apply(out.r.G[[2]],2,sd))
  ymin.r<-output.r.means[1:ind]-output.r.sd[1:ind]*t.crit[1:ind]
  ymax.r<-output.r.means[1:ind]+output.r.sd[1:ind]*t.crit[1:ind]
  bias.r=((output.r.means-out.r.G[[1]])/out.r.G[[1]])^2
  imp.r=(output.r.sd/out.r.G[[1]])^2
  innac.r=bias.r+imp.r
  lm.r=lm(log(innac.r)~log(ind.list))$coefficients
  a.expo.r=exp(lm.r[1])
  b.expo.r=lm.r[2]
  
  
  output.vc.means<-as.numeric(apply(out.vc.G[[2]],2,mean))
  output.vc.sd<-as.numeric(apply(out.vc.G[[2]],2,sd))
  ymin.vc<-output.vc.means[1:ind]-output.vc.sd[1:ind]*t.crit[1:ind]
  ymax.vc<-output.vc.means[1:ind]+output.vc.sd[1:ind]*t.crit[1:ind]
  bias.vc=((output.vc.means-out.vc.G[[1]])/out.vc.G[[1]])^2
  imp.vc=(output.vc.sd/out.vc.G[[1]])^2
  innac.vc=bias.vc+imp.vc
  lm.vc=lm(log(innac.vc)~log(ind.list))$coefficients
  a.expo.vc=exp(lm.vc[1])
  b.expo.vc=lm.vc[2]
  
  output.r2.means<-as.numeric(apply(out.r2.G[[2]],2,mean))
  output.r2.sd<-as.numeric(apply(out.r2.G[[2]],2,sd))
  ymin.r2<-output.r2.means[1:ind]-output.r2.sd[1:ind]*t.crit[1:ind]
  ymax.r2<-output.r2.means[1:ind]+output.r2.sd[1:ind]*t.crit[1:ind]
  bias.r2=((output.r2.means-out.r2.G[[1]])/out.r2.G[[1]])^2
  imp.r2=(output.r2.sd/out.r2.G[[1]])^2
  innac.r2=bias.r2+imp.r2
  lm.r2=lm(log(innac.r2)~log(ind.list))$coefficients
  a.expo.r2=exp(lm.r2[1])
  b.expo.r2=lm.r2[2]
  
  output.f.means<-as.numeric(apply(out.f.G[[2]],2,mean))
  output.f.sd<-as.numeric(apply(out.f.G[[2]],2,sd))
  ymin.f<-output.f.means[1:ind]-output.f.sd[1:ind]*t.crit[1:ind]
  ymax.f<-output.f.means[1:ind]+output.f.sd[1:ind]*t.crit[1:ind]
  bias.f=((output.f.means-out.f.G[[1]])/out.f.G[[1]])^2
  imp.f=(output.f.sd/out.f.G[[1]])^2
  innac.f=bias.f+imp.f
  lm.f=lm(log(innac.f)~log(ind.list))$coefficients
  a.expo.f=exp(lm.f[1])
  b.expo.f=lm.f[2]
  
  output.i.means<-as.numeric(apply(out.i.G[[2]],2,mean))
  output.i.sd<-as.numeric(apply(out.i.G[[2]],2,sd))
  ymin.i<-output.i.means[1:ind]-output.i.sd[1:ind]*t.crit[1:ind]
  ymax.i<-output.i.means[1:ind]+output.i.sd[1:ind]*t.crit[1:ind]
  bias.i=((output.i.means-out.i.G[[1]])/out.i.G[[1]])^2
  imp.i=(output.i.sd/out.i.G[[1]])^2
  innac.i=bias.i+imp.i
  lm.i=lm(log(innac.i)~log(ind.list))$coefficients
  a.expo.i=exp(lm.i[1])
  b.expo.i=lm.i[2]
  
  output.e.means<-as.numeric(apply(out.e.G[[2]],2,mean))
  output.e.sd<-as.numeric(apply(out.e.G[[2]],2,sd))
  ymin.e<-output.e.means[1:ind]-output.e.sd[1:ind]*t.crit[1:ind]
  ymax.e<-output.e.means[1:ind]+output.e.sd[1:ind]*t.crit[1:ind]
  bias.e=((output.e.means-out.e.G[[1]])/out.e.G[[1]])^2
  imp.e=(output.e.sd/out.e.G[[1]])^2
  innac.e=bias.e+imp.e
  lm.e=lm(log(innac.e)~log(ind.list))$coefficients
  a.expo.e=exp(lm.e[1])
  b.expo.e=lm.e[2]
  
  
  output.c.means<-as.numeric(apply(out.c.G[[2]],2,mean))
  output.c.sd<-as.numeric(apply(out.c.G[[2]],2,sd))
  ymin.c<-output.c.means[1:ind]-output.c.sd[1:ind]*t.crit[1:ind]
  ymax.c<-output.c.means[1:ind]+output.c.sd[1:ind]*t.crit[1:ind]
  bias.c=((output.c.means-out.c.G[[1]])/out.c.G[[1]])^2
  imp.c=(output.c.sd/out.c.G[[1]])^2
  innac.c=bias.c+imp.c
  lm.c=lm(log(innac.c)~log(ind.list))$coefficients
  a.expo.c=exp(lm.c[1])
  b.expo.c=lm.c[2]
  
  
  
  out.Final=rbind(output.r2.means,output.i.means,output.r.means,output.e.means,output.f.means,output.c.means,output.vc.means,
                  ymin.r2,ymin.i,ymin.r,ymin.e,ymin.f,ymin.c,ymin.vc,
                  ymax.r2,ymax.i,ymax.r,ymax.e,ymax.f,ymax.c,ymax.vc,
                  out.r2.G[[1]],out.i.G[[1]],out.r.G[[1]],out.e.G[[1]],out.f.G[[1]],out.c.G[[1]],out.vc.G[[1]],
                  bias.r2,bias.i,bias.r,bias.e,bias.f,bias.c,bias.vc,
                  imp.r2,imp.i,imp.r,imp.e,imp.f,imp.c,imp.vc,
                  innac.r2,innac.i,innac.r,innac.e,innac.f,innac.c,innac.vc,
                  a.expo.r2,a.expo.i,a.expo.r,a.expo.e,a.expo.f,a.expo.c,a.expo.vc,
                  b.expo.r2,b.expo.i,b.expo.r,b.expo.e,b.expo.f,b.expo.c,b.expo.vc)
  out.Final=as.matrix(out.Final)
  out.Final
} 
stopCluster(cl)
#B= Reduce("+", out.Final2) / length(out.Final2) 
save.image(file='Results.RData')


###MEAN CONDITIONAL EVOLVABILITY
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



###MEAN EVOLVABILITY
mean.ea<-function(G1){
  #evolvability, 
  #Analytical mean E
  
  eg.val<-eigen(G1)$values
  mean.e.a.1<-sum(eg.val)/length(eg.val)
  return(mean.e.a.1)
}


###MEAN FLEXIBILITY
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

### ??????
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



###MEAN RESPONDIBILITY
#mean respondibility - mean of 1000 iterations of random betas
mean.r.hm<-function(G1){
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

###CalcR2
CalcR2<-function(cv){
  cor.m<-cov2cor(cv)
  int<-mean(cor.m[lower.tri(cor.m)==T]^2)
  return(int)
}


###
t.crit.fxn<-function(n,p){
  t.crit<-1:p*NA
  #gives the t.crit values for n-p individuals	
  for (i in n:p){
    t.crit[i]<-abs(qt(.025,i-1))
  }
  return(t.crit)
}	


###
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





###HOW INACCCURATE
# This code has three main inputs  
# 1) Ind= Number of individuals sampled (We don't recommend using for N ind<N traits)
# 2) MI = overall level of integration assumed for the trait in question, as measured by the r2 (if this is unknown, we suggest using 0.04 as the lower bound and 0.5 as the upper bound and selecting the largest inaccuracy among the two cases)
# Integration statistics are estimated less accurately at lower levels of integration (r2=0.02), while evolvability statistics tend to be estimated less accurately at higher levels of integration (r2=0.5)
# 3) Ntrait= number of traits measured (we don't recommend using this code for less than 10 or more than 90 traits)

# The model output contains the estimated inaccuracy (measured as an squared proportion of the true parameter value - see Main text), given a certain number of individuals.
# R=respondability, E=evolvability, C= conditional evolvability, I= integration,R2=coefficient of determination, F=flexibility 

howInaccurate<-function(Ind,MI,Ntrait){
  #For each MI and trait number, a power function of the form D*x^C was fitted to the simulation data, relating sampling effort to Inaccuracy. 
  #Symbolic regressions were then used to search for models that describe the relationship between the exponent (C), the constant (D) and our variables (MI and number of traits) using Eureqa (Schmidt and Lipson 2013).
  # The following equations correspond to the calculation of the exponent(C) and the constant (D) as a function of MI and Ntrait for each of the statistics (DR=respondability constant, CR=respondability exponent, for example)
  DR= 0.822043217555314 + 1.02667649509627*MI/(7429^MI - 1.12649405122359) + 88.1381620057372/(0.850024841665574 + (47.1965992371813*MI)^(47.1965992371813*MI + 1.10014941587615^(1.04766577842044^Ntrait)))
  CR= 1.51028131874091*MI + 0.00276645385268465*Ntrait - 1.14081520578666 - DR*MI^2 - 0.194008409646132*log(DR) - 0.0272909673397902*Ntrait*MI^2
  DC= 1.19085679789253*Ntrait^2#updated
  CC= 17.4834144510899/Ntrait^2 - 2.04072681947218#updated
  DE= 3.0478004253639*MI + 0.0272519337861271*MI*Ntrait#updated
  CE= 0.196276772596919*log(MI) - 0.74172233488427 - 0.197130689419141*log(DE)#updated
  CF= 0.00534287803210608*Ntrait + 6.7227280875058*sqrt(MI) - 2.94464879434845 - 4.83017188052178*MI#updated
  DF= -0.00023332500071299*Ntrait*CF^7/MI#updated
  DR2= 1.09375/MI^2#updated
  CR2= 6.69482484957354e-5*DR2 - 0.774908241545569 - MI - 0.188495982853889*log(DR2)#updated
  CI= (15.9856266849758 + 64.4425720681751*MI)/Ntrait^2 - 2.17006400602668 - 0.0173392558459841*log(MI) - 0.265231382533379*MI*log(MI)#updated
  DI= -5.63497435204805 - 11.3730576374658*log(MI)#updated
  
  #The following equations use the exponents and constants calculated above to calculate the amount of inaccuracy (Inacc) 
  Inacc_R=DR*(Ind)^CR
  Inacc_C=DC*(Ind)^CC
  Inacc_E=DE*(Ind)^CE
  Inacc_I=DI*(Ind)^CI
  Inacc_F=DF*(Ind)^CF
  Inacc_R2=DR2*(Ind)^CR2
  
  
  Stats=as.data.frame(rbind(Inacc_R,Inacc_E,Inacc_C,Inacc_I,Inacc_R2,Inacc_F))
  return(Stats)
}


###HOW MANY
# This code has three main inputs 
# 1) Ac= maximum inaccuracy wanted (squared deviation from the mean. e.g., 0.05)
# 2) MI = overall level of integration assumed for the trait in question, as measured by the r2 (if this is unknown, we suggest using 0.04 as the lower bound and 0.5 as the upper bound and selecting the largest N among the two cases)
# Integration statistics are estimated less accurately at lower levels of integration (r2=0.02), while evolvability statistics tend to be estimated less accurately at higher levels of integration (r2=0.5)
# 3) Ntrait= number of traits measured (we don't recommend using this code for less than 10 or more than 90 traits)

# The model output contains the recommended sample sizes (# of individuals) required to calculate each of the statistics with the selected level of inaccuracy
# R=respondability, E=evolvability, C= conditional evolvability, I= integration,R2=coefficient of determination, F=flexibility 

howmany<-function(Ac,MI,Ntrait){
  #For each MI and trait number, a power function of the form D*x^C was fitted to the simulation data, relating sampling effort to Inaccuracy. 
  #Symbolic regressions were then used to search for models that describe the relationship between the exponent (C), the constant (D) and our variables (MI and number of traits) using Eureqa (Schmidt and Lipson 2013).
  # The following equations correspond to the calculation of the exponent(C) and the constant (D) as a function of MI and Ntrait for each of the statistics (DR=respondability constant, CR=respondability exponent, for example)
  DR= 0.822043217555314 + 1.02667649509627*MI/(7429^MI - 1.12649405122359) + 88.1381620057372/(0.850024841665574 + (47.1965992371813*MI)^(47.1965992371813*MI + 1.10014941587615^(1.04766577842044^Ntrait)))
  CR= 1.51028131874091*MI + 0.00276645385268465*Ntrait - 1.14081520578666 - DR*MI^2 - 0.194008409646132*log(DR) - 0.0272909673397902*Ntrait*MI^2
  DC= 1.19085679789253*Ntrait^2#updated
  CC= 17.4834144510899/Ntrait^2 - 2.04072681947218#updated
  DE= 3.0478004253639*MI + 0.0272519337861271*MI*Ntrait#updated
  CE= 0.196276772596919*log(MI) - 0.74172233488427 - 0.197130689419141*log(DE)#updated
  CF= 0.00534287803210608*Ntrait + 6.7227280875058*sqrt(MI) - 2.94464879434845 - 4.83017188052178*MI#updated
  DF= -0.00023332500071299*Ntrait*CF^7/MI#updated
  DR2= 1.09375/MI^2#updated
  CR2= 6.69482484957354e-5*DR2 - 0.774908241545569 - MI - 0.188495982853889*log(DR2)#updated
  CI= (15.9856266849758 + 64.4425720681751*MI)/Ntrait^2 - 2.17006400602668 - 0.0173392558459841*log(MI) - 0.265231382533379*MI*log(MI)#updated
  DI= -5.63497435204805 - 11.3730576374658*log(MI)#updated
  
  #The following equations use the exponents and constants calculated above to calculate the amount of individuals that should be measured to obtain a certain amount of inaccuracy
  IndR=(Ac/DR)^(1/CR)
  IndC=(Ac/DC)^(1/CC)
  IndE=(Ac/DE)^(1/CE)
  IndI=(Ac/DI)^(1/CI)
  IndF=(Ac/DF)^(1/CF)
  IndR2=(Ac/DR2)^(1/CR2)
  
  if(is.nan(IndR)==TRUE) IndR=NA else  if(IndR<=Ntrait+1) IndR=Ntrait+1
  if(is.nan(IndE)==TRUE) IndE=NA else if(IndE<=Ntrait+1) IndE=Ntrait+1
  if(is.nan(IndC)==TRUE) IndE=NA else if(IndC<=Ntrait+1) IndC=Ntrait+1
  if(is.nan(IndR2)==TRUE) IndR2=NA else if(IndR2<=Ntrait+1) IndR2=Ntrait+1
  if(is.nan(IndI)==TRUE) IndI=NA else if(IndI<=Ntrait+1) IndI=Ntrait+1
  if(is.nan(IndF)==TRUE) IndF=NA else if(IndF<=Ntrait+1) IndF=Ntrait+1
  
  Stats=as.data.frame(rbind(IndR,IndE,IndC,IndI,IndR2,IndF))
  return(Stats)
}




