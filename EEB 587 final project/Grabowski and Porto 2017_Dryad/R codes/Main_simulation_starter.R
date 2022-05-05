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
