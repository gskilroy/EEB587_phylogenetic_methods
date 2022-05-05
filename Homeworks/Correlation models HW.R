library(geiger)
library(ape)
tree.primates <- read.tree(text="((((Homo:0.21,Pongo:0.21):0.28,Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);")
X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
pic.X <- pic(X, tree.primates)
pic.Y <- pic(Y, tree.primates)


require("corHMM")
?corHMM
data(primates)
ls()
print(primates)
require(phytools)



primates$trait[which(grepl("Hylobates",primates$trait[,1])),2]<-1

trait1<-primates$trait[,2]
names(trait1)<-primates$trait[,1]
primates$tree <- ape::multi2di(primates$tree)
plotSimmap(make.simmap(primates$tree, trait1), pts=FALSE, fsize=0.8)
rate.mat.er<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=2, model="ER")
print(rate.mat.er)

#what does this matrix mean?
#This matrix illustrates a single trait with two discrete states. There is a single rate category by which the trait can change from state 1 to state 2. The 'NAs' function as place holders indicating that the same character state cannot change to itself. 



pp.er<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=rate.mat.er,node.states="marginal")
print(pp.er)


#what do these results mean? 
#The rate parameter is the same for both character states.



pp.ard<-corHMM(primates$tree,primates$trait[,c(1,2)],rate.cat=1,rate.mat=NULL, model= "ARD", node.states="marginal")
print(pp.ard)

#what do these results mean? 
#In this "ARD" model the rate parameter estimate for the two character states differs slightly. 


#which model is better?
#the first model has a lower AIC score and should be considered to be a better fit for this dataset. 




rate.mat.er.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ER")
print(rate.mat.er.4state)

fourstate.trait<-rep(NA,Ntip(primates$tree))
for(i in sequence(Ntip(primates$tree))) {
  if(primates$trait[i,2]==0 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-0
  }
  if(primates$trait[i,2]==0 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-1
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==0) {
    fourstate.trait[i]<-2
  }
  if(primates$trait[i,2]==1 && primates$trait[i,3]==1) {
    fourstate.trait[i]<-3
  }
}
fourstate.data<-data.frame(Genus_sp=primates$trait[,1], T1=fourstate.trait)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, model="ER", node.states="marginal"))
print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat=rate.mat.er.4state, node.states="marginal", model="ARD"))
rate.mat.ard.4state<-corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=1, nstates=4, model="ARD")
print(rate.mat.ard.4state)



rate.mat.gtr.4state<-rate.mat.ard.4state
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(1,4))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(2,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(3,8))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(4,6))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(5,7))
rate.mat.gtr.4state<-corHMM:::rate.par.eq(rate.mat.gtr.4state, c(6,7))
print(rate.mat.gtr.4state)

print(rayDISC(primates$tree, fourstate.data, ntraits=1, rate.mat= rate.mat.gtr.4state, node.states="marginal", model="ARD"))


print(corHMM:::rate.mat.maker(rate.cat=1, hrm=FALSE, ntraits=2, nstates=2, model="ARD"))
rate.mat.pag94<-corHMM:::rate.par.drop(rate.mat.ard.4state, drop.par=c(3,5,8,10))
print(rate.mat.pag94)




#Construct a model to test if state 1 can never be lost
legend_rate <- getStateMat4Dat(primates$trait,model = "ARD")
RateCat1 <- equateStateMatPars(legend_rate, c(1:4))
RateCat1 <- dropStateMatPars(legend_rate$rate.mat, c(1,4))

#I'm struggling with figuring out how to properly construct this model. I will continue working on this over the weekend. 


#Experiment with the effects of frequencies at the root.



#Create and use a model to see if transitions from 00 go to 11 only via 01





