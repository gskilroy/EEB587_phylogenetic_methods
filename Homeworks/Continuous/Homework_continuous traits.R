save(list=ls(), file="ForBrian.rda")

library(ape)
library(geiger)
library(OUwie)
library(phytools)

setwd("/Users/gracekilroy/Continuous")

eval = TRUE

tree<-read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/ilhabela/master/workingFiles/continuousModels/anolis.phy"))

tree<-read.tree("anole.tre")
plot(tree)

continuous.data <- read.csv(file="continuous_data.csv", stringsAsFactors = FALSE)
row.names(continuous.data) <- continuous.data$Species
continuous.data.SVL <- continuous.data$SVL
names(continuous.data.SVL) <- rownames(continuous.data)

str(continuous.data)


cleaned.continuous <-geiger::treedata(phy=tree, data = continuous.data.SVL, sort=TRUE)
visualize.continuous <- geiger::treedata(phy = tree, data = continuous.data)


str(continuous.data.SVL)
print(continuous.data.SVL)



tree.CMap <- contMap(tree,continuous.data.SVL, type = "fan")
phenogram(tree,continuous.data.SVL)





#What is the rate of evolution of your trait on the tree?
BM1 <- geiger::fitContinuous(tree, continuous.data.SVL, model = "BM")
print(paste("The rate of evolution for SLV", "0.136160", "in units of", "meters^2/myr"))
#The rate of evolution for SLV 0.136160 in units of meters^2/myr


OU1 <- fitContinuous(tree, continuous.data.SVL, model="OU")
plot(tree, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU", alpha=0, sigsq=0.136160)
plot(ou.tree)

#How are the trees different? - Using the 'rescale' function changes the branch lengths substantially between the standard OU tree and the scaled OU tree. This is due to a difference in the rate of evolution between the models. 


  
help(package=geiger)

aic.scores <- c(13.400807,15.400807)
names(aic.scores) <- c("BM","OU")
aicw(aic.scores)

AIC.BM1 <- 13.400807
AIC.OU1 <- 15.400807
delta.AIC.BM1 <- 0
delta.AIC.OU1 <- 2

help(package="OUwie")




#######pp.OUwie <- OUwie(phy = tree, data = continuous.data.SVL, model = "OUM", algorithm = "invert", quiet = TRUE)
#come back to this later?


one.discrete.char <- continuous.data.SVL
  
reconstruction.info <- ace(one.discrete.char, tree, type="discrete", method="ML", CI=TRUE)
reconstruction.info$lik.anc
best.states <- colnames(reconstruction.info$lik.anc)[apply(reconstruction.info$lik.anc, 1, which.max)]


plotTree(tree,node.numbers=TRUE)

tree<-paintSubTree(tree,node=185,state="1")
tree<-paintSubTree(tree,node=176,state="2")
tree<-paintSubTree(tree,node=162,state="3")
tree<-paintSubTree(tree,node=151,state="4")
tree<-paintSubTree(tree,node=143,state="5")
tree<-paintSubTree(tree,node=136,state="6")
tree<-paintSubTree(tree,node=123,state="7")
tree<-paintSubTree(tree,node=105,state="8")

plotSimmap(tree,pts=FALSE, node.numbers=TRUE)
phenogram(tree,continuous.data.SVL)


labeled.tree <- 
  
nodeBased.OUMV <- OUwie(tree, cleaned.continuous, model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)



#What do the numbers mean?
  
  Now run all OUwie models:
  
  models <- c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")
results <- lapply(models, RunSingleOUwieModel, phy=tree, data=trait)

AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)


print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model
We get SE for the optima (see nodeBased.OUMV$theta) but not for the other parameters. Let’s see how hard they are to estimate. First, look at ?OUwie.fixed to see how to calculate likelihood at a single point.

?OUwie.fixed
Next, keep all parameters but alpha at their maximum likelihood estimates (better would be to fix just alpha and let the others optimize given this constraint, but this is harder to program for this class). Try a range of alpha values and plot the likelihood against this.

alpha.values<-seq(from= _______________ , to= _______________ , length.out=50)
Keep it simple (and slow) and do a for loop:
  
  likelihood.values <- rep(NA, length(alpha.values))
for (iteration in sequence(length(alpha.values))) {
  likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik
}

plot(x= _______________ , y= _______________, xlab="_______________", ylab="_______________", type="l", bty="n")
points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")
A rule of thumb for confidence for likelihood is all points two log likelihood units worse than the best value. Draw a dotted line on the plot to show this

abline(h=_______________, lty="dotted") #Two log-likelihood
Now, let’s try looking at both theta parameters at once, keeping the other parameters at their MLEs

require("akima")
nreps<-400
theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (iteration in sequence(nreps)) {
  likelihood.values[iteration] <- OUwie.fixed(tree, trait, model="OUMV", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
}
Think of how long that took to do 400 iterations. Now remember how long the search took (longer).

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))
We are interpolating here: contour wants a nice grid. But by centering our simulations on the MLE values, we made sure to sample most thoroughly there

interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400))

contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)

points(x=trait$X[which(trait$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=trait$X[which(trait$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis
The below only works if the discrete trait rate is low, so you have a good chance of estimating where the state is. If it evolves quickly, hard to estimate where the regimes are, so some in regime 1 are incorrectly mapped in regime 2 vice versa. This makes the models more similar than they should be. See Revell 2013, DOI:10.1093/sysbio/sys084 for an exploration of this effect.

library(phytools)
trait.ordered<-data.frame(trait[,2], trait[,2],row.names=trait[,1])
trait.ordered<- trait.ordered[tree$tip.label,]
z<-trait.ordered[,1]
names(z)<-rownames(trait.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,trait,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
print(best)
