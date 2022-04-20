library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
devtools::install_github("thej022214/hisse")
library(hisse)

?sim.bd.taxa

my.tree <- TreeSim::sim.bd.taxa(n=30, numbsim=1, lambda=0.1, mu=0)[[1]]
ape::ltt.plot(my.tree)
ape::ltt.plot(my.tree, log="y")
yule.trees <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=0.1, mu=0, complete=FALSE)
ape::mltt.plot(yule.trees, log = "y")


bd.trees <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=1, mu=.9, complete=FALSE)
ape::mltt.plot(bd.trees, log="y", legend=FALSE)

depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range))
plot(x=c(0, -1*max.depth), y=c(1, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)


depth.range <- range(unlist(lapply(yule.trees,ape::branching.times)), unlist(lapply(bd.trees,ape::branching.times)))
max.depth <- sum(abs(depth.range))
plot(x=c(0, -5), y=c(20, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)





#What happens if speciation rate is much higher than extinction rate? 
eval=TRUE
my.trees.1 <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=1, mu=.2, complete=FALSE)
ape::mltt.plot(my.trees.1, log="y", legend=FALSE)



#How does the simulation change with different values, but keeping their difference constant?
eval=TRUE
my.trees.2 <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=.8, mu=.4, complete=FALSE)
ape::mltt.plot(my.trees.2, log="y", legend=FALSE)

my.trees.3 <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=.5, mu=.1, complete=FALSE)
ape::mltt.plot(my.trees.3, log="y", legend=FALSE)



#If their sum is constant?
eval=TRUE
my.trees.4 <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=.8, mu=.2, complete=FALSE)
ape::mltt.plot(my.trees.4, log="y", legend=FALSE)

my.trees.5 <- TreeSim::sim.bd.taxa(n=30, numbsim=10, lambda=.6, mu=.4, complete=FALSE)
ape::mltt.plot(my.trees.5, log="y", legend=FALSE)




speciation.rates <- c(0.1, 0.1, 0.1, 0.2)
extinction.rates <- rep(0.03, 4)
transition.rates <- c(0.01,0.01,0, 0.01, 0, 0.01, 0.01,0,0.01, 0,0.01,0.01)
pars <- c(speciation.rates, extinction.rates, transition.rates)
phy <- tree.musse(pars, max.taxa=50, x0=1, include.extinct=FALSE)
sim.dat.true <- data.frame(names(phy$tip.state), phy$tip.state)
sim.dat <- sim.dat.true
# Now to hide the "hidden" state
sim.dat[sim.dat[,2]==3,2] = 1
sim.dat[sim.dat[,2]==4,2] = 2
# and convert states 1,2 to 0,1
sim.dat[,2] = sim.dat[,2] - 1

plot(phy)
knitr::kable(cbind(sim.dat, true.char=sim.dat.true$phy.tip.state))








turnover.anc = c(1,1,2,2)
eps.anc = c(1,1,2,2)



trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))

trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

# Now we want three specific rates:
trans.rates.nodual.threerates <- trans.rates.nodual
# Set all transitions from 0->1 to be governed by a single rate:
to.change <- cbind(c(1,3), c(2,4))
trans.rates.nodual.threerates[to.change] = 1
# Now set all transitions from 1->0 to be governed by a single rate:
to.change <- cbind(c(2,4), c(1,3))
trans.rates.nodual.threerates[to.change] = 2
# Finally, set all transitions between the hidden state to be a single rate (essentially giving
# you an estimate of the rate by which shifts in diversification occur:
to.change <- cbind(c(1,3,2,4), c(3,1,4,2))
trans.rates.nodual.threerates[to.change] = 3
trans.rates.nodual.threerates


pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)


load("testrecon1.rda")
class(pp.recon)
pp.recon


plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE)


plot.hisse.states(pp.recon, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))

pp.recon$aic

pp.recon = MarginRecon(phy, sim.dat, f=c(1,1), hidden.states=TRUE, pars=pp$solution, aic=pp$aic, n.cores=2)



hisse.results.list = list()
load("testrecon1.rda")
hisse.results.list[[1]] = pp.recon
load("testrecon2.rda")
hisse.results.list[[2]] = pp.recon
load("testrecon3.rda")
hisse.results.list[[3]] = pp.recon
# Now supply the list the plotting function
plot.hisse.states(hisse.results.list, rate.param="net.div", show.tip.label=FALSE, rate.range=c(0,0.072))


files = system("ls -1 | grep .rda", intern=TRUE)
# Create an empty list object
hisse.results.list = list()
# Now loop through all files, adding the embedded pp.recon object in each
for(i in sequence(length(files))){
  load(files[i])
  hisse.results.list[[i]] = pp.recon
  rm(pp.recon)
}






















##turnover.anc has a single free parameter for both 0A and 1A state combinations; there is a single free parameter for extinction fraction. This is equivalent to a BiSSE model with a fixed turnover and extinction rates across the observed states 0 an 1. Now, say we want to include separate turnover rates for both states 0A and 1A:
turnover.anc = c(1,1,0,0)
eps.anc = c(1,1,0,0)


##include separate turnover rates for both states 0A and 1A:
turnover.anc = c(1,2,0,0)

##full hisse model (corresponds to four separate net turnover rates for 1=0A, 2=1A, 3=0B, 4=1B):
turnover.anc = c(1,2,3,4)


##Yule (pure birth rate):
eps.anc = c(0,0,0,0)



#To generate the index matrix describing the free parameters in the transition model, use the TransMatMaker() function:
trans.rates = TransMatMaker.old(hidden.states=TRUE)
trans.rates


#drop dual transitions between both the observed trait and the hidden trait:
trans.rates.nodual = ParDrop(trans.rates, c(3,5,8,10))
trans.rates.nodual


#set parameter 6 equal to parameter 1:
trans.rates.nodual.equal16 = ParEqual(trans.rates.nodual, c(1,6))
trans.rates.nodual.equal16


#set all rates to be equal:
trans.rates.nodual.allequal = ParEqual(trans.rates.nodual, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))
trans.rates.nodual.allequal

#other way to do this:
trans.rates.nodual.allequal = trans.rates.nodual
trans.rates.nodual.allequal[!is.na(trans.rates.nodual.allequal) & !trans.rates.nodual.allequal == 0] = 1
trans.rates.nodual.allequal

#to run a BiSSE model in HiSSE, the matrix set up would look like this:
trans.rates.bisse = TransMatMaker.old(hidden.states=FALSE)
trans.rates.bisse

#Whatever transition matrix is designed, it is supplied to the trans.rate= argument in the hisse() call:
pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal)

#to output net diversification: 
pp = hisse.old(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.allequal, output.type="net.div")



















