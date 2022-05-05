devtools::install_github("bomeara/phybase")
install.packages("devtools")
library(rotl)
library(ape)
phy <- get_study_tree("ot_485", "tree1")
plot(phy, cex=0.3)
library(geiger)
phy <- drop.random(phy, Ntip(phy) - 10)
plot(phy)
axisPhylo()
library(phybase)
gene.tree <- phybase::sim.coaltree.phylo(phy, pop.size=1e-12)
plot(gene.tree)
library(phytools)
plot(cophylo(phy, gene.tree, cbind(sort(phy$tip.label), sort(gene.tree$tip.label))))
species.tree <- rcoal(7)
species.tree$edge.length <- species.tree$edge.length / (10*max(branching.times(species.tree)))
gene.tree <- phybase::sim.coaltree.phylo(species.tree)
plot(cophylo(species.tree, gene.tree, cbind(sort(species.tree$tip.label), sort(gene.tree$tip.label))))
tip.rows <- which(species.tree$edge[,2]<=Ntip(species.tree))
species.tree2 <- species.tree
species.tree2$edge.length[tip.rows] <- 100 + species.tree2$edge.length[tip.rows]
gene.tree2 <- phybase::sim.coaltree.phylo(species.tree2)
plot(cophylo(species.tree2, gene.tree2, cbind(sort(species.tree2$tip.label), sort(gene.tree2$tip.label))))
species.tree2.clado <- compute.brlen(species.tree2)
gene.tree2.clado <- compute.brlen(gene.tree2)
plot(cophylo(species.tree2.clado, gene.tree2.clado, cbind(sort(species.tree2.clado$tip.label),
                                                          sort(gene.tree2.clado$tip.label))))
devtools::install_github("bomeara/taxon2tree")
1
devtools::install_github("bomeara/taxon2tree")

install.packages("devtools")
devtools::install_github("bomeara/taxon2tree")
install.packages("phylotaR")
c
install_github("ropensci/phylotaR")
devtools::install_github("ropensci/phylotaR")
devtools::install_github("bomeara/taxon2tree")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")
BiocManager::install("Biostrings", ask=FALSE) 
devtools::install_github("bomeara/taxon2tree")

library(taxon2tree)
BatchRun("Myrmecocystus")
install.packages("visNetwork")
BatchRun("Myrmecocystus", raxml="/usr/local/bin/ncbi_path")
BatchRun("Myrmecocystus", raxml="/usr/local/bin/raxml-ng")
getwd()


BatchRun("Myrmecocystus")
install.packages("ape","phangorn", "phylotaR", "rentrez", "targets")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
install.packages("ape")
install.packages("phangorn")
install.packages("phylotaR")
install.packages("rentrez")
install.packages("targets")
library(taxon2tree)

BatchRun("vertebrates_fgf.fasta.nex")
library(ape)
library(phangorn)
library(phylotaR)
library(rentrez)
library(targets)
library(rotl)
force=TRUE  
devtools::install_github("bomeara/phybase", force = TRUE) 


BatchRun("Myrmecocystus", raxml="/usr/local/bin/raxml-ng")
library(Biostrings)

blastSeqKK <- function (x, database = "nr", hitListSize = "10", 
                        filter = "L", expect = "10", program = "blastn",
                        attempts = 10) {
  baseUrl <- "http://www.ncbi.nlm.nih.gov/blast/Blast.cgi"
  query <- paste("QUERY=", as.character(x), "&DATABASE=", database, 
                 "&HITLIST_SIZE=", hitListSize, "&FILTER=", filter, "&EXPECT=", 
                 expect, "&PROGRAM=", program, sep = "")
  url0 <- sprintf("%s?%s&CMD=Put", baseUrl, query)
  results <- tempfile()
  Sys.sleep(5)
  require(XML)
  post <- htmlTreeParse(url0, useInternalNodes = TRUE)
  x <- post[["string(//comment()[contains(., \"QBlastInfoBegin\")])"]]
  rid <- sub(".*RID = ([[:alnum:]]+).*", "\\1", x)
  rtoe <- as.integer(sub(".*RTOE = ([[:digit:]]+).*", "\\1", 
                         x))
  url1 <- sprintf("%s?RID=%s&FORMAT_TYPE=XML&CMD=Get", baseUrl, 
                  rid)
  Sys.sleep(rtoe)
  .tryParseResult <- function(url, attempts){
    for (i in 1:(attempts+1)) {
      result <- tryCatch({
        xmlTreeParse(url, useInternalNodes=TRUE,
                     error = xmlErrorCumulator(immediate=FALSE))
      }, error=function(err) NULL)
      if (!is.null(result)) return(result)
      Sys.sleep(10)
    }
    stop(paste("no results after ", attempts, 
               " attempts; please try again later", sep = ""))
  }
  result <- .tryParseResult(url1, attempts)
  qseq <- xpathApply(result, "//Hsp_qseq", xmlValue)
  hseq <- xpathApply(result, "//Hsp_hseq", xmlValue)
  require(Biostrings)
  res <- list()
  for (i in seq_len(length(qseq))) {
    res[i] <- DNAMultipleAlignment(c(hseq[[i]], qseq[[i]]), 
                                   rowmask = as(IRanges(), "NormalIRanges"), colmask = as(IRanges(), 
                                                                                          "NormalIRanges"))
  }
  res
}

devtools::install_bioc("Biostrings")
  
devtools::install_github("mhahsler/rBLAST")
library(Biostrings)
library(rBLAST)
? rBLAST
? blast


install.packages("devtools")
devtools::install_github("bomeara/taxon2tree", force=TRUE)


getwd()
setwd("/Users/gracekilroy/Desktop/EEB587_phylogenetic_methods/Tree building")

library(taxon2tree)
BatchRun("Myrmecocystus", ncbi_path="/Users/gracekilroy/Downloads/ncbi/bin")






if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")


BiocManager::install("Biostrings", force=TRUE)
library(taxon2tree)
library("ape")
library("phangorn")
library("phylotaR")
library("rentrez")
library("targets")



devtools::install_github("bomeara/phybase", force = TRUE)

getwd()
setwd("/Users/gracekilroy/Taxon2Tree")

BatchRun("Myrmecocystus", ncbi_path="/Users/gracekilroy/Downloads/ncbi/bin")



setwd("/Users/gracekilroy/T2T")

BatchRun("Petrogale", ncbi_path="/Users/gracekilroy/Downloads/ncbi/bin", mintaxa = 5)


library(taxon2tree)


library(ape)
Petrogale.tree<-system.file("RAxML_bestTree.combined", package="ape")
plot("RAxML_bestTree.combined")


plot_phylota_treemap(phylota = Petrogale.tree, cids = , txids = NULL)


Petrogale_tree <- read.tree(file = "RAxML_bestTree.combined")
plot(Petrogale_tree)


library(phylotaR)
?phylotaR
??phylotaR



?r8s.phylo
install.packages("r8s.phylo")
install.packages("treepl")


BatchRun("Simiiformes", ncbi_path="/Users/gracekilroy/Downloads/ncbi/bin", mintaxa = 5)




##HOMEWORK - DISCCRETE TRAITS


getwd()
setwd("/Users/gracekilroy/Discrete")

library(devtools)
install_github("bomeara/geiger-v2")
devtools::install_github("bomeara/geiger-v2", force = TRUE)
library(geiger)

tree <- read.tree(file = "primate_newick")
plot(tree)
discrete.data <- read.csv(file="primate_data.csv", stringsAsFactors=FALSE)

?geiger::treedata

CleanData <- geiger::treedata(phy = tree, data = discrete.data, sort = TRUE)
print(ape::Ntip(CleanData$phy))
print(nrow(CleanData$data))

tree$tip.label
discrete.data


VisualizeData <- geiger::treedata(phy = tree, data = discrete.data)


  tree<-read.tree(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/squamate.phy"))
plot(tree)
discrete.data<-read.csv(text=RCurl::getURL("https://raw.githubusercontent.com/lukejharmon/pcm/master/datafiles/brandley_table.csv"), stringsAsFactors=FALSE)
rownames(discrete.data) <- gsub(" ", "_", discrete.data$Species)

discrete.data.corrected<-as.data.frame(t(discrete.data))


cleaned.discrete <- geiger::treedata(phy=tree, data = discrete.data, sort=TRUE)



library(phangorn)

?phangorn::phyDat

#parsimony
cleaned.discrete.phyDat <- phangorn::phyDat(CleanData$data, type = "USER", levels = c(0,1,2)) 
anc.p <- phangorn::ancestral.pars(CleanData$phy, cleaned.discrete.phyDat)
plotAnc(tree, anc.p, 1)

#likelihood
anc.ml <- ancestral.pml(pml(tree, cleaned.discrete.phyDat), type = "ml")
plotAnc(tree,anc.ml,1)  




transition.rate <- corHMM(phy = CleanData$phy, data = CleanData$data, rate.cat = 1)








#class example
unif_nums <- rnorm(1000)
hist(unif_nums)


unif_nums <- runif(1000)
hist(unif_nums)
unifsum <- function(n=1000) {
  sum(runif(n))
}

many <- replicate(500, unifsum())
hist(many)


unif_nums <- rexp(1000)
hist(unif_nums)
unifsum <- function(n=1000) {
  sum(rexp(n))
}

many <- replicate(500, unifsum())
hist(many)




phy <- ape::rcoal(10)
plot(phy)
trait <- sim.char(phy, par=0.1, nsim=1, model="BM", root=0)[,,1]
fitContinuous(phy,trait)
help(packag="phytools")
contMap(phy, trait)


forwardSim <- function(x) {
  return(x+(.01*x*runif(1)))
}


x_vector <- rep(NA, 1000)
x_vector[1] <- 1
for (i in sequence(length(x_vector)-1)) {
  x_vector[i+1] <- forwardSim(x_vector[i])
}

plot(x_vector)
plot(log(x_vector))



###send to Brian
save(list=ls(), file="ForBrian.rda")
####

tree



#FINAL PROJECT NOTES
##develop one comparison method
#have a good question
#work to solve it
#paper came out a few years ago that had a radically new dino phylogeny
#simulations next week will make this easier
#make a github project and put everything there

#how to traits effect branch lengths?
#beast - rate estimate for each branch

#do a power analysis with made up data
#this is helpful to use for DDRIG proposal











##Diversificaiton HW notes
library(ape)
library(TreeSim)
library(geiger)
library(diversitree)
devtools::install_github("thej022214/hisse")
library(hisse)


my.tree <- TreeSim::sim.bd.taxa(n=300, numbsim=1, lambda=0.1, mu=0)[[1]]
ape::ltt.plot(my.tree)
ape::ltt.plot(my.tree, log="y")
yule.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=0.1, mu=0, complete=FALSE)
ape::mltt.plot(yule.trees, log = "y")


bd.trees <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.9, complete=FALSE)
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
plot(x=c(0, -5), y=c(200, ape::Ntip(yule.trees[[1]])), log="y", type="n", bty="n", xlab="Time", ylab="N")
colors=c(rgb(1,0,0,0.5), rgb(0, 0, 0, 0.5))
list.of.both <- list(bd.trees, yule.trees)
for (i in sequence(2)) {
  tree.list <- list.of.both[[i]]
  for (j in sequence(length(tree.list))) {
    ape::ltt.lines(tree.list[[j]], col=colors[[i]])
  }
}
legend("topleft", legend=c("Birth Death", "Yule"), fill=colors)


?sim.bd.taxa


#What happens if speciation rate is much higher than extinction rate? 
eval=TRUE
my.trees.1 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=1, mu=.2, complete=FALSE)
ape::mltt.plot(my.trees.1, log="y", legend=FALSE)



#How does the simulation change with different values, but keeping their difference constant?
eval=TRUE
my.trees.2 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=.8, mu=.4, complete=FALSE)
ape::mltt.plot(my.trees.2, log="y", legend=FALSE)

my.trees.3 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=.5, mu=.1, complete=FALSE)
ape::mltt.plot(my.trees.3, log="y", legend=FALSE)



#If their sum is constant?
eval=TRUE
my.trees.4 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=.8, mu=.2, complete=FALSE)
ape::mltt.plot(my.trees.4, log="y", legend=FALSE)

my.trees.5 <- TreeSim::sim.bd.taxa(n=300, numbsim=10, lambda=.6, mu=.4, complete=FALSE)
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




##use targets r package to keep everything in line
#can be hard to get into but so much easier once you do

install.packages("targets")






