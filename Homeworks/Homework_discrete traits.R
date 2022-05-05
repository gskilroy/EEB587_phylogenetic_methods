##HOMEWORK - DISCRETE TRAITS

getwd()
setwd("/Users/gracekilroy/Discrete")

library(devtools)
install_github("bomeara/geiger-v2")
devtools::install_github("bomeara/geiger-v2", force = TRUE)
library(geiger)
library(ape)
library(phangorn)

tree <- read.tree(file = "primate_newick")
plot(tree)
tree$tip.label


discrete.data <- read.csv(file="primate_data.csv", stringsAsFactors=FALSE)
row.names(discrete.data) <- discrete.data$Tip.Label
discrete.data$Tip.Label

CleanData <- geiger::treedata(phy = tree, data = discrete.data, sort = TRUE)
VisualizeData <- geiger::treedata(phy = tree, data = discrete.data)



#parsimony
cleaned.discrete.phyDat <- phangorn::phyDat(CleanData$data, type = "USER", levels = c('0', '1', 'Folivore', 'Omnivore', 'Frugivore', 'Polygyny', 'Pair', 'Polygynandry', 'Solitary', 'Harem polygyny', 'unknown', 'Monogamy', 'Spatial polygyny')) 
anc.p <- phangorn::ancestral.pars(CleanData$phy, cleaned.discrete.phyDat)
plotAnc(tree, anc.p, 1)


#likelihood
anc.ml <- ancestral.pml(pml(CleanData$phy, cleaned.discrete.phyDat), type = "ml")
plotAnc(tree,anc.ml,1) 

print(discrete.data)

install.packages("corHMM")

###How can you estimate transition rates between states? Do it.
#You can estimate transition rates between states by using the corHMM package. 
library(corHMM)
?corHMM
LegendAndRateMat <- getStateMat4Dat(discrete.data)
RateMat <- LegendAndRateMat$rate.mat
RateMat
pars2equal <- list(c(1,2), c(3,4), c(5,6))
StateMatA_constrained <- equateStateMatPars(RateMat, pars2equal)
#This set of commands creates a custom model for a estimating rate matrix based on the given data.

transition.rate <- corHMM(phy = CleanData$phy, data = CleanData$data, rate.cat = 7)




###How could you examine if transition rates are equal?
RateCat1 <- getStateMat4Dat(discrete.data)$rate.mat
RateCat1 <- equateStateMatPars(RateCat1, c(1:6))

RateCat2 <- getStateMat4Dat(discrete.data)$rate.mat
RateCat2 <- dropStateMatPars(RateCat2, 3)

StateMats <- list(RateCat1, RateCat2)

FullMat <- getFullMat(StateMats, RateClassMat)

plotMKmodel(FullMat, rate.cat = 2, display = "row", text.scale = 0.7)



###Think about the Lewis (2001) MKV model. Are your traits all variable? Will using this make sense for your data? Try using it. Do results change?
#The discrete traits listed in this dataset are variable across all of the taxa. This can be problematic rates of evolution can be overestimated when only all traits in a dataset are variable across taxa. Lewis (2001) addresses this problem with the MKV model which computes a conditional likelihood based on only variable characters.  

lewis.model <- rayDISC(phy = tree, data = discrete.data, rate.mat = NULL, node.states = "none", lewis.asc.bias = TRUE, root.p = "yang")



###How could you test order of state evolution?
#Using lambda as described by Pagel describeds the best fit distribution of data. Lambda produces values between 0 and 1 where 1 represents phylogenetic independence and 1 means that traits are evolving by a Brownian process. Lambda is found in the geiger package.

