###EEB 587 FINAL PROJECT

##Machado et al 2018 - evolution of morphological integetration in Carnivora facial morphology - changes in Canidae lead to increased evolutionary potential of traits



##Question 1 - explore necessary sample size for exploring evolvability from integration magnitudes with number of traits used, number of "species (aka age groups)"
#previous research has shown that in carniova facial morphology changes correlate 

#power analysis
#Grabowski and Porto 2017

##Question 2 - effect of standardizing by size in integration analysis - compare non-scaled and scaled to evolutionary tree 
#Machado et al 2019 - measuring magnitude of morphological integration




#Process-

#what is the sample size necessary to evaluate changes in magnitudes of integration how does standardizing by size effect

#STEP 1 - replicate Machado et al 2018 (used Carniva), make sure R code is working and methods make sense
#STEP 2 - explore standardizing Machado's data by size, how does evolvability change, compare with known tree
#STEP 3 - theoretical aspect - lay out my research questions and how the above results inform my research questions/methods. Use Grabowski and Porto findings to run a power analysis for my data informed by Machado Carnivora results


#create variancce covariance matrix

#pull in magnitude data from various papers for taxa
#compare against tree











##Rolian 2009 paper - looked at integration and evolvability in primate hands and feet

#apply to outside dataset, Rolian 2009???? - can I get this? or ask Ben for data? 









##Grabowski and Porto r functions:
#howmany.R
#howInaccurate.R


#step 1 - simulate population based on each covariance matrix
  #10,000 draws from a multivariate normal distribution based on the simulated matrices with null mean
  #samples were drawn for all sample sizes capable of producing full rank matrices, a covariance matrix was estimated fro each, statistics for each were calculated
  #repeated 100 times
  #sample size increased by 1 and procedure repeated

#step 2 - calculate parameter values of each statistic based on known population covariance matrix
#step 3 - calculate covariance matrices based on random sample of individuals from the main population 
#step 4 - calculate statistics of interest for each of the sampled matrices and compare the values of each to the known parameter values


#additional step for calculating evolvability statistic (rely on the average response of covariance matrices to simulated selection vectors):
  #create random selection vectors by drawing from a random normal distribution with a mean of 0 and a SD of 1, normalized to length, apply to the covariance matrix for each sample using the equations in Table 1 of Hansen & Houle 2008 to calculate statistic of interest
  #Matrix inversion in highly multidimensional systems is time consuming - results for mean conditional evolvability and mean integration for PAT 1000 were produced using analytical approximations from Hansen & Houle rather than simulation 



#Determining sample size
  #Grabowski study assumes researcher is studying 10 traits with specific magnitudes of integration 
  #Using different numbers of traits (10 to 100) and different magnitudes of MI (r2 0.02 to 0.5)
  #a power function of the form ax^b was fitted to the simulation data for each MI and trait number - this relates sampling effort to inaccuracy, used EUREQA (Schmidt & Lipson 2013)
    #a = constant
    #b = exponent
    #x = variables (MI and number of traits)
  #



#changes in magnitudes of integration correlate with potential changes in evolvability  - look at how  standardizing by size effects these results. Compare 

#looking at changing magnitudes of integration correlated with evolvability - this has been looked at a bit in non-human primates






num_ind <- 150
num_traits <- 30




library(ape)
library(phangorn)

tree <- rtree(10)
allDescendants(tree)
Descendants(tree, ((1+ape::Ntip(tree)):ape::Nnode(tree)), type="tips")



install.packages("phylobase")
library(phylobase)







####
primate_tree <- read.tree(file = "primate_newick")
plot(primate_tree)











