#SIMULATION

#Some taxa are endangered, assume this is binary (endangere or not)
#Endangered taxa have an extinction risk, expected lifespan is exponential, based on the extinction rate
#Non-endangered taxa have an extinction risk too but smaller
#Assume no speciation (optimism-free model)
#Parameter we care about is number of species going extinct
#Assume it is completely random
#latin hypercube sampling is nice, we're going to do a grid sample instead for simplicity
#Plot as a heatmap, height/color is number of species going extinct. X axis is extinction_rate_EN, Y axis is number_EN
#Exponential: N(t) = N(0) * exp(r*T)
#Solving for r:
#     N(t)/N(0) = exp(r*T)
#     0.75 = exp(r*T) #assume 75% survive
#     0.75 = exp(r*100)
#     ln(0.75) = r * 100
#     r = ln(0.75) / 100
log(0.75)/100
#     The slowest rate of extinction is -0.003 

#     r = ln(0.01) / 100
log(0.01)/100
#     Highest rate of extinction could be -0.046


#for species of least concern, have it range from 50% to all survive


#parameters (need to vary these):
extinction_rate_EN <- NA
extinction_rate_LC <- NA
number_all_species <- 100
proportion_EN <- NA
total_time <- 100

extinction_prob_EN_possibilities <- seq(from=0.75, to=0.01, length.out=11)
extinction_prob_LC_possibilities <- seq(from=0.5, to=0.99, length.out=11)
proportion_EN_possibilites <- seq(from=-0, to=1, length.out=11)

all_parameters <- expand.grid(extinction_prob_EN = extinction_prob_EN_possibilities, extinction_prob_LC = extinction_prob_LC_possibilities, proportion_EN = proportion_EN_possibilites)

plot(all_parameters$extinction_prob_EN, all_parameters$extinction_prob_LC)


all_parameters$number_all_species <- number_all_species
all_parameters$total_time <- total_time


RunSimulation <- function (parameter_vector){
}


RunSimulation <- function (parameter_vector) {
  phy <- ape::rcoal(n=parameter_vector$number_all_species)
  taxa <- data.frame(taxon_name=phy$tip.label, endangered = FALSE, extant = TRUE)
  number_EN <- round(parameter_vector$number_all_species*parameter_vector$proportion_EN)
  print(ape::Ntip(phy))
  unlucky_taxa <- sample.int(n=parameter_vector$number_all_species, size=number_EN, replace=FALSE)
  taxa$endangered[unlucky_taxa] <- TRUE
  elapsed_time <- 0
  extinction_rate_EN <- -log(parameter_vector$extinction_prob_EN) / parameter_vector$total_time
  extinction_rate_LC <- -log(parameter_vector$extinction_prob_LC) / parameter_vector$total_time
  taxa$extinction_rate <- extinction_rate_LC
  taxa$extinction_rate[which(taxa$endangered==TRUE)] <- extinction_rate_EN
  taxa$extinction_time <- NA
  while(elapsed_time <- parameter_vector$total_time) {
    taxa_alive <- subset(taxa, extant==TRUE)
    wait_time <- rexp (1, rate=extinction_rate_EN*sum(taxa_alive$endangered) + extinction_rate_LC*sum(!taxa_alive))
    who_died <- sample(x=taxa_alive$taxon_name, size = 1, prob=taxa_alive$extinction_rate)
    taxa[which(taxa$taxaon_name==who_died)]$extant <- FALSE
    elapsed_time <- elapsed_time + wait_time
    taxa$extinction_time[which(taxa$taxon_name==who_died)] <- elapsed_time
    if(sum(taxa$extant)==0) {
      print(paste("Everything died. Nice job. Time =",elapsed_time))
      break()
    }
  }
  parameters$number_EN_extant <- sum(subset(taxa, endangered==TRUE)$extant)
  parameters$number_LC_extant <- sum(subset(taxa, endangered==FALSE)$extant)
  parameters$number_EN_extinct <- sum(!subset(taxa, endangered==TRUE)$extant)
  parameters$number_LC_extinct <- sum(!subset(taxa, endangered==FALSE)$extant)
  
  return(list(taxa=taxa, phy=phy, results=parameter_vector))
}


trial <- RunSimulation(parameter_vector=all_parameters[500,])


all_results <- data.frame()
for (parameter_row in sequence(nrow(all_parameters))) {
  print(paste("Now start", parameter_row, "of", nrow(all_parameters)))
  for (replicate_number in sequence(number_of_replicates)) {
    single_run <- RunSimulation(parameter_vector = all_parameters[parameter_row,])
    single_run_result <- single_run$result
    single_run_result$replicate_number <- replicate_number
    all_results <- rbind(all_results, single_run_result)
  }
}
