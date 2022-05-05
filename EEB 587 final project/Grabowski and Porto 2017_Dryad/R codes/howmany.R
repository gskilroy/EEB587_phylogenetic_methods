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