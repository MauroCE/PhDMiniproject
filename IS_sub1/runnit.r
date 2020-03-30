##First read in the arguments listed at the command line
args=commandArgs()
#this will work for older versions of R as well (newer versions allow argument TRUE, and only look for arguments after --args)

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
check = pmatch("--args",args)
if(is.na(check) | check == length(args)){
    print("No arguments supplied.")
    stop()##supply default values

}else{
    print(paste("number of arguments:",length(args)-check))
    print(args[check+1])
}

# define the number of loci to be simulated by ms
nsim_loci <- 100
#NB don't confuse this with Nloci, which is the number of  loci analysed

# set the log ratio current size to ancestral size
# sample from prior
# alpha in ms for the -G parameter 
# is related to my parameterisation as log(R)/tf

#!!!!!!!!!!!!!!!!!NBBB - Because ms uses silly scaling for time (scales by 4N rather than 2N)
#!!!!!!!!!!!! We need to make sure that we input the correct value to growthtest
#!!!!!!!!!!!! e.g. data simulated with tf = 1 for ms, will be inferred to have tf=2 in growthtest 
# so to be able to get EP-ABC and EP-IS to match up we need to set a flag for growthtest so that 
# it multiplies input tf values by 2 internally. E.g. if we want the MAP to be at tf=1 when we simulated with tf = 1 in ms, we want 
#it such that a value of 1 input into growthtest comes back with the MAP density, so it needs to be multiplied by 2. 

#where R =N_c/N_A
#and tf = T/4N_c
#theta = 4N_cmu
#so lets have prior for log(R), log(tf), log(theta) ~ N((0,0,0),diag(1,1,1))
#imagine mu is 10^-8 per base pair and we have 10000 bases.
#so theta = 100 is equivalent to a pop size of 100/4e.-4 = 250000
#R = 20 implies Na = 12500, and tf = 0.1 gives T= 1.0e5 generations
theta = 20
R = 20
tf = 0.05

logtheta = log(theta)
logR = log(R)
logtf = log(tf)

alpha = logR/exp(logtf)
library(mnormt)
source("is_routine2.R") #this one uses IS-corrected widened prior
report_counter = 1 #this is just to keep a track of reporting evolution of posterior (in ep_IS.R)
system("rm -f abcfiles/*")
for(mainiter in 1:20){
	system("rm -f genetree/*")
	source("gen_test_data.r") 
	source("ep_IS3a.R")
	if(mainiter == 1){
		write(c(log(theta),log(R),log(tf),final.out),file="results.txt",ncol=length(final.out)+3)
	}else{
		write(c(log(theta),log(R),log(tf),final.out),file="results.txt",ncol=length(final.out)+3,append=T)
	}


}
