num_monomorphic = 0
sstable = matrix(ncol=16,nrow=nsim_loci)
for(j in 1:nsim_loci){

	num_attempts = 1
	repeat{
		# perform the ms command
		system(paste("./ms 20",1," -t", theta, "-G",round(alpha,5),"-eG",tf, 0, "> testdata1"))
		
		#attempt to read the ms output
		data = try(read.table(file="testdata1",colClasses="factor",skip=6))
		if(class(data) != "try-error")break
		if(num_attempts > 10)stop("gen_test_data: too many monomorphic loci...")
		num_attempts = num_attempts + 1
	}
	num_monomorphic = num_monomorphic + num_attempts - 1
	#convert to a genetree input file
	freq = as.data.frame(table(data))
	n = length(freq[,2])
	mydata = data.frame(rep(0,n),freq[,2],rep(":",n),freq[,1])
	sink(paste("genetree/gtree_file",j,sep=""),append=F,split=F)
	write.table(format(mydata,justify="right"),row.names=F,col.names=F,quote=F)
	sink()
	
	#write files for ABC
	system("./sample_stat_mab3U < testdata1 > res1") #unfolded fs
	res1 = scan("res1")
	sstable[j,] = res1[c(1:15,18)]
	#above we drop the 16th and 17th summary stats (median frequency of haplotypes and frequency of rarest haplotype)
	#because these both tend to be uniformly 1 and not very informative once we have long haplotypes. 
	
	#also, let's keep ms output appended in case we want to compute more stats.
	if(j == 1){
		fname = paste("msoutput_",mainiter,sep="")
		cmd = paste("mv testdata1",fname)
	
		system(cmd)
	}else{
		
		system(paste("tail -n +3 testdata1 >>",fname))
	}
}
write(t(sstable),file=paste("abcfiles/target_ss_",mainiter,".txt",sep=""),ncol=16)
system(paste("mv ",fname,"abcfiles/"))
if(num_monomorphic > 0)print(paste("WARNING number of monomorphic loci is: ",num_monomorphic))

