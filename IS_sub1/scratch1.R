
theta = 100 #100
R = 20 #20
tf = 0.1 #0.1
target = paste("genetree/gtree_file",1,sep="")
is.iter = 1
isss = 10000


#write(c(log(theta),log(R),log(tf)),"paramfile",ncol=3)
write(cbind(log(theta),log(R),log(tf)),"paramfile",ncol=3)
cm1 = paste("./is_moments",target,as.integer(is.iter),as.integer(isss),"paramfile")
wvec = system(cm1,intern=T)
wvec = as.numeric(wvec)
print(wvec)

