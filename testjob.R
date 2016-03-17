library(parallel)
library(foreach)
library(doParallel)

# num_cores <- detectCores()-1
num_cores <-2
cl <-makeCluster(num_cores)
registerDoParallel(cl)

sigma2_vals = c(0.1,0.5,1,2,4)
Ns = length(sigma2_vals)
fmix_sigma2_vals = c(0.01,0.05,0.1,0.5,1)
Nfs = length(fmix_sigma2_vals)
mmix_sigma2_vals = c(0.01,0.05,0.1,0.5,1)
Nms = length(mmix_sigma2_vals)
P = Ns*Nfs*Nms

M=as.list(1:P)
dim(M)<-c(Ns,Nfs,Nms)

ind2sub <-function(v,ind){
	# v=c(2,3)
	# ind=5
	num_dim = length(v)
	sub <- array(0,dim=c(1,num_dim))
	P = cumprod(v)
	for(l in num_dim:2){ #work backwards to fill in from the last dimension
		ind = (ind-1)%%P[l]+1 # index in current dimension
		div = P[l-1] #size of current dimension
		sub[l] = floor((ind-1)/div)+1 #number of times ind-1 goes into product of higher dimensions + 1
		}
	sub[1] = (ind-1)%%div+1 #left over after ind-1 goes into number or rows
	return(sub)
}

invisible(
foreach(ind = 1:P) %do% {
	sub = ind2sub(dim(M),ind);
	M[[ind]] = sub;
}
)

ind=sample(1:P,1)
v=ind2sub(dim(M),ind)
cbind(as.vector(v),as.vector(M[[v[1],v[2],v[3]]]))

stopCluster(cl)

quit()
