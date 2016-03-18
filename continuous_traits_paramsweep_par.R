## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)

# num_cores <- detectCores()-1
num_cores <-20
cl <-makeCluster(num_cores)
registerDoParallel(cl)

source('ind2sub.R')
source('glue.R')

## ---- integral -------------------------
int<-function(v){ #manual integration
	l=length(v)
	int_step/2*(sum(2*v)-v[1]-v[l])
} 

## ---- parameters ----------------
step = 0.01 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma2 = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# fmix_sigma2 = #variance of female distribution(s)
# mmix_sigma2 = 0.1 #variance of male distribution(s)
# Tsteps = #how many generations

mrange = seq(-10,10,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-10,10,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)

## ---- dynamics -----------------------

dynamics <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

# t = 1
for(t in 1:Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma2)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults) #normalization factor
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pm_beforemut = matrix(0,Nm)
	for(i in 1:Nm){
		Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	}
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
	Pm[,t+1] = Pm_aftermut
	Pf[,t+1] = Pf_adults
}
pop_dens = list(Pm=Pm,Pf=Pf)
return(pop_dens)
}

## ---- variance_sweep ---------------
mut_prob = 0.01
Tsteps = 1000
pm = 0.6
pf = 0.4

sigma2_vals = c(0.1,0.5,1,2,4)
Ns = length(sigma2_vals)
fmix_sigma2_vals = c(0.01,0.05,0.1,0.5,1)
Nfs = length(fmix_sigma2_vals)
mmix_sigma2_vals = c(0.01,0.05,0.1,0.5,1)
Nms = length(mmix_sigma2_vals)
P = Ns*Nfs*Nms
d = c(Ns,Nfs,Nms)

Pm_keep=as.list(1:P)
dim(Pm_keep)<-d
Pf_keep=as.list(1:P)
dim(Pf_keep)<-d

P_keep<-foreach(ind = 1:P, .combine='glue', .multicombine = TRUE, .init=list(list(),list())) %do% {
	v=ind2sub(d,ind)
	s=v[1]
	f=v[2]
	m=v[3]
	sigma2 = sigma2_vals[s]
	fmix_sigma2 = fmix_sigma2_vals[f]
	f_init = pf*dnorm(frange,fmin,fmix_sigma2)+(1-pf)*dnorm(frange,fmax,fmix_sigma2)
	mmix_sigma2 = mmix_sigma2_vals[m]
	m_init = pm*dnorm(mrange,mmin,mmix_sigma2)+(1-pm)*dnorm(mrange,mmax,mmix_sigma2)
	pop_dens = dynamics()
	pop_dens_last = list(pop_dens$Pm[,Tsteps],pop_dens$Pf[,Tsteps])
}

Pm_keep_hold = P_keep[[1]]
Pf_keep_hold = P_keep[[2]]

for(ind in 1:P){
	Pm_keep[[ind]] = Pm_keep_hold[[ind]]
	Pf_keep[[ind]] = Pf_keep_hold[[ind]]
}

save(Pm_keep=Pm_keep,Pf_keep=Pf_keep,file='/homes/ebrush/priv/song_learning_evolution/paramsweep_par.Rdata')

stopCluster(cl)

quit()




