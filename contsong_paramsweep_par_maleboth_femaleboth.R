## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)
library(abind)

# num_cores <- detectCores()-1
num_cores <-10
cl <-makeCluster(num_cores)
registerDoParallel(cl)

source('ind2sub.R')
source('glue.R')
source('int.R')
source('range_setup.R')



## ---- dynamics -----------------------

dynamics_maleboth_femaleboth <-function(){

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

Pm = array(0,dim=c(Nm,Tsteps+1)) #probability of male songs and preferences over time
Pm[,1] = m_init

t = 1
perc = 0
while(t <= Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)	
	fixed_weight[fixed_weight==0] = minweight
	fixed_weight = fixed_weight/sum(fixed_weight)
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pf_adults[apply(pxy,2,sum)==0]=0
	Pm_beforemut = 1/2*(apply(pxy,1,sum)+apply(pxy,2,sum))/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	Pf_new = Pm_beforemut	
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_new
		t = t+1
		} else{
			Pm[,(t+1):(Tsteps+1)] = Pm_aftermut
			Pf[,(t+1):(Tsteps+1)] = Pf_new
			t = Tsteps+1
			}
}
pop_dens = list(Pm=Pm[,store:Tsteps],Pf=Pf[,store:Tsteps])
return(pop_dens)
}

## ---- variance_sweep ---------------
Tsteps = 5000
store = Tsteps-200
pm = 0.6
pf = 0.6

sigma_vals = c(0.1,0.2,0.4,0.8,1,1.5)
Ns = length(sigma_vals)
f_sigma_vals = c(0.1,0.25,0.5,0.75,1)
Nfs = length(f_sigma_vals)
m_sigma_vals = c(0.1,0.5)
Nms = length(m_sigma_vals)
mut_prob_vals = c(0.01)
Nmp = length(mut_prob_vals)
P = Ns*Nfs*Nms*Nmp
d = c(Ns,Nfs,Nms,Nmp)

# # Pm_keep=as.list(1:P)
# dim(Pm_keep)<-d
# Pf_keep=as.list(1:P)
# dim(Pf_keep)<-d

# P_keep<-foreach(ind = 1:P, .combine='glue', .multicombine = TRUE, .init=list(list(),list())) %dopar% {
	# v=ind2sub(d,ind)
	# s=v[1]
	# f=v[2]
	# m=v[3]
	# p=v[4]
	# sigma = sigma_vals[s]
	# f_sigma = f_sigma_vals[f]
	# f_init = pf*dnorm(frange,fmin,f_sigma)+(1-pf)*dnorm(frange,fmax,f_sigma)
	# m_sigma = m_sigma_vals[m]
	# m_init = pm*dnorm(mrange,mmin,m_sigma)+(1-pm)*dnorm(mrange,mmax,m_sigma)
	# mut_prob = mut_prob_vals[p]
	# pop_dens = dynamics_malefath_femaleboth()
	# # pop_dens_last = list(pop_dens$Pm[,Tsteps],pop_dens$Pf[,Tsteps])
# }

# Pm_keep_hold = P_keep[[1]]
# Pf_keep_hold = P_keep[[2]]

# for(ind in 1:P){
	# Pm_keep[[ind]] = Pm_keep_hold[[ind]]
	# Pf_keep[[ind]] = Pf_keep_hold[[ind]]
# }

Pm_onepop=as.list(1:P)
dim(Pm_onepop)<-d
Pf_onepop=as.list(1:P)
dim(Pf_onepop)<-d
pf = 1
pm = 1

P_onepop<-foreach(ind = 1:P, .combine='glue', .multicombine = TRUE, .init=list(list(),list())) %dopar% {
	v=ind2sub(d,ind)
	s=v[1]
	f=v[2]
	m=v[3]
	p=v[4]
	sigma = sigma_vals[s]
	f_sigma = f_sigma_vals[f]
	m_sigma = m_sigma_vals[m]
	f_init = dnorm(frange,fmin,f_sigma)
	f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
	f_init = f_init/sum(f_init)
	f_init = pf*f_init+(1-pf)*rev(f_init)
	m_init = dnorm(mrange,mmin,m_sigma)
	m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
	m_init = m_init / sum(m_init)
	m_init = pf*m_init+(1-pf)*rev(m_init)
	mut_prob = mut_prob_vals[p]
	pop_dens = dynamics_maleboth_femaleboth()
	# pop_dens_last = list(pop_dens$Pm[,Tsteps],pop_dens$Pf[,Tsteps])
}

Pm_onepop_hold = P_onepop[[1]]
Pf_onepop_hold = P_onepop[[2]]

for(ind in 1:P){
	Pm_onepop[[ind]] = Pm_onepop_hold[[ind]]
	Pf_onepop[[ind]] = Pf_onepop_hold[[ind]]
}

Date=Sys.Date()
# save(Pm_keep=Pm_keep,Pf_keep=Pf_keep,Pm_onepop=Pm_onepop,Pf_onopop=Pf_onepop,sigma_vals=sigma_vals,f_sigma_vals,m_sigma_vals,mut_prob_vals=mut_prob_vals,file='/homes/ebrush/priv/song_learning_evolution/song_learning_paramsweep_par.Rdata')
save(Pm_onepop=Pm_onepop,Pf_onepop=Pf_onepop,sigma_vals=sigma_vals,f_sigma_vals,m_sigma_vals,mut_prob_vals=mut_prob_vals,Tsteps=Tsteps,mrange,Nm,frange,Nf,step,int_step,file=paste('/homes/ebrush/priv/song_learning_evolution/song_learning_maleboth_femaleboth_',substr(Date,1,4),'_',substr(Date,6,7),'_',substr(Date,9,10),'.Rdata',sep=''))


stopCluster(cl)

quit()




