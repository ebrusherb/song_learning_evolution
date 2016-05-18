## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)
library(abind)

# num_cores <- detectCores()-1
num_cores <-20
cl <-makeCluster(num_cores)
registerDoParallel(cl)

source('ind2sub.R')
source('glue.R')
source('int.R')


## ---- parameters ----------------
step = 0.05 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma2 = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# fmix_sigma2 = #variance of female distribution(s)
# mmix_sigma2 = 0.1 #variance of male distribution(s)
# Tsteps = #how many generations

mrange = seq(-7.5,7.5,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-7.5,7.5,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)

nonzero_thresh = 1e-5
perc_thresh = 1e-10

## ---- dynamics -----------------------

dynamics_malefath_femaleboth <-function(){

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

Pm = array(0,dim=c(Nm,Nf,Tsteps+1)) #probability of male songs and preferences over time
Pm[,,1] = outer(m_init,f_init,'*')

t = 1
perc = 0
while(t <= Tsteps){
	Pm_adults = Pm[,,t]
	m_adult_songs = apply(Pm_adults,1,int)
	divisor = m_adult_songs
	divisor[m_adult_songs==0] = 1
	m_adults_prefs_normed = Pm_adults / matrix(rep(divisor,times = Nf),nrow=Nm)
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	for(j in 1:Nf){
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma2)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*m_adult_songs) #normalization factor
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*m_adult_songs/z
			}
	}
	# m_prefs.rows = split(m_adults_prefs_normed,row(m_adults_prefs_normed))
	# pxy.rows = split(pxy,row(pxy))
	# total_breeding_mat = lapply(seq_along(pxy.rows),function(i){t(outer(pxy.rows[[i]],m_prefs.rows[[i]],'/'))})
	total_breeding_mat = array(0,dim=c(Nm,Nf,Nf))
	for(i in 1:Nm){
		total_breeding_mat[i,,] = t(outer(pxy[i,],m_adults_prefs_normed[i,],'*'))
	}
	
	Pm_beforemut = array(0,dim=c(Nm,Nf))
	for(i in 1:Nm){
		Pm_beforemut[i,] = 1/2*apply(total_breeding_mat[i,,],1,int)+1/2*apply(total_breeding_mat[i,,],2,int)
	}
	Pm_aftermut = array(0,dim=c(Nm,Nf))
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*abind(Pm_beforemut[2:Nm,],array(0,Nf),along=1) +		mut_prob/2*abind(array(0,Nf),Pm_beforemut[1:Nm-1,],along=1) #and then they change their songs
	pref_breeding_mat = apply(total_breeding_mat,c(2,3),int)
	Pf_new = 1/2*apply(pref_breeding_mat,1,int) + 1/2*apply(pref_breeding_mat,2,int)
	nonzero = which(m_adult_songs>nonzero_thresh)
	perc = max(abs(range(apply(Pm_aftermut[nonzero,],1,int)/m_adult_songs[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_new
		t = t+1
		} else{
			Pm[,,(t+1):(Tsteps+1)] = Pm_aftermut
			Pf[,,(t+1):(Tsteps+1)] = Pf_new
			t = Tsteps+1
			}
}
Pm = apply(Pm,c(1,3),int)
pop_dens = list(Pm=Pm[,(Tsteps-200):Tsteps],Pf=Pf[,(Tsteps-200):Tsteps])
return(pop_dens)
}

## ---- variance_sweep ---------------
Tsteps = 10000
pm = 0.6
pf = 0.6

sigma2_vals = c(0.001,0.01,0.1,1)
Ns = length(sigma2_vals)
fmix_sigma2_vals = c(0.01,0.1,1)
Nfs = length(fmix_sigma2_vals)
mmix_sigma2_vals = c(0.01,0.1,1)
Nms = length(mmix_sigma2_vals)
mut_prob_vals = c(0,0.01,0.1)
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
	# sigma2 = sigma2_vals[s]
	# fmix_sigma2 = fmix_sigma2_vals[f]
	# f_init = pf*dnorm(frange,fmin,fmix_sigma2)+(1-pf)*dnorm(frange,fmax,fmix_sigma2)
	# mmix_sigma2 = mmix_sigma2_vals[m]
	# m_init = pm*dnorm(mrange,mmin,mmix_sigma2)+(1-pm)*dnorm(mrange,mmax,mmix_sigma2)
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
	sigma2 = sigma2_vals[s]
	fmix_sigma2 = fmix_sigma2_vals[f]
	mmix_sigma2 = mmix_sigma2_vals[m]
	f_init = pf*dnorm(frange,fmin,fmix_sigma2)+(1-pf)*dnorm(frange,fmax,fmix_sigma2)
	m_init = pm*dnorm(mrange,mmin,mmix_sigma2)+(1-pm)*dnorm(mrange,mmax,mmix_sigma2)
	mut_prob = mut_prob_vals[p]
	pop_dens = dynamics_malefath_femaleboth()
	# pop_dens_last = list(pop_dens$Pm[,Tsteps],pop_dens$Pf[,Tsteps])
}

Pm_onepop_hold = P_onepop[[1]]
Pf_onepop_hold = P_onepop[[2]]

for(ind in 1:P){
	Pm_onepop[[ind]] = Pm_onepop_hold[[ind]]
	Pf_onepop[[ind]] = Pf_onepop_hold[[ind]]
}

Date=Sys.Date()
# save(Pm_keep=Pm_keep,Pf_keep=Pf_keep,Pm_onepop=Pm_onepop,Pf_onopop=Pf_onepop,sigma2_vals=sigma2_vals,fmix_sigma2_vals,mmix_sigma2_vals,mut_prob_vals=mut_prob_vals,file='/homes/ebrush/priv/song_learning_evolution/song_learning_paramsweep_par.Rdata')
save(Pm_onepop=Pm_onepop,Pf_onepop=Pf_onepop,sigma2_vals=sigma2_vals,fmix_sigma2_vals,mmix_sigma2_vals,mut_prob_vals=mut_prob_vals,Tsteps=Tsteps,mrange,Nm,frange,Nf,step,int_step,file=paste('/homes/ebrush/priv/song_learning_evolution/song_learning_malefath_femaleboth_',substr(Date,1,4),'_',substr(Date,6,7),'_',substr(Date,9,10),'.Rdata',sep=''))


stopCluster(cl)

quit()




