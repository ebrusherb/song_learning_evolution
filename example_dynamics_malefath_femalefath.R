## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)

# num_cores <- detectCores()-1
# num_cores <-10
# cl <-makeCluster(num_cores)
# registerDoParallel(cl)

source('ind2sub.R')
source('glue.R')
source('int.R')
source('range_setup.R')

## ---- dynamics -----------------------

dynamics_malefath_femalefath <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

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
	song_beforemut = apply(pxy,1,sum)/sum(Pf_adults)	
	song_aftermut = (1-mut_prob)*song_beforemut + mut_prob/2*c(song_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,song_beforemut[1:Nm-1]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(song_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,t+1] = song_aftermut
		Pf[,t+1] = song_beforemut
		t = t+1
		} else{
			Pm[,(t+1):(Tsteps+1)] = song_aftermut
			Pf[,(t+1):(Tsteps+1)] = song_beforemut
			t = Tsteps+1
			}
}
pop_dens = list(Pm=Pm[,(store):Tsteps],Pf=Pf[,(store):Tsteps])
return(pop_dens)
}

Tsteps = 5000
store = 1
pf = 1
pm = 1
store = 1
# s = 1
# f = Nfs
# m = 1
# p = 2
# sigma = sigma_vals[s]
# f_sigma = f_sigma_vals[f]
# m_sigma = m_sigma_vals[m]
# mut_prob = mut_prob_vals[p]
sigma = 0.1
f_sigma = 1
m_sigma = 0.5
mut_prob = 0.01
f_init = dnorm(frange,fmin,f_sigma)
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,m_sigma)
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

pop_dens = dynamics_malefath_femalefath()

Date=Sys.Date()
save(sigma,f_sigma,m_sigma,mut_prob,pop_dens,file=paste('/homes/ebrush/priv/song_learning_evolution/song_learning_malefath_femalefath_example_',substr(Date,1,4),'_',substr(Date,6,7),'_',substr(Date,9,10),'.Rdata',sep=''))

