setwd('/Users/eleanorbrush/Desktop/song_learning_evolution-af45d67bc0ab76fabc6d9e6af02ec1735362ef29/')
source('int.R')
source('range_setup_old.R')

## ---- dynamics -----------------------

dynamics_old <-function(){
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
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
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
	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pf_adults[apply(pxy,2,sum)==0]=0
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_adults
		t = t+1
		} else{
			Pm[,(t+1):(Tsteps+1)] = Pm_aftermut
			Pf[,(t+1):(Tsteps+1)] = Pf_adults
			t = Tsteps+1
			}
}
pop_dens = list(Pm=Pm[,(store):Tsteps],Pf=Pf[,(store):Tsteps])
return(pop_dens)
}

Tsteps = 20
store = 1
pf = 1
pm = 1
mut_prob = 0

sigma = 0.4
f_sigma = 1
m_sigma = 0.01

f_init = dnorm(frange,fmin,f_sigma)
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,m_sigma)
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)
fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)	
fixed_weight[fixed_weight==0] = minweight
fixed_weight = fixed_weight/sum(fixed_weight)

p1 = dynamics_old()

sigma2 = sigma^2
sigmay2 = f_sigma^2
sigmax2 = m_sigma^2
steps = Tsteps

trait_chunk_num = Nm
# source("/Users/eleanorbrush/Documents/research/song_learning_evolution/range_setup.R")
source('/Users/eleanorbrush/Documents/research/song_learning_evolution/discretize_and_dynamics.R')

p2 = dynamics()

source('/Users/eleanorbrush/Documents/research/song_learning_evolution/recursion_all.R')
r = recursion_all(sigmax2,sigmay2,sigma2,0)