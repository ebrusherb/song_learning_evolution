#yes when the initial song distribution is overdispersed

sigmay2_init = 1
sigmax2_init = 1e-4
sigmax2_init = 5e-3
sigma2 = 0.2
mut_prob = 0
pf = 1
pm = 1
steps = 20
store = 1

trait_chunk_num = 351
source("/Users/eleanorbrush/Documents/research/song_learning_evolution/range_setup.R")
source('/Users/eleanorbrush/Desktop/song_learning_evolution-af45d67bc0ab76fabc6d9e6af02ec1735362ef29/range_setup_old.R')
source('/Users/eleanorbrush/Documents/research/song_learning_evolution/discretize_and_dynamics.R')

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
minweight = 10^(-320)
continuous_weight[continuous_weight==0] = minweight

fixed_weight = continuous_weight/sum(continuous_weight)

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init[f_init<1e-8] = 1e-8
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)

m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
# m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init[m_init<1e-8] = 1e-8
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	# Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_adults
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_adults
			t = steps+1
			}
}

plot(log(Pm[,1]));points(log(Pm[,steps]),col='red')

# ##########
# rm(list=ls())
# trait_chunk_num = 1501
# source("/Users/eleanorbrush/Documents/research/song_learning_evolution/range_setup.R")
# l1 = lapply(ls(),function(x) get(x))
# names(l1) = setdiff(ls(),c(ls(pattern='l1')))
# rm(list=setdiff(ls(),c(ls(pattern='l1'),ls(pattern='trait_chunk_num'))))
# source('/Users/eleanorbrush/Desktop/song_learning_evolution-af45d67bc0ab76fabc6d9e6af02ec1735362ef29/range_setup_old.R')
# l2 = lapply(setdiff(ls(),c(ls(pattern='l1'))),function(x) get(x))
# names(l2) = setdiff(ls(),c(ls(pattern='l1'),ls(pattern='l2')))