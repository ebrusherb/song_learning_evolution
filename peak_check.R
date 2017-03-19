sigmay2_init = 4
sigmax2_init = 2
sigma2 = 0.01
mut_prob = 0
pf = 1
pm = 1
steps = 1000
store = 1

trait_chunk_num = 5
source("/Users/eleanorbrush/Documents/research/song_learning_evolution/range_setup.R")
source('discretize_and_dynamics.R')

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
continuous_weight[continuous_weight==0] = minweight

fixed_weight = continuous_weight/sum(continuous_weight)

pref_chunk_num = 15
fixed_weight = discretize(continuous_weight)
fixed_weight = fixed_weight / sum(fixed_weight)

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
j = midpt
# j = min(which(mrange>0))
s = sign(midpt-j+0.5)
toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
hold = f_init
hold[toreplace[1]:toreplace[2]] = f_init[pull[1]:pull[2]]
if(toreplace[1]>1){
	hold[1:(toreplace[1]-1)] = f_init[(pull[2]+1):Nf]}
f_init = hold

m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

##### good:
alpha = 1
fixed_weight = rep(1,trait_chunk_num)
fixed_weight[midpt] = 1+alpha
fixed_weight = fixed_weight / sum(fixed_weight)
f_init = c(0.1,0.2,0.4,0.2,0.1)
m_init = c(0.4,0.2,0.2,0.1,0.1)

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
	Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
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