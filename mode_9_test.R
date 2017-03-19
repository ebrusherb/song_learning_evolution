trait_chunk_num = 151
sigmay2 = 3
sigmay2_init = sigmay2
sigmax2 = 0.8
sigmax2_init = sigmax2
pf = 1
pm = 1 
rho = 0.6
minweight = 10^(-320)
mut_prob = 0

steps = 25
store = 1
source('range_setup.R')


f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

sigma2 = 0.5
continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = array(0,c(Nm,Nf,steps+1)) #probability of female preferences over time
Pf[,,1] = sapply(matrix(f_init,nrow=1),function(x) x*m_init,simplify=TRUE)

reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,,t]
	pxy = array(0,c(Nm,Nm,Nf)) #probability of a x / (x,y) pair
	Pf_new = array(0,c(Nm,Nf))
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
		}
		song_pairs = unlist(lapply(split(pxy[,,j],reverse_diagonal),sum))
		song_pairs = song_pairs[seq(1,2*Nm-1,by=2)]+1/2*c(song_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,song_pairs[seq(2,2*Nm-1,by=2)])
		Pf_new[,j] = song_pairs
	}
	### by plotting apply(pxy,1,sum) against what i know variance of successful mating males is I see there's something wrong with pxy
	Pf_adults[,apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_new = apply(Pf_new,1,sum)/sum(Pf_adults)
	Pf_new = Pf_new/sum(Pf_adults) 	
	nonzero = which(Pm_new>nonzero_thresh)
	perc = max(abs(range(Pm_new[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_new
		Pf[,,t+1] = Pf_new
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_new
			Pf[,,(t+1):(steps+1)] = Pf_new
			print(c(t,'problem'))
			t = steps+1
			}
}