source('ind2sub.R')
source('glue.R')
source('int.R')
source('range_setup_symm.R')

## ---- dynamics -----------------------

dynamics <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)	
	fixed_weight[fixed_weight==0] = minweight
	fixed_weight = fixed_weight/sum(fixed_weight)
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
	Pm_aftermut[which(Pm_aftermut<1e-100)] = 0
	Pf_adults[which(Pf_adults<1e-100)] = 0
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


Tsteps = 100
store = 1
Tmin = Tsteps

sigma = 0.1
f_sigma = 1
m_sigma = .01
mut_prob = 0.01
pf = 1
pm = 1
f_init = dnorm(frange,fmin,f_sigma)
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,m_sigma)
# m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

pop_dens = dynamics()
Pm = pop_dens$Pm
Pf = pop_dens$Pf

preferences = array(0,dim=c(Nm,Nf,Tmin))
female_tots = matrix(0,Nf,Tmin)
pxy = array(0,dim=c(Nm,Nf,Tmin))
growth = array(0,dim=c(Nm,Tmin))
growth_desperate = array(0,dim=c(Nm,Tmin))
growth_analyt = array(0,dim=c(Nm,Tmin))
var_m = vector(length=Tmin)
m = m_sigma^2

for(t in 1:Tmin){
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
		z = sum(weight*Pm[,t])
		preferences[,j,t] = weight #preference given by each female to each male
		
		female_tots[j,t] = z #total preferences given by each female
		if(z>z_thresh){
			pxy[,j,t] = Pf[j,t]*weight*Pm[,t]/z
			}
	}
	zero = which(Pm[,t]==0)
	growth[,t] =  apply(pxy[,,t],1,sum)/Pm[,t]
	growth[zero,t] = 0
	growth_desperate[,t] = apply(preferences[,,t]/matrix(rep(female_tots[,t],times=Nm),nrow=Nm,byrow=TRUE),1,sum)
	growth_desperate[zero,t] = 0
	growth_analyt[,t] = (m[t]+sigma^2)/sqrt(sigma^2*(m[t]+sigma^2)+f_sigma^2*m[t])*exp(-1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m[t]+sigma^2))))
	m = c(m,m[t]*(sigma^4+(sigma^2+f_sigma^2)*m[t])/(m[t]+sigma^2)^2)
	var_m[t] = sum((mrange-(mmin))^2*Pm[,t])
}

Pm_analyt = array(0,dim=c(Nm,Tmin))
Pm_analyt[,1] = Pm[,1]
for(t in 2:Tmin){
	Pm_analyt[,t] = growth_analyt[,t-1]*Pm_analyt[,t-1]
}

female_tots_analyt = dnorm(frange,mean=fmin,sd=sqrt(m_sigma^2+sigma^2))
# growth_analyt = (m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))
growth_analyt_log = log((m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2))-(1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))
