source('range_setup_asymm.R')

sigma = 0.8
f_sigma = .1
m_sigma = .5
mut_prob = 0
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

Pm = matrix(0,Nm,2) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,2) #probability of female preferences over time
Pf[,1] = f_init

preferences = array(0,dim=c(Nm,Nf,2))
female_tots = matrix(0,Nf,2)
pxy = array(0,dim=c(Nm,Nf,2))
pxy2 = array(0,dim=c(Nm,Nf,2))
growth = array(0,dim=c(Nm,2))
growth_desperate = array(0,dim=c(Nm,2))

t = 1
perc = 0
fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
fixed_weight = fixed_weight/sum(fixed_weight)
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

female_tots_analyt = step*dnorm(frange,mean=mmin,sd=sqrt(m_sigma^2+sigma^2))
Pm_analyt =step*dnorm(mrange,mean=(sigma^2*mmin+m_sigma^2*fmin)/(m_sigma^2+sigma^2),sd=sqrt(m_sigma^2*(m_sigma^2*f_sigma^2+(m_sigma^2+sigma^2)*sigma^2)/(m_sigma^2+sigma^2)^2))
growth_mean = (fmin*(m_sigma^2+sigma^2)-mmin*f_sigma^2)/(m_sigma^2+sigma^2-f_sigma^2)
growth_analyt = sqrt((m_sigma^2+sigma^2)^2)/sqrt(f_sigma^2*m_sigma^2+sigma^2*(m_sigma^2+sigma^2))*exp(-(mrange-growth_mean)^2/(2*(sigma^2+(f_sigma^2*(m_sigma^2+sigma^2)/(m_sigma^2+sigma^2-f_sigma^2)))))*exp(-(mmin-fmin)^2/(2*(f_sigma^2-(m_sigma^2+sigma^2))))

par(mfrow=c(1,3))
plot(female_tots_analyt,female_tots[,1]);abline(0,1,col='red')
plot(growth_analyt,growth[,1]);abline(0,1,col='red')
plot(log(Pm_analyt),log(apply(pxy[,,1],1,sum)));abline(0,1,col='red')