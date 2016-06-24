source('range_setup_symm.R')

sigma = 0.1
f_sigma = 1
m_sigma = .01
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

female_tots_analyt = dnorm(frange,mean=fmin,sd=sqrt(m_sigma^2+sigma^2))
i = 650
j = 660; 
x = mrange[i]
y = frange[j]
pxy_analyt = step^2*1/(2*pi)*(sqrt(m_sigma^2+sigma^2))/(f_sigma*m_sigma*sigma)*exp(-1/2*(1/f_sigma^2-1/(m_sigma^2+sigma^2))*(1+y)^2)*exp(-(1+x)^2/(2*m_sigma^2))*exp(-(x-y)^2/(2*sigma^2))
# growth_analyt = step*1/sqrt(2*pi)*sqrt(m_sigma^2+sigma^2)/(f_sigma*sigma)*exp(-1/2*(1/f_sigma^2-1/(m_sigma^2+sigma^2))*(1+y)^2)*exp(-(x-y)^2/(2*sigma^2))
# growth_analyt_range =  step*1/sqrt(2*pi)*sqrt(m_sigma^2+sigma^2)/(f_sigma*sigma)*exp(-1/2*(1/f_sigma^2-1/(m_sigma^2+sigma^2))*(1+frange)^2)*exp(-(x-frange)^2/(2*sigma^2))
growth_analyt = sqrt(m_sigma^2*sigma^2*(m_sigma^2+sigma^2))/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(1+x)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))

# x = -2 
# m_sigma = .001;sigma = .001;f_sigma = 1;v =(-(x-frange)^2/(2*sigma^2))+(-1/2*(1/f_sigma^2-1/(m_sigma^2+sigma^2))*(1+frange)^2);plot(v);1/f_sigma^2-1/(m_sigma^2+sigma^2);sigma_new = 1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))*sigma^2/((1/f_sigma^2-1/(m_sigma^2+sigma^2))+sigma^2);sigma_new

i=270
x=mrange[i]
j = 500
y = frange[j]
v=Pf[,t]*preferences[i,,t]/female_tots[,t]
check = 1/sqrt(2*pi*f_sigma^2)*exp(-(fmin-y)^2/(2*f_sigma^2))*1/sqrt(2*pi*sigma^2)*exp(-(x-y)^2/(2*sigma^2))/(1/sqrt(2*pi*(m_sigma^2+sigma^2))*exp(-(fmin-y)^2/(2*(m_sigma^2+sigma^2))))
check = sqrt(m_sigma^2+sigma^2)/sqrt(2*pi*f_sigma^2*sigma^2)*exp(-1/2*(1/f_sigma^2-1/(m_sigma^2+sigma^2))*(fmin-y)^2)*exp(-1/2/sigma^2*(x-y)^2)
mua=fmin
mub=x
sigma2a=1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))
sigma2b=sigma^2
muab=(mua*sigma^2+mub/(1/f_sigma^2-1/(m_sigma^2+sigma^2)))/(1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))+sigma^2)
sigma2ab=sigma^2*f_sigma^2*(m_sigma^2+sigma^2)/(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)
check = sqrt(m_sigma^2+sigma^2)/sqrt(2*pi*f_sigma^2*sigma^2)*exp(-(y-mua)^2/(2*sigma2a))*exp(-(y-mub)^2/(2*sigma2b));print(check)
check = sqrt(m_sigma^2+sigma^2)/sqrt(2*pi*f_sigma^2*sigma^2)*exp(-(y-muab)^2/(2*sigma2ab))*exp(-(mua-mub)^2/(2*(sigma2a+sigma2b)));print(check)
check = sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*1/sqrt(2*pi)*exp(-(y-muab)^2/(2*sigma2ab))*exp(-(mua-mub)^2/(2*(sigma2a+sigma2b)));print(check)
check = sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*1/sqrt(2*pi)*exp(-1/2*(y-muab)^2/sigma2ab)*exp(-1/2*(x-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))));print(check)
check = sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*sqrt(sigma2ab)*exp(-1/2*(x-mmin)^2/(sigma2a+sigma2b))*1/sqrt(2*pi*sigma2ab)*exp(-1/2*(y-muab)^2/sigma2ab);print(check)
check_vec = sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*sqrt(sigma2ab)*exp(-1/2*(x-mmin)^2/(sigma2a+sigma2b))*1/sqrt(2*pi*sigma2ab)*exp(-1/2*(frange-muab)^2/sigma2ab)
print(c(v[j],step*check))
print(c(int(v),step*sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*sqrt(sigma2ab)*exp(-1/2*(mmin-x)^2/(sigma2a+sigma2b))))
print(c(int(v),step*sqrt(m_sigma^2+sigma^2)/sqrt(f_sigma^2*sigma^2)*sqrt(f_sigma^2*sigma^2*(m_sigma^2+sigma^2))/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(x-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))))
print(c(int(v),step*(m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(x-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))))
growth_analyt = (m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))
growth_analyt = (m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2)*exp(-1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))
growth_analyt_log = log((m_sigma^2+sigma^2)/sqrt(sigma^2*(m_sigma^2+sigma^2)+f_sigma^2*m_sigma^2))-(1/2*(mrange-mmin)^2/(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2))))
# plot(v,step*check_vec);abline(0,1)
print(sigma^2+1/(1/f_sigma^2-1/(m_sigma^2+sigma^2)))
plot(log(growth[,t]),log(growth_analyt));abline(0,1,col='red')

m = m_sigma
for(T in 2:1000){
	mnew = m[T-1]*(sigma^4+(sigma^2+f_sigma^2)*m[T-1])/(m[T-1]+sigma^2)^2
	m = c(m,mnew)
}
m2 = m_sigma
for(T in 2:1000){
	mnew = m2[T-1]^4*(f_sigma^2+(m2[T-1]+sigma^2)*m2[T-1]*sigma^2)/(m2[T-1]+sigma^2)^2
	m2 = c(m2,mnew)
}
# plot((m));abline(h=(f_sigma^2-sigma^2),col='red')

# m_sigma_predicted = array(0,dim(Pm_onepop))

# for(i in 1:prod(dim(Pm_onepop))){
	# sub = ind2sub(dim(Pm_onepop),i)
	# sigma = sigma_vals[sub[1]]
	# f_sigma = f_sigma_vals[sub[2]]
	# m_sigma = m_sigma_vals[sub[3]]
	# m = m_sigma
	# for(T in 2:5000){
		# m = c(m,m[T-1]*(sigma^4+(sigma^2+f_sigma^2)*m[T-1])/(m[T-1]+sigma^2)^2)
	# }
	# m_sigma_predicted[i] = m[length(m)]
# }

# m_sigma_predicted_malefath_femalefath = array(0,dim(Pm_onepop))

# for(i in 1:prod(dim(Pm_onepop))){
	# sub = ind2sub(dim(Pm_onepop),i)
	# sigma = sigma_vals[sub[1]]
	# f_sigma = f_sigma_vals[sub[2]]
	# m_sigma = m_sigma_vals[sub[3]]
	# m = m_sigma
	# for(T in 2:5000){
		# m = c(m,m[T-1]*(sigma^4+(sigma^2+m[T-1])*m[T-1])/(m[T-1]+sigma^2)^2)
	# }
	# m_sigma_predicted_malefath_femalefath[i] = m[length(m)]
# }