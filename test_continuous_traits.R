Tsteps = 100
mut_delta = 0 #how to implement mutations of different sizes?
step = 0.01
int_step = step
sigma2 = .5
mut_prob = 0.1
alpha = 0.5
mix_sigma2 = 0.01

#########
int<-function(v){
	l=length(v)
	int_step/2*(sum(2*v)-v[1]-v[l])
}


mrange = seq(-5,6,by=step)
Nm = length(mrange)
mrange_orig = seq(-1,1,by=step)
frange = seq(-5,6,by=step)
Nf = length(frange)
frange_orig = seq(-1,1,by=step)

m0 = which(mrange==-1)
m1 = which(mrange==1)
Pm = matrix(0,Nm,Tsteps+1)
breaks = runif(n=length(mrange_orig)-1)
breaks = sort(breaks)
m_init = c(breaks[1],diff(c(breaks,1)))
m_init = runif(n=length(mrange_orig))
m_init = m_init/int(m_init)
Pm[m0:m1,1] = m_init
# Pm[m0] = 0.6
# Pm[m1] = 0.4
# Pm[,1] = Pm[,1]/int(Pm[,1])


f0 = which(frange==-1)
f1 = which(frange==1)
Pf = matrix(0,Nf,Tsteps+1)
breaks = runif(n=length(mrange_orig)-1)
breaks = sort(breaks)
f_init = c(breaks[1],diff(c(breaks,1)))
Pf[f0:f1,1] = f_init
Pf[f0,1] = .4
Pf[f1,1] = .6
Pf[,1] = Pf[,1]/int(Pf[,1])
# p = .4
# Pf[,1] = p*dnorm(frange,-1,mix_sigma2)+(1-p)*dnorm(frange,1,mix_sigma2)

t = 1
for(t in 1:Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf)

	for(j in 1:Nf){
		y = frange[j]
		weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults)
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pm_beforemut = matrix(0,Nm)
	for(i in 1:Nm){
		Pm_beforemut[i] = int(pxy[i,])
	}
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + mut_prob/2*c(0,Pm_beforemut[1:Nm-1])
	Pm[,t+1] = Pm_aftermut
	Pf[,t+1] = Pf_adults
}


# #########
# Pm = matrix(0,2,Tsteps+1)
# Pf = matrix(0,2,Tsteps+1)
# Pm[,1] = c(.6,.4)
# Pf[,1] = c(.4,.6)
# alpha = 0
# for(t in 1:Tsteps){
	# Pm_adults = Pm[,t]
	# Pf_adults = Pf[,t]
	# pxy = matrix(0,2,2)
	# weight = c(1,1)+c(alpha,0)
	# z = sum(weight*Pm_adults)
	# pxy[,1] = Pf_adults[1]*weight*Pm_adults/z
	# weight = c(1,1)+c(0,alpha)
	# z = sum(weight*Pm_adults)
	# pxy[,2] = Pf_adults[2]*weight*Pm_adults/z
	# Pm_beforemut = matrix(0,2)
	# for(i in 1:2){
		# Pm_beforemut[i] = sum(pxy[i,])
	# }
	# Pm_aftermut = matrix(0,2)
	# Pm_aftermut = matrix(c((1-mut_prob),mut_prob,mut_prob,1-mut_prob),2,2)%*%Pm_beforemut
	# Pm[,t+1] = Pm_aftermut
	# Pf[,t+1] = Pf_adults
	# }
