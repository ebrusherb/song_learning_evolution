## ---- parameters ----------------
step = 0.01 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma2 = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# fmix_sigma2 = #variance of female distribution(s)
# mmix_sigma2 = 0.1 #variance of male distribution(s)
# Tsteps = #how many generations

mrange = seq(-10,10,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-10,10,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)


## ---- integral -------------------------
int<-function(v){ #manual integration
	l=length(v)
	int_step/2*(sum(2*v)-v[1]-v[l])
}

x = runif(
 10)


## ---- dynamics -----------------------

dynamics <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

# t = 1
for(t in 1:Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma2)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults) #normalization factor
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pm_beforemut = matrix(0,Nm)
	for(i in 1:Nm){
		Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	}
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) 
		+ mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
	Pm[,t+1] = Pm_aftermut
	Pf[,t+1] = Pf_adults
}
pop_dens = list(Pm=Pm,Pf=Pf)
return(pop_dens)
}

## ---- peaks -----------------
sigma2 = 0.1
mut_prob = 0.01
fmix_sigma2 = 1
Tsteps = 60

m_init = array(0, dim = c(Nm,1))
m_init[m0] = 0.6
m_init[m1] = 0.4
m_init[,1] = m_init[,1]/int(m_init[,1])

p = .4
f_init = p*dnorm(frange,-1,fmix_sigma2)+(1-p)*dnorm(frange,1,fmix_sigma2)

pop_dens = dynamics()
Pm = pop_dens$Pm
Pf = pop_dens$Pf

## ---- explanation -----------------
v = Pm[,1]
Tsteps = 1 
m_init = v
m_init = m_init/int(m_init)

p = .4
f_init = p*dnorm(frange,-1,fmix_sigma2)+(1-p)*dnorm(frange,1,fmix_sigma2)
# f_init = .3*dnorm(frange,-1,fmix_sigma2)+.3*dnorm(frange,0,fmix_sigma2)+.4*dnorm(frange,1,fmix_sigma2)

Pm2 = matrix(0,Nm,Tsteps+1)
Pm2[,1] = m_init

Pf2 = matrix(0,Nf,Tsteps+1)
Pf2[,1] = f_init

pop_dens2 = dynamics()
Pm2 = pop_dens2$Pm
Pf2 = pop_dens2$Pf

# break down the mating and preference probabilities
c = 1e-10 #recognition cutoff
d = sqrt(-2*sigma2*log(sqrt(2*pi*sigma2)*c)) 
#^distance / difference at which a female can recognize a male
recognized = array(0,dim=c(Nm,1))

for(i in 1:Nm){
	x = mrange[i]
	y1 = x - d
	y2 = x + d
	w1 = which(frange>=x-d)
	w2 = which(frange<=x+d)
	recognized[i] = int(Pf2[intersect(w1,w2),1]
	*dnorm(x-frange[intersect(w1,w2)],mean=0,sd=sqrt(sigma2))) 
	#^how many females recognize each male, weighted by their preferences
}

preferences = array(0,dim=c(Nm,Nf))
female_tots = matrix(0,Nf)

for(j in 1:Nf){
	y = frange[j]
	# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
	weight = dnorm(mrange,mean=y,sd=sqrt(sigma2))
	# weight = matrix (0,Nf,1)
	# weight[c(f0,x1)] = 1
	# weight[j] = 1+alpha
	z = int(weight*Pm2[,Tsteps])
	# if(z!=0){
		# pxy[,j] = Pf2[j,1]*weight*Pm2[,Tsteps]/z
		# }
	preferences[,j] = weight #preference given by each female to each male
	female_tots[j] = z #total preferences given by each female
}

## ---- playing -----------------
sigma2 = 1
mut_prob = 0.01
fmix_sigma2 = 1
mmix_sigma2 = 0.05
Tsteps = 200

#option: initial male distribution uniformly randomly distributed
m_init = runif(n=length(mrange_orig))
m_init = c(array(0,dim=c(m0-1,1)),m_init,array(0,dim=c(Nm-m1,1)))
m_init = m_init/int(m_init)
# # option: initial male distribution is two delta functions
# m_init = array(0, dim = c(Nm,1))
# m_init[m0] = 0.6
# m_init[m1] = 0.4
# m_init[,1] = m_init[,1]/int(m_init[,1])
# # option: initial male distribution is a mix of two normal distributions
# p = 0.6
# m_init = p*dnorm(mrange,-1,mmix_sigma2)+(1-p)*dnorm(mrange,1,mmix_sigma2)

# #option: initial female distribution uniformly randomly distributed
# f_init = runif( n=length(frange_orig)-1)
# f_init = c(array(0,dim=c(f0-1,1)),f_init,array(0,dim=c(Nf-f1,1)))
# f_init = f_init/int(f_init)
# #option: initial female distribution is two delta functions
# f_init = array(0, dim=c(Nf,1))
# f_init[f0,1] = .4
# f_init[f1,1] = .6
# f_init = f_init/int(f_init[,1])
#option: initial female distribution is mix of two normal distributions
p = .4
f_init = p*dnorm(frange,-1,fmix_sigma2)+(1-p)*dnorm(frange,1,fmix_sigma2)
# f_init = .3*dnorm(frange,-1,fmix_sigma2)+.3*dnorm(frange,0,fmix_sigma2)+.4*dnorm(frange,1,fmix_sigma2)

Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

# t = 1
for(t in 1:Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf)
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma2))
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
	Pm_aftermut = (1-mut_prob)*Pm_beforemut 
	+ mut_prob/2*c(Pm_beforemut[2:Nm],0) + mut_prob/2*c(0,Pm_beforemut[1:Nm-1])
	Pm[,t+1] = Pm_aftermut
	Pf[,t+1] = Pf_adults
}
		

