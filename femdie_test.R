dynamics_full <-function(){
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
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma)*exp(-(mrange-y)^2/(2*sigma))
		weight = dnorm(mrange,mean=y,sd=(sigma)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults) #normalization factor
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pm_beforemut = apply(pxy,1,int)
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
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
pop_dens = list(Pm=Pm[,(1):Tsteps],Pf=Pf[,(1):Tsteps])
return(pop_dens)
}

dynamics_femdie <-function(){
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
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma)*exp(-(mrange-y)^2/(2*sigma))
		weight = dnorm(mrange,mean=y,sd=(sigma)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults) #normalization factor
		if(z>1e-20){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pf_adults[apply(pxy,2,int)==0]=0
	Pm_beforemut = apply(pxy,1,int)/int(Pf_adults)
	Pf_adults = Pf_adults/int(Pf_adults)
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
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
pop_dens = list(Pm=Pm[,(1):Tsteps],Pf=Pf[,(1):Tsteps])
return(pop_dens)
}

dynamics_intfix <-function(){
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
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma)*exp(-(mrange-y)^2/(2*sigma))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma)) #female preference function	# weight = matrix (0,Nf,1)
		weight = weight/sum(weight)
		minweight = 10^max(floor(log(min(weight[which(weight>0)]),base=10)),-320)
		weight[weight==0] = minweight
		# weight = 0.5*(diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mean=y,sd=sigma))-diff(1-pnorm(c(mrange[1]-step/2,mrange+step/2),mean=y,sd=sigma)))
		# weight = -diff(1-pnorm(c(mrange[1]-step/2,mrange+step/2),mean=y,sd=sigma))
		# weight = weight + 10^(floor(log(min(weight[which(weight>0)]),base=10))-1)
		# # n = sum(weight)
		# weight = 0.5*(weight + c(rev(weight[1:(2*j-1)]),weight[(2*j):Nf])) 
		# weight = weight + 10^(floor(log(min(weight[which(weight>0)]),base=10))-1)
		# weight = weight / sum(weight) * n
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = sum(weight*Pm_adults) #normalization factor
		if(z>1e-20){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
		# j = ceiling(Nf/2)
		# y = frange[j]
		# weight = diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mean=y,sd=sigma))
		# # n = sum(weight)
		# weight = 0.5*(weight + rev(weight))
		# weight = weight + 10^(floor(log(min(weight[which(weight>0)]),base=10))-1)
		# # weight = weight / sum(weight) * n
		# # weight[c(f0,x1)] = 1
		# # weight[j] = 1+alpha
		# z = sum(weight*Pm_adults) #normalization factor
		# if(z>1e-20){
			# pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			# }
	# for(j in (ceiling(Nf/2)+1):Nf){
		# y = frange[j]
		# # weight = 1/sqrt(2*pi*sigma)*exp(-(mrange-y)^2/(2*sigma))
		# # weight = dnorm(mrange,mean=y,sd=sqrt(sigma)) #female preference function	# weight = matrix (0,Nf,1)
		# weight = diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mean=y,sd=sigma))
		# # n = sum(weight)
		# weight = 0.5*(weight + c(weight[1:(2*j-Nf-1)],rev(weight[(2*j-Nf):Nf])))
		# weight = weight + 10^(floor(log(min(weight[which(weight>0)]),base=10))-1)
		# # weight = weight / sum(weight) * n
		# # weight[c(f0,x1)] = 1
		# # weight[j] = 1+alpha
		# z = sum(weight*Pm_adults) #normalization factor
		# if(z>1e-20){
			# pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			# }
	# }

	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pf_adults[apply(pxy,2,sum)==0]=0
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
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
pop_dens = list(Pm=Pm[,(1):Tsteps],Pf=Pf[,(1):Tsteps])
return(pop_dens)
}

Tsteps = 10

sigma = 0.001
f_sigma = 0.5
m_sigma = 0.01
mut_prob = 0
f_init = dnorm(frange,fmin,f_sigma)
m_init = dnorm(mrange,mmin,m_sigma)

pop_dens = dynamics_full()
pop_dens2 = dynamics_femdie()

f_init = diff(pnorm(c(frange[1]-step/2,frange+step/2),fmin,f_sigma))
f_init = 0.5*(f_init+c(rev(f_init[1:(2*m0-1)]),f_init[(2*m0):Nm]))
m_init = diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mmin,m_sigma))
m_init = 0.5*(m_init+c(rev(m_init[1:(2*m0-1)]),m_init[(2*m0):Nm]))

pop_dens3 = dynamics_intfix()

f_init = dnorm(frange,fmin,f_sigma)
m_init = dnorm(mrange,mmin,m_sigma)
m_init = m_init + 1e-300
m_init = m_init / int(m_init)

pop_dens4 = dynamics_femdie()

Tsteps = 100

sigma = 0.3162
f_sigma = 1
m_sigma = 0.01
mut_prob = 0.01

f_init = dnorm(frange,fmin,f_sigma)
m_init = dnorm(mrange,mmin,m_sigma)

pop_dens = dynamics_full()
pop_dens2 = dynamics_femdie()

f_init = diff(pnorm(c(frange[1]-step/2,frange+step/2),fmin,f_sigma))
f_init = 0.5*(f_init+c(rev(f_init[1:(2*m0-1)]),f_init[(2*m0):Nm]))
m_init = diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mmin,m_sigma))
m_init = 0.5*(m_init+c(rev(m_init[1:(2*m0-1)]),m_init[(2*m0):Nm]))

pop_dens3 = dynamics_intfix()

f_init = dnorm(frange,fmin,f_sigma)
m_init = dnorm(mrange,mmin,m_sigma)
m_init = m_init + 1e-300
m_init = m_init / int(m_init)

pop_dens4 = dynamics_femdie()


f_init = diff(pnorm(c(frange[1]-step/2,frange+step/2),fmin,f_sigma))
f_init = 0.5*(f_init+c(rev(f_init[1:(2*m0-1)]),f_init[(2*m0):Nm]))
f_init = f_init + 10^(floor(log(min(f_init[which(f_init>0)]),base=10))+1)
f_init = f_init / sum(f_init)

m_init = diff(pnorm(c(mrange[1]-step/2,mrange+step/2),mmin,m_sigma))
m_init = 0.5*(m_init+c(rev(m_init[1:(2*m0-1)]),m_init[(2*m0):Nm]))
m_init = m_init + 1e-300
m_init = m_init / sum(m_init)

pop_dens5 = dynamics_intfix()


