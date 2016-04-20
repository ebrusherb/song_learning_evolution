## ---- parameters ----------------
step = 0.01 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma2 = #variance of female preference function
mut_prob = 0.01 #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# fmix_sigma2 = #variance of female distribution(s)
# mmix_sigma2 = 0.1 #variance of male distribution(s)
# Tsteps = #how many generations

mrange = seq(-7.5,7.5,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-7.5,7.5,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)

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
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
	Pm[,t+1] = Pm_aftermut
	Pf[,t+1] = Pf_adults
}
pop_dens = list(Pm=Pm,Pf=Pf)
return(pop_dens)
}

## ---- variance_sweep ---------------
Tsteps2 = 1000
Tsteps = Tsteps2
pm = 1
pf = 1


sigma2 = 0.01
mmix_sigma2 = 0.1

fmix_sigma2_test = c(0.001,0.005,0.0075,0.01,0.05,0.1,0.25)
Nftest = length(fmix_sigma2_test)
Pm_test = as.list(1:Nftest)
Pf_test = as.list(1:Nftest)

for(f in 1:Nftest){
	fmix_sigma2 = fmix_sigma2_test[f]
	f_init = pf*dnorm(frange,fmin,fmix_sigma2)+(1-pf)*dnorm(frange,fmax,fmix_sigma2)
	m_init = pm*dnorm(mrange,mmin,mmix_sigma2)+(1-pm)*dnorm(mrange,mmax,mmix_sigma2)
	pop_dens = dynamics()
	Pm_test[[f]] = pop_dens$Pm	
	Pf_test[[f]] = pop_dens$Pf
}

Tsteps = d[2]

sixpal = brewer.pal(6,'Set1')

f=1;
plot(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_test[[f]][,Tsteps2],type='l',col=sixpal[f],ylim=c(0,0.5),xlim=c(-1.5,1.5))
for(f in 2:Nftest){
	lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_test[[f]][,Tsteps2],type='l',col=sixpal[f])
}
legend(1,.4,legend=fmix_sigma2_test,col=sixpal,lty=1,bty='n')

par(mfrow=c(2,3))
for(f in 1:Nftest){
	plot(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_test[[f]][,Tsteps2],type='l',col=sixpal[1],ylim=c(0,0.5),xlim=c(-1.5,1.5),main=fmix_sigma2_test[f])
	lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pf_test[[f]][,Tsteps2],col='black')
}

par(mfrow=c(1,1))
f=1
plot(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_test[[f]][,Tsteps2],type='l',col=sixpal[f],ylim=c(0,0.5),xlim=c(-1.5,1.5),main=fmix_sigma2_test[f])
lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pf_test[[f]][,Tsteps2],col=sixpal[f],type='o',lty=2)
for(f in 2:Nftest){
	lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pf_test[[f]][,Tsteps2],col=sixpal[f])
	lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pf_test[[f]][,Tsteps2],col=sixpal[f],type='o',lty=2)
}
legend(1,.4,legend=fmix_sigma2_test,col=sixpal,lty=1,bty='n')

four_var_test_m = array(NA,dim=c(Nftest,1))
four_var_test_f = array(NA,dim=c(Nftest,1))
for(f in 1:Nftest){
	eq_pop = Pm_test[[f]][,Tsteps2]
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourierCoefficients[1] / length(subset)
	l = lm(log(coeffs[2:13]) ~ poly(1:12,2,raw=TRUE));
	four_var_test_m[f] = -summary(l)$coefficients[3,1]
	fmix_sigma2 = fmix_sigma2_test[f]
	f_init = pf*dnorm(frange,fmin,fmix_sigma2)+(1-pf)*dnorm(frange,fmax,fmix_sigma2)
	fourComps = fft(f_init);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourierCoefficients[1] / length(subset)
	l = lm(log(coeffs[2:13]) ~ poly(1:12,2,raw=TRUE));
	four_var_test_f[f] = -summary(l)$coefficients[3,1]
}

save(Pm_test=Pm_test,Pf_test=Pf_test,sigma2=sigma2,fmix_sigma2_test=fmix_sigma2_test,mmix_sigma2,mut_prob,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/effect_of_f_wo_peaks.Rdata')