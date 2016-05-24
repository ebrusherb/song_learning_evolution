s=5
f=5
m=1

sigma2 = sigma2_vals[s]
fmix_sigma2 = fmix_sigma2_vals[f]
mmix_sigma2 = mmix_sigma2_vals[m]

i = 236

pf = 0.6
pm =0.4

mend = Pm_onepop[[i]][,Tend]
fend= Pf_onepop[[i]][,Tend]

morig = mend
forig = fend

mpeak_diff = m1-m0
fpeak_diff = f1-f0

mshift = matrix(0,Nm)
mshift[mpeak_diff:Nm] = mend[1:(Nm-mpeak_diff+1)]
fshift = matrix(0,Nf)
fshift[fpeak_diff:Nf] = fend[1:(Nf-fpeak_diff+1)]

m_init = pm*morig+(1-pm)*mshift
f_init = pf*forig + (1-pf)*fshift

Tsteps = 100
pop_dens = dynamics()

Pm=pop_dens$Pm
Pf=pop_dens$Pf

plot(Pm[,1],type='l',col='red',ylim=c(0,max(Pm)),xlab='Song',ylab='Density',main=paste("Pref width = ",sigma2_vals[s],"Male var = ",mmix_sigma2_vals[m]))
lines(Pf[,1],col='blue')
lines(Pm[,Tsteps+1])

