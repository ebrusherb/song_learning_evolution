s=2
f=5
m=4

sigma2 = sigma2_vals[s]
fmix_sigma2 = fmix_sigma2_vals[f]
mmix_sigma2 = mmix_sigma2_vals[m]

mend = Pm_onepop[[s,f,m]][,Tsteps]
fend= Pf_onepop[[s,f,m]][,Tsteps]

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

