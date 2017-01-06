## ---- parameters ----------------
trait_step = 0.1 #step size of trait space
int_step = trait_step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# f_sigma = #variance of female distribution(s)
# m_sigma = 0.1 #variance of male distribution(s)
# steps = #how many generations

mrange = seq(-8,8,by=trait_step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=trait_step) 
frange = mrange #range of female preferences
# frange = seq(-7.5,7.5,by=trait_step) 
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=trait_step)
midpt = ceiling(Nf/2)

nonzero_thresh = 1e-5
perc_thresh = 1e-10
z_thresh = 0
