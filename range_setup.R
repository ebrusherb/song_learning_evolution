#### NEED TO CHANGE SETUP SO -1 is always included in mrange/frange
## ---- parameters ----------------
# trait_step = 0.1 #step size of trait space
# int_step = trait_step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# f_sigma = #variance of female distribution(s)
# m_sigma = 0.1 #variance of male distribution(s)
# steps = #how many generations

# mrange = seq(-8,8,by=trait_step) #range of male songs
mrange = rev(seq(-1,-15,length.out=ceiling(trait_chunk_num/2)))
trait_step = diff(mrange)[1]
mrange = c(mrange,seq(-1+trait_step,13,length.out=floor(trait_chunk_num/2)))
Nm = length(mrange) 
if(Nm%%2==0){
	mrange = c(mrange[1]-trait_step,mrange)
	Nm = Nm+1
}
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=trait_step) 

midpt = ceiling(Nm/2)

nonzero_thresh = 1e-5
perc_thresh = 1e-10
z_thresh = 0
