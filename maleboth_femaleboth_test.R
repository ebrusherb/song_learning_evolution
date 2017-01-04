step = 0.01 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# f_sigma = #variance of female distribution(s)
# m_sigma = 0.1 #variance of male distribution(s)
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
midpt =ceiling(Nf/2)

nonzero_thresh = 1e-5
perc_thresh = 1e-10


sigma = 0.5
f_sigma = 1
m_sigma =0.1
mut_prob =0.1
pf=1
pm=1

f_init = dnorm(frange,fmin,f_sigma)
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,m_sigma)
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

# # k=100
# t0 = proc.time()

# total_breeding_mat = array(0,dim=c(Nm,Nf,Nf))
	# for(i in (m0-k):(m0+k)){
		# total_breeding_mat[i,,] = t(outer(pxy[i,],m_adults_prefs_normed[i,],'*'))
	# }
	
# print(proc.time()-t0)

# t0 = proc.time()

# total_breeding_mat = array(0,dim=c(Nm,Nf,Nf))

# for(i in (m0-k):(m0+k)){
	# for(j in 1:Nf){
		# total_breeding_mat[i,j,] = pxy[i,j] * 
	# }
# }

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

Pm = array(0,dim=c(Nm,Nf,Tsteps+1)) #probability of male songs and preferences over time
Pm[,,1] = outer(m_init,f_init,'*')

t = 1
perc = 0


Pm_adults = Pm[,,t]
	m_adult_songs = apply(Pm_adults,1,sum)
	divisor = m_adult_songs
	divisor[m_adult_songs==0] = 1
	m_adults_prefs_normed = Pm_adults / matrix(rep(divisor,times = Nf),nrow=Nm)
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	total_breeding_mat = array(0,dim=c(Nm,Nf,Nf))
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)	
	fixed_weight[fixed_weight==0] = minweight
	fixed_weight = fixed_weight/sum(fixed_weight)
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*m_adult_songs)
		if(z>1e-20){
			pxy[,j] = Pf_adults[j]*weight*m_adult_songs/z
			}
	}