problem = c()
mint = c()
fint = c()
for(i in 1:prod(dim(Pm_onepop))){
	if(abs(sum(Pm_onepop[[i]][,Tend])-1)>1e-10){
		problem=c(problem,i)
	}
	mint = c(mint,sum(Pm_onepop[[i]][,Tend]))
	fint = c(fint,sum(Pf_onepop[[i]][,Tend]))
}

####WHY AREN"T THE INTEGRALS EQUAL TO ONE?!?!?!?!? OH CRAP OH CRAP 

# plot(600:700,Pm_onepop[[1]][600:700,Tend],t='n',ylim=c(0,10))

# for(w in intersect(problem,which(pmat>0))){
	# lines(Pm_onepop[[w]][,Tend],t='o')
# }

# lines(dnorm(mrange,-1,0.001),col='red')

# i = which.min(mint)

# Tsteps =10

# sigma = smat[i]
# f_sigma = fmat[i]
# m_sigma = mmat[i]
# mut_prob = pmat[i]
# f_init = dnorm(frange,fmin,f_sigma)
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
# f_init = f_init/sum(f_init)
# m_init = dnorm(mrange,mmin,m_sigma)
# m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
# m_init = m_init / sum(m_init)

# pop_dens = dynamics()