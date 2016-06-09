d = dim(Pm_onepop)
P = prod(d)
perc = array(NA,dim=d)
t=Tend-1
thresh = 1e-10
w = c()
for(i in 1:P){
	perc[i] = max(abs(range(Pm_onepop[[i]][which(Pm_onepop[[i]][,t]>thresh),t+1]/Pm_onepop[[i]][which(Pm_onepop[[i]][,t]>thresh),t],na.rm=TRUE)-c(1,1)))
	sub = ind2sub(d,i)
	f = sub[2]
	m = sub[3]
	if(perc[i]>1e-4){w=c(w,i)}
}

# w = which(perc>=1e-4)
# w = intersect(w,union(which(as.vector(mmat)==0.001),which(as.vector(fmat)==0.001)))
print(length(w))


# plot(mrange,mrange,ylim=c(0,20),t='n')
# for( i in w){
	# lines(mrange,Pm_onepop[[i]][,200])
# }

# plot(mrange,mrange,ylim=c(0,20),t='n')
# for(i in setdiff(1:P,w)){
	# lines(mrange,Pm_onepop[[i]][,200])
# }

# Tsteps = 20000
# store = 1
# s = 1
# f = Nfs
# m = 1
# p = 2
# sigma = sigma_vals[s]
# f_sigma = f_sigma_vals[f]
# m_sigma = m_sigma_vals[m]
# f_init = dnorm(frange,fmin,f_sigma)
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
# f_init = f_init/sum(f_init)
# f_init = pf*f_init+(1-pf)*rev(f_init)
# m_init = dnorm(mrange,mmin,m_sigma)
# m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
# m_init = m_init / sum(m_init)
# m_init = pf*m_init+(1-pf)*rev(m_init)
# mut_prob = mut_prob_vals[p]
# pop_dens = dynamics()