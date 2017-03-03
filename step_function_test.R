steps = 100
sigma2 = 0.5
sigmay2_init = 1.1
sigmax2_init = 0.5
trait_chunk_num = 301
pref_chunk_num = 301
source("range_setup.R")

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

var_step <-function(k1,k2,delta){
	fixed_weight1 = rep(0,Nm)
	fixed_weight1[midpt+(-k2:k2)] = 1
	fixed_weight1[midpt+(-k1:k1)] = 1+delta
	fixed_weight1 = fixed_weight1/sum(fixed_weight1)
	v = sum((mrange+1)^2*fixed_weight1)-sigma2
	return(v)
}

k1_vals = seq(1,19,by=2)
k2_vals = k1_vals + 17
x1 = length(k1_vals)
x2 = length(k2_vals)
equilibrium_step = as.list(rep(NA,(x1*x2*2)))
dim(equilibrium_step) = c(x1,x2,2)

for(i in 1:x1){
	for(j in 1:x2){
		k1 = k1_vals[i]
		k2 = k2_vals[j]
		if(sign(var_step(k1,k2,0))!=sign(var_step(k1,k2,10000))){
			delta = uniroot(var_step,k1=k1,k2=k2,lower=0,upper=10000)$root
			fixed_weight1 = rep(0,Nm)
			fixed_weight1[midpt+(-k2:k2)] = 1
			fixed_weight1[midpt+(-k1:k1)] = 1+delta
			fixed_weight1 = fixed_weight1/sum(fixed_weight1)
			
			continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
			fixed_weight2 = continuous_weight/sum(continuous_weight)
			
			fixed_weight = fixed_weight1
			p1 = dynamics()
			
			fixed_weight = fixed_weight2
			p2 = dynamics()	
			
			equilibrium_step[[i,j,1]] = p1$Pm[,steps]
			equilibrium_step[[i,j,2]] = p2$Pm[,steps]	
			}
	}
}
