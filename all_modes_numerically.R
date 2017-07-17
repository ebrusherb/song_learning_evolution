trait_chunk_num = 91
sigmay2 = 0.1
sigmax2 = 0.7
sigma2 = 1
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob = 0.0

steps = 50
store = 1
source('range_setup.R')

m_init = dnorm(mrange,-1,sqrt(sigmax2))
m_init = m_init/sum(m_init)

f_init = dnorm(frange,-1,sqrt(sigmay2))
f_init = f_init/sum(f_init)

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

d3 = dynamics_mode3()
v3 = apply(d3$Pm,2,function(v) sum((mrange+1)^2*v))

d7 = dynamics_mode7()
v7 = apply(d7$Pf,2,function(v) sum((mrange+1)^2*v))

d2 = dynamics_mode2()
v2 = apply(d2$Pm,2,function(v) sum((mrange+1)^2*v))

r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,0)