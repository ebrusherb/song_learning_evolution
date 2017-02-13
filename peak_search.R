sigmay2_init = 4
sigmax2_init = 2
sigma2 = 0.01
mut_prob = 0
pf = 1
pm = 1
steps = 1000
store = 1

trait_chunk_num = 5
source("/Users/eleanorbrush/Documents/research/song_learning_evolution/range_setup.R")

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
continuous_weight[continuous_weight==0] = minweight

fixed_weight = continuous_weight/sum(continuous_weight)

pref_chunk_num = 15
fixed_weight = discretize(continuous_weight)
fixed_weight = fixed_weight / sum(fixed_weight)

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
j = midpt
# j = min(which(mrange>0))
s = sign(midpt-j+0.5)
toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
hold = f_init
hold[toreplace[1]:toreplace[2]] = f_init[pull[1]:pull[2]]
if(toreplace[1]>1){
	hold[1:(toreplace[1]-1)] = f_init[(pull[2]+1):Nf]}
f_init = hold

m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

##### good:
alpha = 1
fixed_weight = rep(1,trait_chunk_num)
fixed_weight[midpt] = 1+alpha
fixed_weight = fixed_weight / sum(fixed_weight)
f_init = c(0.1,0.2,0.4,0.2,0.1)
m_init = c(0.4,0.2,0.2,0.1,0.1)
####### 

pop_dens = dynamics()
layout(matrix(1:2,ncol=2))
plot(mrange,fixed_weight,t='o')
plot(mrange,pop_dens$Pm[,steps],t='o')


##### explicitly compare preference functions
continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
minweight = 10^(-320)
continuous_weight[continuous_weight==0] = minweight
pref_chunk_num = trait_chunk_num
fixed_weight = discretize(continuous_weight)
fixed_weight = fixed_weight / sum(fixed_weight)
f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

ex = sum(mrange*fixed_weight)
vx = sum((mrange-(ex))^2*fixed_weight)
m4 = sum((mrange-(ex))^4*fixed_weight)
kurt1 = m4 / vx^2

pop_dens = dynamics()
print(round(pop_dens$Pm[,steps],5))

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
minweight = 10^(-1)
continuous_weight[continuous_weight<10^(-10)] = minweight
pref_chunk_num = trait_chunk_num
fixed_weight = discretize(continuous_weight)
fixed_weight = fixed_weight / sum(fixed_weight)
f_init = pop_dens$Pf[,steps]
m_init = pop_dens$Pm[,steps]

ex = sum(mrange*fixed_weight)
vx = sum((mrange-(ex))^2*fixed_weight)
m4 = sum((mrange-(ex))^4*fixed_weight)
kurt2 = m4 / vx^2

pop_dens2 = dynamics()
print(round(pop_dens2$Pm[,steps],5))