setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution/')
source('dynamics_pxy.R')
library(RColorBrewer)
library(pracma)
library(PearsonDS)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 301
sigmay2 = 2
sigmax2 = 0.8
sigma2 = 0.9
mut_prob = 0.0
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)

steps = 5000
store = 1
source('range_setup.R')

k1 = 45
k2 = 45

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

f_init = dnorm(frange,fmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init1 = f_init/sum(f_init)

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

v = c(1,sigmay2,0)
p = solve(m,v)
p[1] = 0

f_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
f_init = c(f_init,rev(f_init[1:(length(f_init)-1)]))
f_init2 = f_init

excess_kurt = -0.1
f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
f_init3 = f_init/sum(f_init)

excess_kurt = 0.01
f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
f_init4 = f_init/sum(f_init)

z = dnorm(mrange,mean=-1,sd=sqrt(sigmax2+sigma2))