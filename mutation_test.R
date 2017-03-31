setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
source('saveit.R')

trait_chunk_num = 301
sigmay2 = 2
sigmax2 = 0.8
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob = 0.0

steps = 20000
store = 1
source('range_setup.R')

f_init = dnorm(frange,fmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)

sigma2_vals = seq(0.1,1.9,length.out=10)

xs = length(sigma2_vals)

k1_vals = c(7,7,15,21,35)
k2_vals = c(35,25,21,21,35)

i = 6

sigma2 = sigma2_vals[i]
	
continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

j = 3

k1 = k1_vals[j]
k2 = k2_vals[j]

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

v = c(1,sigmax2,0)
p = solve(m,v)
p[1] = 0
		
m_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
m_init = c(m_init,rev(m_init[1:(length(m_init)-1)]))

mut_prob = 0
p1_step_song_dist = dynamics_bothmut()	

mut_prob = 0.01
p2_step_song_dist = dynamics_bothmut()

keep = seq(1,steps+1,by=200)
p1_step_song_dist$Pm=p1_step_song_dist$Pm[, keep]
p2_step_song_dist$Pm=p2_step_song_dist$Pm[, keep]
p1_step_song_dist$Pf=p1_step_song_dist$Pf[, keep]
p2_step_song_dist$Pf=p2_step_song_dist$Pf[, keep]
p1_step_song_dist$pxy=p1_step_song_dist$pxy[,, keep]
p2_step_song_dist$pxy=p2_step_song_dist$pxy[,, keep]
p1_step_song_dist$z=p1_step_song_dist$z[, keep]
p2_step_song_dist$z=p2_step_song_dist$z[, keep]

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

trait_chunk_num = 3

sigma2_vals = seq(0.1,1.9,length.out=10)

xs = length(sigma2_vals)

k1_vals = c(35,25,45,21,15)
k2_vals = c(35,35,45,21,21)

i = 6

sigma2 = sigma2_vals[i]
	
continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

j = 3
 
k1 = k1_vals[j]
k2 = k2_vals[j]

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

mut_prob = 0
p1_step_pref_dist = dynamics_bothmut()	

mut_prob = 0.01
p2_step_pref_dist = dynamics_bothmut()

keep = seq(1,steps+1,by=200)
p1_step_pref_dist$Pm=p1_step_pref_dist$Pm[, keep]
p2_step_pref_dist$Pm=p2_step_pref_dist$Pm[, keep]
p1_step_pref_dist$Pf=p1_step_pref_dist$Pf[, keep]
p2_step_pref_dist$Pf=p2_step_pref_dist$Pf[, keep]
p1_step_pref_dist$pxy=p1_step_pref_dist$pxy[,, keep]
p2_step_pref_dist$pxy=p2_step_pref_dist$pxy[,, keep]
p1_step_pref_dist$z=p1_step_pref_dist$z[, keep]
p2_step_pref_dist$z=p2_step_pref_dist$z[, keep]

saveit(p1_step_song_dist=p1_step_song_dist,p2_step_song_dist=p2_step_song_dist,p1_step_pref_dist=p1_step_pref_dist,p2_step_pref_dist=p2_step_pref_dist,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/mutation_test.Rdata')