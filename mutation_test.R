setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
source('saveit.R')
source('recursion_all.R')
library(RColorBrewer)
library(pracma)
library(PearsonDS)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 281
sigmay2 = 2
sigmax2 = 0.8
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob_nonzero = 0.001
steps = 20000
store_vec = c(1,100,10000,20000)
source('range_setup.R')

sigma2 = 1.1
	
continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

f_init = dnorm(frange,fmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init_norm = f_init/sum(f_init)

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init_norm = m_init/sum(m_init)

k1 = 15
k2 = 25

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

v = c(1,sigmax2,0)
p = solve(m,v)
p[1] = 0
		
m_init_step = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
m_init_step = c(m_init_step,rev(m_init_step[1:(length(m_init_step)-1)]))

k1 = 45
k2 = 45

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

v = c(1,sigmay2,0)
p = solve(m,v)
p[1] = 0

f_init_step = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
f_init_step = c(f_init_step,rev(f_init_step[1:(length(f_init_step)-1)]))

f_init = f_init_norm
m_init = m_init_norm
mut_prob = mut_prob_nonzero
p2 = dynamics_memory(store=TRUE)

f_init = f_init_norm
m_init = m_init_step

mut_prob = 0
p1_step_song_dist = dynamics_memory(store=TRUE)	

mut_prob = mut_prob_nonzero
p2_step_song_dist = dynamics_memory(store=TRUE)

f_init = f_init_step
m_init = m_init_norm

mut_prob = 0
p1_step_pref_dist = dynamics_memory(store=TRUE)	

mut_prob = 0.01
p2_step_pref_dist = dynamics_memory(store=TRUE)


saveit(p2=p2,p1_step_song_dist=p1_step_song_dist,p2_step_song_dist=p2_step_song_dist,p1_step_pref_dist=p1_step_pref_dist,p2_step_pref_dist=p2_step_pref_dist,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,store_vec=store_vec,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/mutation_test.Rdata')

load('/Users/eleanorbrush/Documents/research/song_learning_evolution/mutation_test.Rdata')


r = recursion_all(sigmax2,sigmay2,sigma2,rho)

t_toplot = c(1:4)

p = c()

for(t in store_vec[t_toplot]){
	p = cbind(p,dnorm(mrange,mean=-1,sd=sqrt(r[3,1,t]))/sum(dnorm(mrange,mean=-1,sd=sqrt(r[3,1,t]))))
}

col_vec = brewer.pal(9,'Greys')[c(4,5,7,9)]
# col_vec = c(brewer.pal(9,'Set1')[c(3,1:2)],'black')
col_vec = brewer.pal(4,'Reds')
# col_vec = brewer.pal(9,'YlOrRd')[c(1,3,5,7)]
# col_vec = brewer.pal(4,'Spectral')
col_vec = brewer.pal(11,'RdBu')[c(11,9,3,2)]
lwd = 2
marg = c(0.38,0.3,0.0,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 5
ylim = c(0,0.1)
t = dim(p2$Pm_store)[2]
w = which(p2$Pm_store[,t]>1e-10)

pdf('/Users/eleanorbrush/Desktop/mutation_sensitivity.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=FALSE))

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p[w,t],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)
legend(-4.5,0.1,legend=store_vec[t_toplot],lty=1,col=col_vec,lwd=lwd,bty='n')

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p2$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p1_step_song_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p2_step_song_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p1_step_pref_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
for(t in 1:length(t_toplot)){
	lines(mrange[w]+1,p2_step_pref_dist$Pm_store[w,t_toplot[t]],lwd=lwd,col=col_vec[t])
}
mtext('Song, x',side=1,line=1.9,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)

dev.off()