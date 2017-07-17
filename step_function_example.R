setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_by_mode.R')
source('saveit.R')
library(RColorBrewer)

trait_chunk_num = 75
source('range_setup.R')
mrange=c(rev(seq(0,-10,by=-0.25)),c(seq(0.25,8,by=0.25)))
Nm = length(mrange)
midpt = which(mrange==-1)
frange=mrange
Nf = Nm
sigmay2 = 2
sigmax2 = 0.1
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob = 0.01

sigma2 = 1.5

k1 = 9
k2 = 12

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

v = c(1,sigma2,0)
p = solve(m,v)
p[1] = 0


fixed_weight1 = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
fixed_weight1 = c(fixed_weight1,rev(fixed_weight1[1:(length(fixed_weight1)-1)]))

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight2 = continuous_weight/sum(continuous_weight)


m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

f_init = dnorm(frange,fmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)

fixed_weight = fixed_weight2


w = setdiff(which(round(mrange,3)>=-7),which(round(mrange,3)>5))

fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12
col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.7,0.5,0.02,0.15)
omarg = c(0.1,0.4,0.35,0.0)

width = 3.8
height = 3.25

pdf(file='/Users/eleanorbrush/Desktop/step_function_example.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.5,0))

plot(mrange[w]+1,fixed_weight2[w],t='o',xlim=range(mrange[w]+1)+c(-.01,.01),yaxt='n',xlab='',ylab='',lwd=lwd,col=col_vec[1],ylim=range(c(fixed_weight1,fixed_weight2))+c(-.00,0),xaxt='n')
points(mrange[w]+1,fixed_weight1[w],t='o',col=col_vec[2],lwd=lwd)
axis(2,at=c(p[2],p[3],c(0,0.04,0.08)),labels=c(expression(p[1]),expression(p[2]),c(0,0.04,0.08)))
axis(1,at=c(-4,-1,0,1,4),labels=c(-4,-1,0,1,4))
mtext('Difference in preference and song, y-x',side=1,line=2,at=-0.3,cex=largefontsize/smallfontsize)
mtext('Preference',side=2,line=1.7,at=p[3]/2,cex=largefontsize/smallfontsize)

dev.off()