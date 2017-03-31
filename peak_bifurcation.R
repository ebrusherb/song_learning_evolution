setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
library(RColorBrewer)
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 101
sigmay2 = 2
sigmax2 = 0.8
sigma2 = 0.9
mut_prob = 0.01
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
minprob = 0

steps = 30
store = 1
source('range_setup.R')

k1 = 45
k2 = 45

k1 = 11
k2 = 6

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

v = c(1,sigmay2,minprob)
p = solve(m,v)
p[1] = minprob

f_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
f_init = c(f_init,rev(f_init[1:(length(f_init)-1)]))
f_init2 = f_init

f_init = f_init1
p1 = dynamics_bothmut()

f_init = f_init2
p2 = dynamics_bothmut()

# saveit(p1=p1,p2=p2,k1=k1,k2=k2,sigmay2=sigmay2,sigmax2=sigmax2,sigma2=sigma2,steps=steps,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/peak_example.Rdata')

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.45,0.43,0.02,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 4.5

# pdf('/Users/eleanorbrush/Desktop/peak_example.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:4,ncol=2,byrow=TRUE))

t=20
w1 = which(p1$Pm[,t]>1e-9)
w2 = which(p2$Pm[,t]>1e-15)

plot(mrange[w1]+1,f_init1[w1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(f_init1,f_init2)))
lines(mrange[w1]+1,f_init2[w1],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(range(c(f_init1,f_init2))),cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p1$Pm[w1,steps+1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(p1$Pm[w1,steps+1],p2$Pm[,steps+1])))
points(mrange+1,p2$Pm[,steps+1],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Song',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=1.7,at=mean(range(c(p1$Pm[w1,steps],p2$Pm[,steps]))),cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p1$z[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(p1$z[w1,t],p2$z[w1,t])))
points(mrange[w1]+1,p2$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='')
mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext('Z',side=2,line=1.7,at=mean(range(c(p1$z[w1,t],p2$z[w2,t]))),cex=largefontsize/smallfontsize)

# plot(mrange[w1]+1,(apply(p1$pxy[w1,,t],1,sum)-p1$Pm[w1,t])/p1$Pm[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range((apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],na.rm=TRUE),yaxt='n')
# axis(2,at=c(-.1,0,0.2,.4),labels=c(-.1,0,0.2,.4))
# lines(mrange[w2]+1,(apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],t='l',lwd=lwd,col='black,xlab='',ylab='')
# mtext('Song',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# mtext('Percent increase',side=2,line=1.7,at=mean(range((apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],na.rm=TRUE)),cex=largefontsize/smallfontsize)

plot(mrange[w1]+1,p2$z[w1,t]-p1$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='',ylim=range(p2$z[,t]-p1$z[,t],na.rm=TRUE))
# axis(2,at=0.0005*2*seq(-6,6,by=1))#,labels=c(expression(-1 %*% 10^{-3}),expression(-5 %*% 10^{-4}),0,expression(5 %*% 10^{-4}),expression(1 %*% 10^{-3})))
mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
mtext('Difference in Z',side=2,line=1.7,at=mean(range(p2$z[,t]-p1$z[,t],na.rm=TRUE)),cex=largefontsize/smallfontsize)

# dev.off()

print(c(sum((mrange+1)^4*p1$Pf[,1])/sigmay2^2,sum((mrange+1)^4*p2$Pf[,1])/sigmay2^2))