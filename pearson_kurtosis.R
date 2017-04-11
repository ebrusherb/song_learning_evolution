setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_pxy.R')
source('saveit.R')
library(RColorBrewer)
library(pracma)
library(PearsonDS)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

trait_chunk_num = 173
sigmay2 = 2
sigmax2 = 0.8
sigma2 = 1.1
mut_prob = 0.0
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
minprob = 0

steps = 30000
store = steps+1
source('range_setup_long.R')

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

f_init = dnorm(frange,fmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init1 = f_init/sum(f_init)

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
fixed_weight = continuous_weight/sum(continuous_weight)

f_init = f_init1
p1 = dynamics_memory()

excess_kurt_vals = c(-0.5,-0.4,-0.3,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2)
mut_prob_vals = c(0,0.001,0.01,0.05)
xk = length(excess_kurt_vals)
xm = length(mut_prob_vals)

equilibrium = as.list(1:(xk*xm))
dim(equilibrium)=c(xk,xm)

for(k in 1:xk){
	for(m in 1:xm){
	excess_kurt = excess_kurt_vals[k]
	mut_prob = mut_prob_vals[m]
	
	f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
	f_init2 = f_init/sum(f_init)
	
	f_init = f_init2
	p2 = dynamics_memory()
	
	equilibrium[[k,m]] = p2$Pm
	}
}


# saveit(p1=p1,equilibrium=equilibrium,excess_kurt_vals=excess_kurt_vals,mut_prob_vals=mut_prob_vals,sigmax2=sigmax2,sigma2=sigma2,sigmay2=sigmay2,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/pearson_kurtosis.Rdata')

#### 
# col_vec = brewer.pal(9,'Set1')[-c(6,7)]
# lwd = 2
# marg = c(0.45,0.43,0.02,0.15)
# omarg = c(0.03,1,0.35,0.0)

# width = 6.5
# height = 4.5

# # pdf('/Users/eleanorbrush/Desktop/peak_example.pdf',width=width,height=height,family=fontfamily)

# par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
# layout(matrix(1:4,ncol=2,byrow=TRUE))

# t=20
# w1 = which(p1$Pm[,t]>1e-9)
# w2 = which(p2$Pm[,t]>1e-15)

# plot(mrange[w1]+1,f_init1[w1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(f_init1,f_init2)))
# lines(mrange[w1]+1,f_init2[w1],t='l',lwd=lwd,col='black',xlab='',ylab='')
# mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=1.7,at=mean(range(c(f_init1,f_init2))),cex=largefontsize/smallfontsize)

# plot(mrange[w1]+1,p1$Pm[w1,steps+1],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(p1$Pm[w1,t+1],p2$Pm[,steps+1])))
# points(mrange+1,p2$Pm[,steps+1],t='l',lwd=lwd,col='black',xlab='',ylab='')
# mtext('Song',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=1.7,at=mean(range(c(p1$Pm[w1,steps],p2$Pm[,steps]))),cex=largefontsize/smallfontsize)

# plot(mrange[w1]+1,p1$z[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range(c(p1$z[w1,t],p2$z[w1,t])))
# points(mrange[w1]+1,p2$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='')
# mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# mtext('Z',side=2,line=1.7,at=mean(range(c(p1$z[w1,t],p2$z[w2,t]))),cex=largefontsize/smallfontsize)

# # plot(mrange[w1]+1,(apply(p1$pxy[w1,,t],1,sum)-p1$Pm[w1,t])/p1$Pm[w1,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=range((apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],na.rm=TRUE),yaxt='n')
# # axis(2,at=c(-.1,0,0.2,.4),labels=c(-.1,0,0.2,.4))
# # lines(mrange[w2]+1,(apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],t='l',lwd=lwd,col='black,xlab='',ylab='')
# # mtext('Song',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# # mtext('Percent increase',side=2,line=1.7,at=mean(range((apply(p2$pxy[w2,,t],1,sum)-p2$Pm[w2,t])/p2$Pm[w2,t],na.rm=TRUE)),cex=largefontsize/smallfontsize)

# plot(mrange[w1]+1,p2$z[w1,t]-p1$z[w1,t],t='l',lwd=lwd,col='black',xlab='',ylab='',ylim=range(p2$z[,t]-p1$z[,t],na.rm=TRUE))
# # axis(2,at=0.0005*2*seq(-6,6,by=1))#,labels=c(expression(-1 %*% 10^{-3}),expression(-5 %*% 10^{-4}),0,expression(5 %*% 10^{-4}),expression(1 %*% 10^{-3})))
# mtext('Preference',side=1,line=1.5,at=0,cex=largefontsize/smallfontsize)
# mtext('Difference in Z',side=2,line=1.7,at=mean(range(p2$z[,t]-p1$z[,t],na.rm=TRUE)),cex=largefontsize/smallfontsize)

# # dev.off()

# print(c(sum((mrange+1)^4*p1$Pf[,1])/sigmay2^2,sum((mrange+1)^4*p2$Pf[,1])/sigmay2^2))

col_vec = brewer.pal(9,'Greys')[c(5,7,9)]
col_vec = c('black',brewer.pal(9,'Set1')[1])
lwd = 2
marg = c(0.38,0.3,0.0,0.15)
omarg = c(0.03,1,0.35,0.0)
ylim = c(0,0.22)

width = 6.5
height = 4.5

w = 87+(-27:27)

pdf('/Users/eleanorbrush/Desktop/pearson_kurtosis.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=TRUE))

for(i in 4:9){
	plot(mrange[w]+1,mrange[w],t='n',ylim=ylim,xlab='',ylab='')
	for(j in 1:2){
		lines(mrange[w]+1,equilibrium[[i,j]][w],lwd=lwd,col=col_vec[j])
	}
	mtext('Song',side=1,line=1.7,at=0,cex=largefontsize/smallfontsize)
	mtext('Frequency',side=2,line=1.7,at=mean(ylim),cex=largefontsize/smallfontsize)
	if(i ==4){		legend(-7.5,0.23,legend=mut_prob_vals[1:2],lty=1,col=col_vec[1:2],lwd=lwd,bty='n')}
}

dev.off()