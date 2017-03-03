setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12 

########### equilibrium variance 
sigmax2_eq3 <- function(sigmax2,sigmay2,sigma2){
	sigmay2-sigma2
}

sigmax2_mode7 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[7,1,steps])
}

sigmax2_eq9 <- function(sigmax2,sigmay2,sigma2){
	(3*sigmay2-5*sigma2+sqrt(9*sigmay2^2-30*sigma2*sigmay2+sigma2^2))/6
}

sigmay2_max = 5
sigma2_max = 2
sigma2_min = 0.01
rho = 0.6
steps_long = 100
steps_short = 7

sigmax2 = 0.8
sigmay2 = 4
sigma2 = 1 


sigmay2_vals = matrix(seq(0,sigmay2_max,length.out=500),nrow=1)
sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=500),nrow=1)

sigmax2_sigma2 = matrix(NA,nrow=3,ncol=length(sigma2_vals))
sigmax2_sigma2[1,which(sigma2_vals<sigmay2)] = apply(matrix(sigma2_vals[which(sigma2_vals<sigmay2)],nrow=1),2,sigmax2_eq3,sigmax2=sigmax2,sigmay2=sigmay2)
sigmax2_sigma2[2,which(sigma2_vals<3/(5+2*sqrt(6))*sigmay2)] = apply(matrix(sigma2_vals[which(sigma2_vals<3/(5+2*sqrt(6))*sigmay2)],nrow=1),2,sigmax2_eq9,sigmax2=sigmax2,sigmay2=sigmay2)
steps=steps_long
hold = apply(sigma2_vals,2,sigmax2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)
w = which(log(hold)>0)
steps=steps_short
sigmax2_sigma2[3,w] = apply(matrix(sigma2_vals[w],nrow=1),2,sigmax2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)

sigmax2_sigmay2 = matrix(NA,nrow=3,ncol=length(sigmay2_vals))
sigmax2_sigmay2[1,which(sigmay2_vals>sigma2)] = apply(matrix(sigmay2_vals[which(sigmay2_vals>sigma2)],nrow=1),2,sigmax2_eq3,sigmax2=sigmax2,sigma2=sigma2)
sigmax2_sigmay2[2,which(sigmay2_vals>(5+2*sqrt(6))/3*sigma2)] = apply(matrix(sigmay2_vals[which(sigmay2_vals>(5+2*sqrt(6))/3*sigma2)],nrow=1),2,sigmax2_eq9,sigmax2=sigmax2,sigma2=sigma2)
steps=steps_long
hold = apply(sigmay2_vals,2,sigmax2_mode7,sigmax2=sigmax2,sigma2=sigma2)
w = which(log(hold)>0)
steps=steps_short
sigmax2_sigmay2[3,w] = apply(matrix(sigmay2_vals[w],nrow=1),2,sigmax2_mode7,sigmax2=sigmax2,sigma2=sigma2)


width = 6
height = 3
marg = c(0.53,0.5,0.04,0.1)
omarg = c(0.4,0.5 ,0.3,0.0)
ltys = c(4,5,1)

pdf(file='/Users/eleanorbrush/Desktop/equilibrium_variance.pdf',width=width,height=height,family=fontfamily)

par(mfrow=c(1,2),ps=smallfontsize,mai=marg,oma=omarg)

plot(sigmay2_vals,sigmay2_vals,t='n',ylim=range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigmay2_vals,sigmax2_sigmay2[i,],lty=ltys[i])
}
mtext(expression(sigma[y]^2),side=1,at=mean(sigmay2_vals),line=2.2,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2^*),side=2,at=mean(range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE)))),line=1.6,cex=largefontsize/smallfontsize)

legend(-0.3,4,legend=c('','',''),lty=ltys,col='black',bty='n')

plot(sigma2_vals,sigma2_vals,t='n',ylim=range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigma2_vals,sigmax2_sigma2[i,],lty=ltys[i])
}
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2.1,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE)))),line=1.6,cex=largefontsize/smallfontsize)

dev.off()

########## transient variance

sigmax2_mode1 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[1,1,steps])
}

sigmax2_mode2 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[2,1,steps])
}

sigmax2_mode8 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[8,1,steps])
}

sigmax2_max = 2
sigmay2_max = 5
sigma2_max = 2
sigma2_min = 0.01
rho = 0.6
steps_long = 100
steps_short = 10
steps = steps_short

sigmax2 = 0.8
sigmay2 = 4
sigma2 = 1 

sigmax2_vals = matrix(seq(0,sigmax2_max,length.out=500),nrow=1)
sigmay2_vals = matrix(seq(0,sigmay2_max,length.out=500),nrow=1)
sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=500),nrow=1)

sigmax2_sigmay2 = matrix(NA,nrow=4,ncol=length(sigmay2_vals))
sigmax2_sigmay2[1,] = apply(sigmay2_vals,2,sigmax2_mode1,sigmax2=sigmax2,sigma2=sigma2)
sigmax2_sigmay2[2,] = apply(sigmay2_vals,2,sigmax2_mode2,sigmax2=sigmax2,sigma2=sigma2)
sigmax2_sigmay2[3,] = apply(sigmay2_vals,2,sigmax2_mode8,sigmax2=sigmax2,sigma2=sigma2)

sigmax2_sigma2 = matrix(NA,nrow=4,ncol=length(sigma2_vals))
sigmax2_sigma2[1,] = apply(sigma2_vals,2,sigmax2_mode1,sigmax2=sigmax2,sigmay2=sigmay2)
sigmax2_sigma2[2,] = apply(sigma2_vals,2,sigmax2_mode2,sigmax2=sigmax2,sigmay2=sigmay2)
sigmax2_sigma2[3,] = apply(sigma2_vals,2,sigmax2_mode8,sigmax2=sigmax2,sigmay2=sigmay2)

width = 6
height = 3
marg = c(0.53,0.5,0.04,0.1)
omarg = c(0.35,0.45 ,0.3,0.0)
ltys = c(4,1,5,2)
pdf(file='/Users/eleanorbrush/Desktop/transient_variance.pdf',width=width,height=height,family=fontfamily)

par(mfrow=c(1,2),ps=smallfontsize,mai=marg,oma=omarg)

plot(sigmay2_vals,sigmay2_vals,t='n',ylim=range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:3){
	lines(sigmay2_vals,sigmax2_sigmay2[i,],lty=ltys[i])
}
mtext(expression(sigma[y]^2),side=1,at=mean(sigmay2_vals),line=2.1,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)
legend(-0.3,4,legend=c('','','',""),lty=ltys,col='black',bty='n')

plot(sigma2_vals,sigma2_vals,t='n',ylim=range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:3){
	lines(sigma2_vals,sigmax2_sigma2[i,],lty=ltys[i])
}
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

dev.off()