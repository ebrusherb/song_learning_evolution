setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12 
width = 4
height = 4


#################### vector fields of dynamics for various modes

marg = c(0.7,0.7,0.3,0.2)
steps = 1
sigma2 = 1

vectorfield_modified <- function (fun, xlim, ylim, n = 16, scale = 0.05, col = "green", 
    ...) 
{
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.numeric(ylim), 
        length(ylim) == 2)
    xpts <- linspace(xlim[1], xlim[2], n)
    ypts <- linspace(ylim[1], ylim[2], n)
    vectors <- multi.outer(fun,xpts,ypts)
    vectors <-aperm(vectors,c(2,3,1))
    dx = vectors[,,1]
    dy = vectors[,,2]
    plot(xlim, ylim, type = "n")
    grid()
    quiver(x, y, dx, dy, scale = scale, col = col, ...)
}

myvectorfield <- function (fun, xlim, ylim, n = 10, C=100, scale = 0.01, col = "black", length = 0.1, add_contour=FALSE, contour_fun = Q_fun, contour_levels = c(4), contour_n=100,
    ...) 
{
    stopifnot(is.numeric(xlim), length(xlim) == 2, is.numeric(ylim), 
        length(ylim) == 2)
    xpts <- linspace(xlim[1], xlim[2], n)
    ypts <- linspace(ylim[1], ylim[2], n)
    vectors <- multi.outer(fun,xpts,ypts)
    vectors <-C*aperm(vectors,c(3,2,1))
    M <- meshgrid(xpts,ypts)
    x = M$X
    y = M$Y
    dx = vectors[,,1]
    dy = vectors[,,2]
    plot(xlim, ylim, type = "n",xlab="",ylab="")
    grid()
    quiver(x, y, dx, dy, scale = scale , length = length, col = col, ...)
    if(add_contour){
    	xpts <- linspace(xlim[1], xlim[2], contour_n)
	    ypts <- linspace(ylim[1], ylim[2], contour_n)
    	contour_mat = outer(xpts,ypts,FUN=contour_fun)
    	contour(xpts,ypts,contour_mat,levels=contour_levels,add=TRUE,col='red',labels=" ")
    }
}


mode2 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[2,1:2,2]-r[2,1:2,1])
}

mode7 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[7,1:2,2]-r[7,1:2,1])
} 

mode8 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[8,1:2,2]-r[8,1:2,1])
} 

mode9 <- function(sigmax2,cov){
	if(sigmax2*sigmay2!=0){
		rho = cov/sqrt(sigmax2*sigmay2) } 
		else { rho = 0 }
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[9,c(1,3),2]-r[9,c(1,3),1])
} 

solveQ <-function(sigmax2){
	p = Re(polyroot(c((-2*sigma2-3*sigmax2)*(sigma2+sigmax2),2*sqrt(sigmax2)*(sigma2+sigmax2),sigmax2)))
	p = p[p>0]
	p^2
}



# # # # # mode 7 
pdf(file='/Users/eleanorbrush/Desktop/mode7_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1
sigmax2_pts = unique(c(seq(0,0.01,length.out=30),linspace(sigmax2_lim[1],sigmax2_lim[2],150)))
Q4 = apply(matrix(sigmax2_pts,nrow=1),2,solveQ)

par(ps=smallfontsize,mai=marg)
myvectorfield(mode7,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)
lines(sigmax2_pts,Q4,col='red')
points(0,0,pch=20,cex=1,col='red')

sigmax2_toplot = c(2,4,1)
sigmay2_toplot = c(3,2,3.45)
steps = 4
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all(sigmax2,sigmay2,sigma2,rho)
	lines(r[7,1,],r[7,2,])	
	k=1
	arrows(r[7,1,steps-k],r[7,2,steps-k],r[7,1,steps+1],r[7,2,steps+1],length=0.05)
}

dev.off()

# # # # # mode 9 
pdf(file='/Users/eleanorbrush/Desktop/mode9_dynamics.pdf',width=width,height=height,family=fontfamily)

sigmay2 = 3.4
sigmax2_lim = c(0,5)
cov_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode9,sigmax2_lim,cov_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
points(0,0,pch=20,cex=1,col='red')
sigmax2_eq = (3*sigmay2-5*sigma2+sqrt(sigma2^2-30*sigma2*sigmay2+9*sigmay2^2))/6
cov_eq = sigmax2_eq*sigmay2/(sigma2+sigmax2_eq)
points(sigmax2_eq,cov_eq,pch=20,cex=1,col='red')
sigmax2_eq = (3*sigmay2-5*sigma2-sqrt(sigma2^2-30*sigma2*sigmay2+9*sigmay2^2))/6
cov_eq = sigmax2_eq*sigmay2/(sigma2+sigmax2_eq)
points(sigmax2_eq,cov_eq,pch=1,cex=1,col='red')

sigmax2_toplot = c(1,1.65,2)
cov_toplot = c(0.03,0.05,4)
steps = 25
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	cov =cov_toplot[i]
	rho = cov/sqrt(sigmax2*sigmay2)
	r = recursion_all(sigmax2,sigmay2,sigma2,rho)
	lines(r[9,1,],r[9,3,])	
	k=1
	arrows(r[9,1,steps-k],r[9,3,steps-k],r[9,1,steps+1],r[9,3,steps+1],length=0.05)
}

mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext('Cov',side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 2 
pdf(file='/Users/eleanorbrush/Desktop/mode2_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode2,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2 = 4
sigmay2 = 2.05
steps = 20 
r = recursion_all(sigmax2,sigmay2,sigma2,rho)
lines(r[2,1,],r[2,2,])
k = 3
arrows(r[2,1,steps-k],r[2,2,steps-k],r[2,1,steps],r[2,2,steps],length=0.05)
points(0,0,pch=20,cex=1,col='red')
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 8
pdf(file='/Users/eleanorbrush/Desktop/mode8_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode8,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2_toplot = c(4,1)
sigmay2_toplot = c(2,3.45)
steps = 4
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all(sigmax2,sigmay2,sigma2,rho)
	lines(r[8,1,],r[8,2,])	
	k=0
	arrows(r[8,1,steps-k],r[8,2,steps-k],r[8,1,steps+1],r[8,2,steps+1],length=0.05)
}
points(0,0,pch=20,cex=1,col='red')
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

################ equilibrium variance of songs and traits for four modes
marg = c(0.53,0.5,0.04,0.1)
omarg = c(0.4,0.4 ,0.3,0.0)


sigmax2_eq3 <- function(sigmax2,sigmay2,sigma2){
	sigmay2-sigma2
}

sigmax2_eq5 <-function(sigmax2,sigmay2,sigma2){
	sigmax2
}

sigmax2_mode7 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[7,1,steps])
}

sigmax2_eq9 <- function(sigmax2,sigmay2,sigma2){
	(3*sigmay2-5*sigma2+sqrt(9*sigmay2^2-30*sigma2*sigmay2+sigma2^2))/6
}

sigmay2_eq3<- function(sigmax2,sigmay2,sigma2){
	sigmay2
}

sigmay2_eq5 <-function(sigmax2,sigmay2,sigma2){
	sigmax2*(sigma2+sigmax2)/(2*sigmax2+sigma2)
}

sigmay2_mode7 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[7,2,steps+1])
}

sigmay2_eq9 <- function(sigmax2,sigmay2,sigma2){
	sigmay2
}

sigmax2_max = 2
sigmay2_max = 5
sigma2_max = 2
sigma2_min = 0.01
rho = 0.6
steps_long = 100
steps_short = 7

sigmax2 = 0.8
sigmay2 = 4
sigma2 = 1 

sigmax2_vals = matrix(seq(0,sigmax2_max,length.out=500),nrow=1)
sigmay2_vals = matrix(seq(0,sigmay2_max,length.out=500),nrow=1)
sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=500),nrow=1)

sigmax2_sigmay2 = matrix(NA,nrow=4,ncol=length(sigmay2_vals))
sigmax2_sigmay2[1,which(sigmay2_vals>sigma2)] = apply(matrix(sigmay2_vals[which(sigmay2_vals>sigma2)],nrow=1),2,sigmax2_eq3,sigmax2=sigmax2,sigma2=sigma2)
sigmax2_sigmay2[2,] = apply(sigmay2_vals,2,sigmax2_eq5,sigmax2=sigmax2,sigma2=sigma2)
sigmax2_sigmay2[4,which(sigmay2_vals>(5+2*sqrt(6))/3*sigma2)] = apply(matrix(sigmay2_vals[which(sigmay2_vals>(5+2*sqrt(6))/3*sigma2)],nrow=1),2,sigmax2_eq9,sigmax2=sigmax2,sigma2=sigma2)
steps=steps_long
hold = apply(sigmay2_vals,2,sigmax2_mode7,sigmax2=sigmax2,sigma2=sigma2)
w = which(log(hold)>0)
steps=steps_short
sigmax2_sigmay2[3,w] = apply(matrix(sigmay2_vals[w],nrow=1),2,sigmax2_mode7,sigmax2=sigmax2,sigma2=sigma2)

sigmax2_sigma2 = matrix(NA,nrow=4,ncol=length(sigma2_vals))
sigmax2_sigma2[1,which(sigma2_vals<sigmay2)] = apply(matrix(sigma2_vals[which(sigma2_vals<sigmay2)],nrow=1),2,sigmax2_eq3,sigmax2=sigmax2,sigmay2=sigmay2)
sigmax2_sigma2[2,]=apply(sigma2_vals,2,sigmax2_eq5,sigmax2=sigmax2,sigmay2=sigmay2)
sigmax2_sigma2[4,which(sigma2_vals<3/(5+2*sqrt(6))*sigmay2)] = apply(matrix(sigma2_vals[which(sigma2_vals<3/(5+2*sqrt(6))*sigmay2)],nrow=1),2,sigmax2_eq9,sigmax2=sigmax2,sigmay2=sigmay2)
steps=steps_long
hold = apply(sigma2_vals,2,sigmax2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)
w = which(log(hold)>0)
steps=steps_short
sigmax2_sigma2[3,w] = apply(matrix(sigma2_vals[w],nrow=1),2,sigmax2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)

sigmay2_sigmax2 = matrix(NA,nrow=4,ncol=length(sigmax2_vals))
sigmay2_sigmax2[1,] = apply(sigmax2_vals,2,sigmay2_eq3,sigmay2=sigmay2,sigma2=sigma2)
sigmay2_sigmax2[2,] = apply(sigmax2_vals,2,sigmay2_eq5,sigmay2=sigmay2,sigma2=sigma2)
sigmay2_sigmax2[4,] = apply(sigmax2_vals,2,sigmay2_eq9,sigmay2=sigmay2,sigma2=sigma2)
steps=steps_long
hold = apply(sigmax2_vals,2,sigmay2_mode7,sigma2=sigma2,sigmay2=sigmay2)
w = which(log(hold)>0)
steps =steps_short
sigmay2_sigmax2[3,w] = apply(matrix(sigmax2_vals[w],nrow=1),2,sigmay2_mode7,sigma2=sigma2,sigmay2=sigmay2)

sigmay2_sigma2 = matrix(NA,nrow=4,ncol=length(sigma2_vals))
sigmay2_sigma2[1,] = apply(sigma2_vals,2,sigmay2_eq3,sigmax2=sigmax2,sigmay2=sigmay2)
sigmay2_sigma2[2,] = apply(sigma2_vals,2,sigmay2_eq5,sigmax2=sigmax2,sigmay2=sigmay2)
sigmay2_sigma2[4,] = apply(sigma2_vals,2,sigmay2_eq9,sigmax2=sigmax2,sigmay2=sigmay2)
steps = steps_long
hold = apply(sigma2_vals,2,sigmay2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)
steps = steps_short
w = which(log(hold)>0)
sigmay2_sigma2[3,w] = apply(matrix(sigma2_vals[w],nrow=1),2,sigmay2_mode7,sigmax2=sigmax2,sigmay2=sigmay2)

ltys = c(4,5,1,2)
pdf(file='/Users/eleanorbrush/Desktop/equilibrium_variance_full.pdf',width=width,height=height,family=fontfamily)

par(mfrow=c(2,2),ps=smallfontsize,mai=marg,oma=omarg)

plot(sigmay2_vals,sigmay2_vals,t='n',ylim=range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigmay2_vals,sigmax2_sigmay2[i,],lty=ltys[i])
}
lines(sigmay2_vals,sigmax2_sigmay2[4,],lty=ltys[4])
mtext(expression(sigma[y]^2),side=1,at=mean(sigmay2_vals),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE)))),line=1.7,cex=largefontsize/smallfontsize)

legend(-0.3,4.5,legend=c('','','',""),lty=ltys,col='black',bty='n')


plot(sigma2_vals,sigma2_vals,t='n',ylim=range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigma2_vals,sigmax2_sigma2[i,],lty=ltys[i])
}
lines(sigma2_vals,sigmax2_sigma2[4,],lty=ltys[4])
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2.3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(c(range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmay2,na.rm=TRUE)))),line=1.6,cex=largefontsize/smallfontsize)

plot(sigmax2_vals,sigmax2_vals,t='n',ylim=range(c(range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigmax2_vals,sigmay2_sigmax2[i,],lty=ltys[i])
}
lines(sigmax2_vals,sigmay2_sigmax2[4,],lty=ltys[4])
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_vals),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(range(c(range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE)))),line=1.7,cex=largefontsize/smallfontsize)


plot(sigma2_vals,sigma2_vals,t='n',ylim=range(c(range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE))),xlab='',ylab='')
for(i in 1:3){
	lines(sigma2_vals,sigmay2_sigma2[i,],lty=ltys[i])
}
lines(sigma2_vals,sigmay2_sigma2[4,],lty=ltys[4])
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2.3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(range(c(range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE)))),line=1.6,cex=largefontsize/smallfontsize)

dev.off()

############### transient variance

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

sigmay2_mode1 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[1,2,steps])
}

sigmay2_mode2 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[2,2,steps])
}

sigmay2_mode8 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[8,2,steps])
}

sigmay2_mode4 <-function(sigmax2,sigmay2,sigma2){
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[4,2,steps])
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

sigmax2_sigmax2 = matrix(NA,nrow=4,ncol=length(sigmax2_vals))
sigmax2_sigmax2[1,] = apply(sigmax2_vals,2,sigmax2_mode1,sigmay2=sigmay2,sigma2=sigma2)
sigmax2_sigmax2[2,] = apply(sigmax2_vals,2,sigmax2_mode2,sigmay2=sigmay2,sigma2=sigma2)
sigmax2_sigmax2[3,] = apply(sigmax2_vals,2,sigmax2_mode8,sigmay2=sigmay2,sigma2=sigma2)

sigmay2_sigmay2 = matrix(NA,nrow=4,ncol=length(sigmay2_vals))
sigmay2_sigmay2[1,] = apply(sigmay2_vals,2,sigmay2_mode1,sigmax2=sigmax2,sigma2=sigma2)
sigmay2_sigmay2[2,] = apply(sigmay2_vals,2,sigmay2_mode2,sigmax2=sigmax2,sigma2=sigma2)
sigmay2_sigmay2[3,] = apply(sigmay2_vals,2,sigmay2_mode8,sigmax2=sigmax2,sigma2=sigma2)
sigmay2_sigmay2[4,] = apply(sigmay2_vals,2,sigmay2_mode4,sigmax2=sigmax2,sigma2=sigma2)

sigmay2_sigma2 = matrix(NA,nrow=4,ncol=length(sigma2_vals))
sigmay2_sigma2[1,] = apply(sigma2_vals,2,sigmay2_mode1,sigmax2=sigmax2,sigmay2=sigmay2)
sigmay2_sigma2[2,] = apply(sigma2_vals,2,sigmay2_mode2,sigmax2=sigmax2,sigmay2=sigmay2)
sigmay2_sigma2[3,] = apply(sigma2_vals,2,sigmay2_mode8,sigmax2=sigmax2,sigmay2=sigmay2)
sigmay2_sigma2[4,] = apply(sigma2_vals,2,sigmay2_mode4,sigmax2=sigmax2,sigmay2=sigmay2)

sigmay2_sigmax2 = matrix(NA,nrow=4,ncol=length(sigmax2_vals))
sigmay2_sigmax2[1,] = apply(sigmax2_vals,2,sigmay2_mode1,sigmay2=sigmay2,sigma2=sigma2)
sigmay2_sigmax2[2,] = apply(sigmax2_vals,2,sigmay2_mode2,sigmay2=sigmay2,sigma2=sigma2)
sigmay2_sigmax2[3,] = apply(sigmax2_vals,2,sigmay2_mode8,sigmay2=sigmay2,sigma2=sigma2)
sigmay2_sigmax2[4,] = apply(sigmax2_vals,2,sigmay2_mode4,sigmay2=sigmay2,sigma2=sigma2)

width = 6
height = 4
marg = c(0.53,0.5,0.04,0.1)
omarg = c(0,0.24 ,0.3,0.0)
ltys = c(4,1,5,2)
pdf(file='/Users/eleanorbrush/Desktop/transient_variance_full.pdf',width=width,height=height,family=fontfamily)

par(mfrow=c(2,3),ps=smallfontsize,mai=marg,oma=omarg)

plot(sigmay2_vals,sigmay2_vals,t='n',ylim=range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:3){
	lines(sigmay2_vals,sigmax2_sigmay2[i,],lty=ltys[i])
}
mtext(expression(sigma[y]^2),side=1,at=mean(sigmay2_vals),line=3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)
legend(-0.3,4,legend=c('','','',""),lty=ltys,col='black',bty='n')

plot(sigma2_vals,sigma2_vals,t='n',ylim=range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:3){
	lines(sigma2_vals,sigmax2_sigma2[i,],lty=ltys[i])
}
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

plot(sigmax2_vals,sigmax2_vals,t='n',ylim=range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:3){
	lines(sigmax2_vals,sigmax2_sigmax2[i,],lty=ltys[i])
}
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_vals),line=3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[x]^2),side=2,at=mean(range(range(sigmax2_sigmay2,na.rm=TRUE),range(sigmax2_sigma2,na.rm=TRUE),range(sigmax2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

plot(sigmay2_vals,sigmay2_vals,t='n',ylim=range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:4){
	lines(sigmay2_vals,sigmay2_sigmay2[i,],lty=ltys[i])
}
mtext(expression(sigma[y]^2),side=1,at=mean(sigmay2_vals),line=3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

plot(sigma2_vals,sigma2_vals,t='n',ylim=range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:4){
	lines(sigma2_vals,sigmay2_sigma2[i,],lty=ltys[i])
}
mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

plot(sigmax2_vals,sigmax2_vals,t='n',ylim=range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE)),xlab='',ylab='')
for(i in 1:4){
	lines(sigmax2_vals,sigmay2_sigmax2[i,],lty=ltys[i])
}
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_vals),line=3,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(range(range(sigmay2_sigmay2,na.rm=TRUE),range(sigmay2_sigma2,na.rm=TRUE),range(sigmay2_sigmax2,na.rm=TRUE))),line=1.6,cex=largefontsize/smallfontsize)

dev.off()