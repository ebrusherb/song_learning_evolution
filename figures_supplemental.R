setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
library(RColorBrewer)
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12 
width = 3.5
height = 3.5


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

mode1 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all(sigmax2,sigmay2,sigma2,rho)
	return(r[1,1:2,2]-r[1,1:2,1])
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
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)
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

mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext('Covariance, Cov',side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

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
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

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
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 1
pdf(file='/Users/eleanorbrush/Desktop/mode1_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode1,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2_toplot = c(4,0.5)
sigmay2_toplot = c(2,3.45)
steps = 10
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all(sigmax2,sigmay2,sigma2,rho)
	lines(r[1,1,],r[1,2,])	
	k=0
	arrows(r[1,1,steps-k],r[1,2,steps-k],r[1,1,steps+1],r[1,2,steps+1],length=0.05)
}
points(0,0,pch=20,cex=1,col='red')
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

# #######  equilibrium and transient variance

recursion_all_end <- function(sigmax2,sigmay2,sigma2,equil=TRUE){
	r = recursion_all(sigmax2,sigmay2,sigma2,rho)
	toreturn = r[,1:2,steps+1]
	if(equil){		
			if(sigma2>sigmay2){
				toreturn[3,1:2] = NA
			}else{toreturn[3,1:2] = c(sigmay2-sigma2,sigmay2)
				}
			if(diff(r[9,1,steps+(0:1)])<0 || r[9,1,steps]<1e-10){
				toreturn[9,1:2] = NA
			}
		} else{
			if(sigma2<sigmay2){
				toreturn[3,1:2] = NA
			} 
			if(diff(r[7,1,steps+(0:1)])>0){
				toreturn[7,1:2] = NA
			} 
			if(diff(r[9,1,steps+(0:1)])>0){
				toreturn[9,1:2] = NA
			}
			}
	return(toreturn)
}

time_to_zero <- function(sigmax2,sigmay2,sigma2){
	r  = recursion_all(sigmax2,sigmay2,sigma2,rho)
	time = array(NA,c(9,2))
	for(i in 1:9){
		for(j in 1:2){
			w = which(r[i,j,]<disappear_thresh)
			if(length(w)==0){
				time[i,j] = NA
			} else {
				time[i,j] = log(w[1])
			}
		}
	}
	return(time)
}


sigmax2_max = 2
sigmax2_min = 0.01
sigmay2_max = 4.5
sigmay2_min = 0.01
sigma2_max = 4.5
sigma2_min = 0.01
rho = 0.6
steps_long = 1000
steps_short = 10
steps_time = 5000
disappear_thresh <- 0.05 #5e-3

sigmax2 = 1
sigmay2 = 4.1
sigma2 = 1

vec_length = 100

#### sigmax2 as function of sigmax2

sigmax2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmax2 <- as.list(1:3)

sigmax2_pos = c(3,9)
sigmax2_zero = c(1,2,3,4,7,8,9)

sigmay2_pos = c(5)
sigmay2_zero = c(1,2,4,7,8)

steps = steps_long

r_sigmax2[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigma2=sigma2,sigmay2=sigmay2,simplify=FALSE)
r_sigmax2[[1]] <- matrix(unlist(lapply(r_sigmax2[[1]],function(y) y[,1])),nrow=9)
r_sigmax2[[1]][setdiff(1:9,sigmax2_pos),] = NA

steps = steps_short

r_sigmax2[[2]] <- sapply(sigmax2_vals,recursion_all_end,sigma2=sigma2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
r_sigmax2[[2]] <- matrix(unlist(lapply(r_sigmax2[[2]],function(y) y[,1])),nrow=9)
r_sigmax2[[2]][setdiff(1:9,sigmax2_zero),] = NA

steps = steps_time 

r_sigmax2[[3]] <- sapply(sigmax2_vals,time_to_zero,sigma2=sigma2,sigmay2=sigmay2,simplify=FALSE)
r_sigmax2[[3]] <- matrix(unlist(lapply(r_sigmax2[[3]],function(y) y[,1])),nrow=9)
r_sigmax2[[3]][setdiff(1:9,sigmax2_zero),] = NA

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.5
height =6

pdf('/Users/eleanorbrush/Desktop/sigmax2_sigmax2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:3,7:9)[order(old_to_new[c(1:3,7:9)])],nrow=3,byrow=TRUE)

for(k in 1:3){
	subset = modes[k,]
	
	ylim = range(r_sigmax2[1:2],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmax2[[1]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmax2[[2]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial song variance, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Song variance, ',sigma[x]^2)),side=2,at=mean(ylim)-0.1*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	ylim = range(r_sigmax2[3],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmax2[[3]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(paste('Initial song variance, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(0.6,2.6,legend=c(expression(paste(sigma[x]^2, "*",', song genetic')),expression(paste(sigma[x]^2~(10),', song genetic')),expression(paste(sigma[x]^2, "*",', song paternally learned')),expression(paste(sigma[x]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}
}


dev.off()

#### sigmax2 as function of sigmay2

sigmay2_vals = matrix(seq(sigmay2_min,sigmay2_max,length.out=vec_length),nrow=1)

r_sigmay2 <- as.list(1:3)

sigmax2_pos = c(3,9)
sigmax2_zero = c(1,2,3,4,7,8,9)

sigmay2_pos = c(5)
sigmay2_zero = c(1,2,4,7,8)

steps = steps_long

r_sigmay2[[1]] <- sapply(sigmay2_vals,recursion_all_end,sigma2=sigma2,sigmax2=sigmax2,simplify=FALSE)
r_sigmay2[[1]] <- matrix(unlist(lapply(r_sigmay2[[1]],function(y) y[,1])),nrow=9)
r_sigmay2[[1]][setdiff(1:9,sigmax2_pos),] = NA

steps = steps_short

r_sigmay2[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigma2=sigma2,sigmax2=sigmax2,equil=FALSE,simplify=FALSE)
r_sigmay2[[2]] <- matrix(unlist(lapply(r_sigmay2[[2]],function(y) y[,1])),nrow=9)
r_sigmay2[[2]][setdiff(1:9,sigmax2_zero),] = NA

steps = steps_time 

r_sigmay2[[3]] <- sapply(sigmay2_vals,time_to_zero,sigma2=sigma2,sigmax2=sigmax2,simplify=FALSE)
r_sigmay2[[3]] <- matrix(unlist(lapply(r_sigmay2[[3]],function(y) y[,1])),nrow=9)
r_sigmay2[[3]][setdiff(1:9,sigmax2_zero),] = NA

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.5
height =6

pdf('/Users/eleanorbrush/Desktop/sigmax2_sigmay2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:3,7:9)[order(old_to_new[c(1:3,7:9)])],nrow=3,byrow=TRUE)

for(k in 1:3){
	subset = modes[k,]
	
	ylim = range(r_sigmay2[1:2],na.rm=TRUE)
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmay2[[1]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmay2[[2]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial preference variance, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Song variance, ',sigma[x]^2)),side=2,at=mean(ylim)-0.1*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(-0.2,3.7,legend=c(expression(paste(sigma[x]^2, "*",', song genetic')),expression(paste(sigma[x]^2~(10),', song genetic')),expression(paste(sigma[x]^2, "*",', song paternally learned')),expression(paste(sigma[x]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}
	
	ylim = range(r_sigmay2[3],na.rm=TRUE)
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmay2[[3]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(paste('Initial preference variance, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
}


dev.off() 

#### sigmay2 as a function of sigma2

sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=vec_length),nrow=1)

r_sigma2 <- as.list(1:3)

sigmax2_pos = c(3,9)
sigmax2_zero = c(1,2,3,4,7,8,9)

sigmay2_pos = c(5)
sigmay2_zero = c(1,2,4,7,8)

steps = steps_long

r_sigma2[[1]] <- sapply(sigma2_vals,recursion_all_end,sigmay2=sigmay2,sigmax2=sigmax2,simplify=FALSE)
r_sigma2[[1]] <- matrix(unlist(lapply(r_sigma2[[1]],function(y) y[,2])),nrow=9)
r_sigma2[[1]][setdiff(1:9,sigmay2_pos),] = NA

steps = steps_short

r_sigma2[[2]] <- sapply(sigma2_vals,recursion_all_end,sigmay2=sigmay2,sigmax2=sigmax2,equil=FALSE,simplify=FALSE)
r_sigma2[[2]] <- matrix(unlist(lapply(r_sigma2[[2]],function(y) y[,2])),nrow=9)
r_sigma2[[2]][setdiff(1:9,sigmay2_zero),] = NA

steps = steps_time 

r_sigma2[[3]] <- sapply(sigma2_vals,time_to_zero,sigmay2=sigmay2,sigmax2=sigmax2,simplify=FALSE)
r_sigma2[[3]] <- matrix(unlist(lapply(r_sigma2[[3]],function(y) y[,2])),nrow=9)
r_sigma2[[3]][setdiff(1:9,sigmay2_zero),] = NA

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'Greens')[4],brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Greens')[7],brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.5
height =4

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigma2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:4,nrow=2,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:2,4:5,7:8)[order(old_to_new[c(1:2,4:5,7:8)])],nrow=2,byrow=TRUE)

for(k in 1:2){
	subset = modes[k,]
	
	ylim = range(r_sigma2[1:2],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigma2[[1]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigma2[[2]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Variance of pref. function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.2,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(ylim)-0.23*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	
	ylim = range(r_sigma2[3],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigma2[[3]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(paste('Variance of pref. function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.3,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(0.2,6.2,legend=c(expression(paste(sigma[y]^2,"*",', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song genetic')),expression(paste(sigma[y]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],dark_col[2],dark_col[3]),bty='n')	
		}
}


dev.off()

#### sigmay2 as a function of sigmax2


sigmax2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmax2 <- as.list(1:3)

sigmax2_pos = c(3,9)
sigmax2_zero = c(1,2,3,4,7,8,9)

sigmay2_pos = c(5)
sigmay2_zero = c(1,2,4,7,8)

steps = steps_long

r_sigmax2[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
r_sigmax2[[1]] <- matrix(unlist(lapply(r_sigmax2[[1]],function(y) y[,2])),nrow=9)
r_sigmax2[[1]][setdiff(1:9,sigmay2_pos),] = NA

steps = steps_short

r_sigmax2[[2]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
r_sigmax2[[2]] <- matrix(unlist(lapply(r_sigmax2[[2]],function(y) y[,2])),nrow=9)
r_sigmax2[[2]][setdiff(1:9,sigmay2_zero),] = NA

steps = steps_time 

r_sigmax2[[3]] <- sapply(sigmax2_vals,time_to_zero,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
r_sigmax2[[3]] <- matrix(unlist(lapply(r_sigmax2[[3]],function(y) y[,2])),nrow=9)
r_sigmax2[[3]][setdiff(1:9,sigmay2_zero),] = NA

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'Greens')[4],brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Greens')[7],brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.5
height =4

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigmax2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:4,nrow=2,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:2,4:5,7:8)[order(old_to_new[c(1:2,4:5,7:8)])],nrow=2,byrow=TRUE)

for(k in 1:2){
	subset = modes[k,]
	
	ylim = range(r_sigmax2[1:2],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmax2[[1]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmax2[[2]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial song variance, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(ylim)-0.23*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==2){
		legend(-0.1,2.9,legend=c(expression(paste(sigma[y]^2,"*",', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song genetic')),expression(paste(sigma[y]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],dark_col[2],dark_col[3]),bty='n')	
		}
	
	ylim = range(r_sigmax2[3],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmax2[[3]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(paste('Initial song variance, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
		
}


dev.off()


#### sigmay2 as a function of sigmay2

sigmay2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmay2 <- as.list(1:3)

sigmax2_pos = c(3,9)
sigmax2_zero = c(1,2,3,4,7,8,9)

sigmay2_pos = c(5)
sigmay2_zero = c(1,2,4,7,8)

steps = steps_long

r_sigmay2[[1]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
r_sigmay2[[1]] <- matrix(unlist(lapply(r_sigmay2[[1]],function(y) y[,2])),nrow=9)
r_sigmay2[[1]][setdiff(1:9,sigmay2_pos),] = NA

steps = steps_short

r_sigmay2[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
r_sigmay2[[2]] <- matrix(unlist(lapply(r_sigmay2[[2]],function(y) y[,2])),nrow=9)
r_sigmay2[[2]][setdiff(1:9,sigmay2_zero),] = NA

steps = steps_time 

r_sigmay2[[3]] <- sapply(sigmay2_vals,time_to_zero,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
r_sigmay2[[3]] <- matrix(unlist(lapply(r_sigmay2[[3]],function(y) y[,2])),nrow=9)
r_sigmay2[[3]][setdiff(1:9,sigmay2_zero),] = NA

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'Greens')[4],brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Greens')[7],brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.5
height =4

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigmay2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:4,nrow=2,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:2,4:5,7:8)[order(old_to_new[c(1:2,4:5,7:8)])],nrow=2,byrow=TRUE)

for(k in 1:2){
	subset = modes[k,]
	
	ylim = range(r_sigmay2[1:2],na.rm=TRUE)
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmay2[[1]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmay2[[2]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial preference variance, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=2.6,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(ylim)-0.23*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(0.2,0.7,legend=c(expression(paste(sigma[y]^2,"*",', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song genetic')),expression(paste(sigma[y]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],dark_col[2],dark_col[3]),bty='n')	
		}
	
	ylim = range(r_sigmay2[3],na.rm=TRUE)
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmay2[[3]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(paste('Initial preference variance, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=2.6,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
		
}


dev.off()

#### effect of step distributions on eq distribution of preferences
col_vec = c(brewer.pal(9,'Set1')[c(1,7,2:4)],brewer.pal(8,'Dark2')[c(1)],brewer.pal(9,'Set1')[c(9,5,8)])
lwd = 2
marg = c(0.45,0.43,0.02,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 5.625

pdf('/Users/eleanorbrush/Desktop/effect_of_step_distribution_on_preferences.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=2))

load('step_song_equilibrium.Rdata')

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

xlim = c(-6,6)
ylim=c(0,.13)
plot(mrange+1,m_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song, x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
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
		
		if(length(which(p<0))==0 && p[3]>=p[2] ){

			m_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
			m_init = c(m_init,rev(m_init[1:(length(m_init)-1)]))
			
			lines(mrange+1,m_init,col=col_vec[j+1],lwd=lwd)
		}
}

for(i in 2){
	plot(sigma2_vals,var_mat[i,,1],ylim=c(0,round(max(var_mat[i,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='')
	mtext(expression(paste('Variance of pref. function, ', sigma^2)),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste('Eq. song var., ',sigma[x]^2, "*")),expression(paste('Eq. pref. var., ',sigma[y]^2, "*")))[[i]],side=2,line=1.7,at=(round(max(var_mat[i,,],na.rm=TRUE),1)+0.1)/2,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+1],t='o',lwd=lwd)}
	}
}

xlim = c(-6,6)
ylim=c(0,.16)
k=5
plot(mrange+1,equilibrium[[1,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[2,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[2,k,j+1]],col=col_vec[j+1],lwd=lwd)
		}
}

load('step_pref_equilibrium.Rdata')

f_init = dnorm(mrange,mmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)

xlim = c(-6,6)
ylim=c(0,.13)
plot(mrange+1,f_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
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
		
		if(length(which(p<0))==0 && p[3]>=p[2] ){

			f_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
			f_init = c(f_init,rev(f_init[1:(length(f_init)-1)]))
			
			lines(mrange+1,f_init,col=col_vec[j+5],lwd=lwd)
		}
}

for(i in 2){
	plot(sigma2_vals,var_mat[i,,1],ylim=c(0,round(max(var_mat[i,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='')
	mtext(expression(paste('Variance of pref. function, ', sigma^2)),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste('Eq. song var., ',sigma[x]^2, "*")),expression(paste('Eq. pref. var., ',sigma[y]^2, "*")))[[i]],side=2,line=1.7,at=(round(max(var_mat[i,,],na.rm=TRUE),1)+0.1)/2,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+5],t='o',lwd=lwd)}
	}
}

xlim = c(-6,6)
ylim=c(0,.16)
k=5
plot(mrange+1,equilibrium[[2,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[2,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[2,k,j+1]],col=col_vec[j+5],lwd=lwd)
		}
}

dev.off()



# source('step_function_example.R')
# sourece('mutation_test.R')
# source('pearson_kurtosis.R)