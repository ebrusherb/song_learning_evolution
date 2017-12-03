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

myred = brewer.pal(9,'Set1')[1]

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
    	contour(xpts,ypts,contour_mat,levels=contour_levels,add=TRUE,col=myred,labels=" ")
    }
}



mode2 <- function(sigmax2,cov){
	if(sigmax2*sigmay2!=0){
		rho = cov/sqrt(sigmax2*sigmay2) } 
		else { rho = 0 }
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[2,c(1,3),2]-r[2,c(1,3),1])
} 

mode5 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[5,1:2,2]-r[5,1:2,1])
} 


mode6 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[6,1:2,2]-r[6,1:2,1])
}

mode8 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[8,1:2,2]-r[8,1:2,1])
} 

mode9 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[9,1:2,2]-r[9,1:2,1])
}

mode11 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[11,1:2,2]-r[11,1:2,1])
} 

mode12 <- function(sigmax2,sigmay2){
	cov = rho*sqrt(sigmax2*sigmay2)
	r =recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	return(r[12,1:2,2]-r[12,1:2,1])
}

solveQ <-function(sigmax2){
	p = Re(polyroot(c((-2*sigma2-3*sigmax2)*(sigma2+sigmax2),2*sqrt(sigmax2)*(sigma2+sigmax2),sigmax2)))
	p = p[p>0]
	p^2
}

# # # # # mode 2
pdf(file='/Users/eleanorbrush/Desktop/mode2_dynamics.pdf',width=width,height=height,family=fontfamily)

sigmay2 = 3.4
sigmax2_lim = c(0,5)
cov_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode2,sigmax2_lim,cov_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
points(0,0,pch=20,cex=1,col=myred)
sigmax2_eq = (3*sigmay2-5*sigma2+sqrt(sigma2^2-30*sigma2*sigmay2+9*sigmay2^2))/6
cov_eq = sigmax2_eq*sigmay2/(sigma2+sigmax2_eq)
points(sigmax2_eq,cov_eq,pch=20,cex=1,col=myred)
sigmax2_eq = (3*sigmay2-5*sigma2-sqrt(sigma2^2-30*sigma2*sigmay2+9*sigmay2^2))/6
cov_eq = sigmax2_eq*sigmay2/(sigma2+sigmax2_eq)
points(sigmax2_eq,cov_eq,pch=1,cex=1,col=myred)

sigmax2_toplot = c(1,1.65,2)
cov_toplot = c(0.03,0.05,4)
steps = 25
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	cov =cov_toplot[i]
	rho = cov/sqrt(sigmax2*sigmay2)
	r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	lines(r[2,1,],r[2,3,])	
	k=1
	arrows(r[2,1,steps-k],r[2,3,steps-k],r[2,1,steps+1],r[2,3,steps+1],length=0.05)
}

mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext('Covariance, C',side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

# # # # # mode 5 
pdf(file='/Users/eleanorbrush/Desktop/mode5_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1
sigmax2_pts = unique(c(seq(0,0.01,length.out=30),linspace(sigmax2_lim[1],sigmax2_lim[2],150)))
Q4 = apply(matrix(sigmax2_pts,nrow=1),2,solveQ)

par(ps=smallfontsize,mai=marg)
myvectorfield(mode5,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)
lines(sigmax2_pts,Q4,col=myred)
points(0,0,pch=20,cex=1,col=myred)

sigmax2_toplot = c(2,4,1)
sigmay2_toplot = c(3,2,3.45)
steps = 4
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	lines(r[5,1,],r[5,2,])	
	k=1
	arrows(r[5,1,steps-k],r[5,2,steps-k],r[5,1,steps+1],r[5,2,steps+1],length=0.05)
}

dev.off()



#### mode 6
pdf(file='/Users/eleanorbrush/Desktop/mode6_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode6,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2_toplot = c(4,0.5)
sigmay2_toplot = c(2,3.45)
steps = 10
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	lines(r[6,1,],r[6,2,])	
	k=0
	arrows(r[6,1,steps-k],r[6,2,steps-k],r[6,1,steps+1],r[6,2,steps+1],length=0.05)
}
points(0,0,pch=20,cex=1,col=myred)
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
points(0,0,pch=20,cex=1,col=myred)
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 9 
pdf(file='/Users/eleanorbrush/Desktop/mode9_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode9,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2 = 4
sigmay2 = 2.05
steps = 20 
r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
lines(r[9,1,],r[9,2,])
k = 3
arrows(r[9,1,steps-k],r[9,2,steps-k],r[9,1,steps],r[9,2,steps],length=0.05)
points(0,0,pch=20,cex=1,col=myred)
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 11
pdf(file='/Users/eleanorbrush/Desktop/mode11_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode11,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2_toplot = c(4)
sigmay2_toplot = c(2)
steps = 10
for(i in 1:length(sigmax2_toplot)){
	sigmax2 = sigmax2_toplot[i]
	sigmay2 =sigmay2_toplot[i]
	r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	lines(r[11,1,],r[11,2,])	
	k=0
	arrows(r[11,1,steps-k],r[11,2,steps-k],r[11,1,steps+1],r[11,2,steps+1],length=0.05)
}
points(0,0,pch=20,cex=1,col=myred)
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

#### mode 12 
pdf(file='/Users/eleanorbrush/Desktop/mode12_dynamics.pdf',width=width,height=height,family=fontfamily)

rho = 1
sigmax2_lim = c(0,5)
sigmay2_lim = c(0,5)
steps = 1

par(ps=smallfontsize,mai=marg)
myvectorfield(mode12,sigmax2_lim,sigmay2_lim,n=12,C=50,scale=0.005,length=0.05,add_contour=FALSE)
sigmax2 = 4
sigmay2 = 2.05
steps = 20 
r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
lines(r[12,1,],r[12,2,])
k = 3
arrows(r[9,1,steps-k],r[9,2,steps-k],r[9,1,steps],r[9,2,steps],length=0.05)
points(0,0,pch=20,cex=1,col=myred)
mtext(expression(paste('Song variance, ',sigma[x]^2)),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(paste('Preference variance, ',sigma[y]^2)),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()

# #######  equilibrium and transient variance


recursion_all_end <- function(sigmax2,sigmay2,sigma2,equil=TRUE){
	r = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	toreturn = r[,1:2,steps+1]
	if(equil){		
			if(sigma2>sigmay2){
				toreturn[3,1:2] = c(NA,sigmay2)
			}else{toreturn[3,1:2] = c(sigmay2-sigma2,sigmay2)
				}
			if(diff(r[2,1,steps+(0:1)])<0 || r[2,1,steps]<1e-10){
				toreturn[2,1:2] = c(NA,sigmay2)
			}
		} else{
			if(sigma2<sigmay2){
				toreturn[3,1:2] = c(NA,sigmay2)
			} 
			if(diff(r[5,1,steps+(0:1)])>0){
				toreturn[5,1:2] = NA
			} 
			if(diff(r[2,1,steps+(0:1)])>0){
				toreturn[2,1:2] = c(NA,sigmay2)
			}
			}
	return(toreturn)
}

time_to_zero <- function(sigmax2,sigmay2,sigma2){
	r  = recursion_all_new_numbers(sigmax2,sigmay2,sigma2,rho)
	time = array(NA,c(12,2))
	for(i in 1:12){
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
steps_short = 25
steps_time = 5000
disappear_thresh <- 0.05 #5e-3

sigmax2 = 1
sigmay2 = 4.1
sigma2 = 1

vec_length = 100

pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system
lwd=2

#### sigmax2 as function of sigmax2

sigmax2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmax2 <- as.list(1:3)

sigmax2_pos = c(2,3)
sigmax2_zero = c(2,3,4,5,6,8,9,11,12)

sigmay2_pos = c(7,10)
sigmay2_zero = c(4,5,6,8,9,11,12)

steps = steps_long

r_sigmax2[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigma2=sigma2,sigmay2=sigmay2,simplify=FALSE)
r_sigmax2[[1]] <- matrix(unlist(lapply(r_sigmax2[[1]],function(y) y[,1])),nrow=12)
r_sigmax2[[1]][setdiff(1:12,sigmax2_pos),] = NA

steps = steps_short

r_sigmax2[[2]] <- sapply(sigmax2_vals,recursion_all_end,sigma2=sigma2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
r_sigmax2[[2]] <- matrix(unlist(lapply(r_sigmax2[[2]],function(y) y[,1])),nrow=12)
r_sigmax2[[2]][setdiff(1:12,sigmax2_zero),] = NA

steps = steps_time 

r_sigmax2[[3]] <- sapply(sigmax2_vals,time_to_zero,sigma2=sigma2,sigmay2=sigmay2,simplify=FALSE)
r_sigmax2[[3]] <- matrix(unlist(lapply(r_sigmax2[[3]],function(y) y[,1])),nrow=12)
r_sigmax2[[3]][setdiff(1:12,sigmax2_zero),] = NA


#####

marg = c(0.4,0.5,0.1,0.05)
omarg = c(0.5,0.7,0.3,0.0)

width = 6.8
height = 6

pdf('/Users/eleanorbrush/Desktop/sigmax2_sigmax2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:8,nrow=4,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(2,3,5,6,8,9,11,12),nrow=4,byrow=TRUE)

for(k in 1:4){
	subset = modes[k,]
	if(k==1){
		ylim = range(r_sigmax2[1:2],na.rm=TRUE)
		}else{
		ylim = c(0,0.1)	
		}
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmax2[[1]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmax2[[2]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial var. of songs, ', sigma[x](0)^2)),side=1,at=mean(sigmax2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of songs, ',sigma[x]^2)),side=2,at=mean(ylim)-0.19*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	ylim = range(r_sigmax2[3],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmax2[[3]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(2^(0:5)),labels=2^(0:5))
	axis(2,at=log(32),labels=32)
	mtext(expression(paste('Initial var. of songs, ', sigma[x](0)^2)),side=1,at=mean(sigmax2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(0.6,2.6,legend=c(expression(paste(sigma[x]^2, "*",', genetic')),expression(paste(sigma[x]^2~(25),', genetic')),expression(paste(sigma[x]^2, "*",', paternally learned')),expression(paste(sigma[x]^2~(25),', paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}
		
}


dev.off()

#### sigmax2 as function of sigmay2

sigmay2_vals = matrix(seq(sigmay2_min,sigmay2_max,length.out=vec_length),nrow=1)

r_sigmay2 <- as.list(1:3)

sigmax2_pos = c(2,3)
sigmax2_zero = c(2,3,4,5,6,8,9,11,12)

sigmay2_pos = c(7,10)
sigmay2_zero = c(4,5,6,8,9,11,12)

steps = steps_long

r_sigmay2[[1]] <- sapply(sigmay2_vals,recursion_all_end,sigma2=sigma2,sigmax2=sigmax2,simplify=FALSE)
r_sigmay2[[1]] <- matrix(unlist(lapply(r_sigmay2[[1]],function(y) y[,1])),nrow=12)
r_sigmay2[[1]][setdiff(1:12,sigmax2_pos),] = NA

steps = steps_short

r_sigmay2[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigma2=sigma2,sigmax2=sigmax2,equil=FALSE,simplify=FALSE)
r_sigmay2[[2]] <- matrix(unlist(lapply(r_sigmay2[[2]],function(y) y[,1])),nrow=12)
r_sigmay2[[2]][setdiff(1:12,sigmax2_zero),] = NA

steps = steps_time 

r_sigmay2[[3]] <- sapply(sigmay2_vals,time_to_zero,sigma2=sigma2,sigmax2=sigmax2,simplify=FALSE)
r_sigmay2[[3]] <- matrix(unlist(lapply(r_sigmay2[[3]],function(y) y[,1])),nrow=12)
r_sigmay2[[3]][setdiff(1:12,sigmax2_zero),] = NA

#####
lwd=2

marg = c(0.4,0.5,0.1,0.05)
omarg = c(0.5,0.7,0.3,0.0)

width = 6.8
height =6

pdf('/Users/eleanorbrush/Desktop/sigmax2_sigmay2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:8,nrow=4,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(2,3,5,6,8,9,11,12),nrow=4,byrow=TRUE)

for(k in 1:4){
	subset = modes[k,]
	if(k==1){
		ylim = c(0,4)
		}else{ylim = c(0,0.1)}
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmay2[[1]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmay2[[2]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial var. of preferences, ', sigma[y](0)^2)),side=1,at=mean(sigmay2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of songs, ',sigma[x]^2)),side=2,at=mean(ylim)-0.19*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(-0.3,4.5,legend=c(expression(paste(sigma[x]^2, "*",', song genetic')),expression(paste(sigma[x]^2~(25),', song genetic')),expression(paste(sigma[x]^2, "*",', song paternally learned')),expression(paste(sigma[x]^2~(25),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}
		
	ylim = range(r_sigmay2[3],na.rm=TRUE)
	
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigmay2[[3]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(2^(0:9)),labels=2^(0:9))
	mtext(expression(paste('Initial var. of preferences, ', sigma[y](0)^2)),side=1,at=mean(sigmay2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
}


dev.off() 

#### sigmay2 as a function of sigma2
steps_long = 1000
steps_short = 10
steps_time = 500

sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=vec_length),nrow=1)

r_sigma2 <- as.list(1:3)

sigmax2_pos = c(2,3)
sigmax2_zero = c(2,3,5,6,8,9,11,12)

sigmay2_pos = c(7,10)
sigmay2_zero = c(4,5,6,8,9,11,12)

steps = steps_long

r_sigma2[[1]] <- sapply(sigma2_vals,recursion_all_end,sigmay2=sigmay2,sigmax2=sigmax2,simplify=FALSE)
r_sigma2[[1]] <- matrix(unlist(lapply(r_sigma2[[1]],function(y) y[,2])),nrow=12)
r_sigma2[[1]][setdiff(1:12,sigmay2_pos),] = NA

steps = steps_short

r_sigma2[[2]] <- sapply(sigma2_vals,recursion_all_end,sigmay2=sigmay2,sigmax2=sigmax2,equil=FALSE,simplify=FALSE)
r_sigma2[[2]] <- matrix(unlist(lapply(r_sigma2[[2]],function(y) y[,2])),nrow=12)
r_sigma2[[2]][setdiff(1:12,sigmay2_zero),] = NA

steps = steps_time 

r_sigma2[[3]] <- sapply(sigma2_vals,time_to_zero,sigmay2=sigmay2,sigmax2=sigmax2,simplify=FALSE)
r_sigma2[[3]] <- matrix(unlist(lapply(r_sigma2[[3]],function(y) y[,2])),nrow=12)
r_sigma2[[3]][setdiff(1:12,sigmay2_zero),] = NA

#####
pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'Greens')[4],brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Greens')[7],brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.8
height =6

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigma2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(4:12,byrow=T,nrow=3)

for(k in 1:3){
	subset = modes[k,]
	if(k==1){
		ylim = c(0,3)
	}
	else{
	ylim = range(r_sigma2[1:2],na.rm=TRUE)
	}
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigma2[[1]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigma2[[2]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Var. of preference function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of preferences, ',sigma[y]^2)),side=2,at=mean(ylim)-0.18*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	
	ylim = range(r_sigma2[3],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigma2[[3]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*2^(0:9)),labels=5*2^(0:9))
	axis(2,at=log(320),labels=320)
	mtext(expression(paste('Var. of preference function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(0.2,6.2,legend=c(expression(paste(sigma[y]^2,"*",', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song genetic')),expression(paste(sigma[y]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],dark_col[2],dark_col[3]),bty='n')	
		}
}


dev.off()

#### sigmay2 as a function of sigmax2

sigmax2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmax2 <- as.list(1:3)

steps = steps_long

r_sigmax2[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
r_sigmax2[[1]] <- matrix(unlist(lapply(r_sigmax2[[1]],function(y) y[,2])),nrow=12)
r_sigmax2[[1]][setdiff(1:12,sigmay2_pos),] = NA

steps = steps_short

r_sigmax2[[2]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
r_sigmax2[[2]] <- matrix(unlist(lapply(r_sigmax2[[2]],function(y) y[,2])),nrow=12)
r_sigmax2[[2]][setdiff(1:12,sigmay2_zero),] = NA

steps = steps_time 

r_sigmax2[[3]] <- sapply(sigmax2_vals,time_to_zero,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
r_sigmax2[[3]] <- matrix(unlist(lapply(r_sigmax2[[3]],function(y) y[,2])),nrow=12)
r_sigmax2[[3]][setdiff(1:12,sigmay2_zero),] = NA

#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.8
height =6

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigmax2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

for(k in 1:3){
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
	mtext(expression(paste('Initial var. of songs, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=2.7,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of preferences, ',sigma[y]^2)),side=2,at=mean(ylim)-0.18*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	
	ylim = range(r_sigmax2[3],na.rm=TRUE)
	
	plot(sigmax2_vals,sigmax2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmax2[[3]][subset[j],]))!=0)){
			lines(sigmax2_vals,r_sigmax2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(2^(0:9)),labels=2^(0:9))
	mtext(expression(paste('Initial var. of songs, ',sigma[x]^2~(0))),side=1,at=mean(sigmax2_vals),line=2.7,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
		
	if(k==1){
		legend(-0.1,2.05,legend=c(expression(paste(sigma[y]^2,"*",', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song obliquely learned')),expression(paste(sigma[y]^2~(10),', song genetic')),expression(paste(sigma[y]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],dark_col[2],dark_col[3]),bty='n')	
		}
	
}


dev.off()


#### sigmay2 as a function of sigmay2

sigmay2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)

r_sigmay2 <- as.list(1:3)

steps = steps_long

r_sigmay2[[1]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
r_sigmay2[[1]] <- matrix(unlist(lapply(r_sigmay2[[1]],function(y) y[,2])),nrow=12)
r_sigmay2[[1]][setdiff(1:12,sigmay2_pos),] = NA

steps = steps_short

r_sigmay2[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
r_sigmay2[[2]] <- matrix(unlist(lapply(r_sigmay2[[2]],function(y) y[,2])),nrow=12)
r_sigmay2[[2]][setdiff(1:12,sigmay2_zero),] = NA

steps = steps_time 

r_sigmay2[[3]] <- sapply(sigmay2_vals,time_to_zero,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
r_sigmay2[[3]] <- matrix(unlist(lapply(r_sigmay2[[3]],function(y) y[,2])),nrow=12)
r_sigmay2[[3]][setdiff(1:12,sigmay2_zero),] = NA


#####
lwd=2

marg = c(0.53,0.5,0.04,0.05)
omarg = c(0.4,0.7,0.3,0.0)

width = 6.8
height =6

pdf('/Users/eleanorbrush/Desktop/sigmay2_sigmay2.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

for(k in 1:3){
	subset = modes[k,]
	if(k>1){
	ylim = range(r_sigmay2[1:2],na.rm=TRUE)
	}else{ylim=c(0,0.012)}
	plot(sigmay2_vals,sigmay2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:3){
		if(length(which(!is.na(r_sigmay2[[1]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigmay2[[2]][subset[j],]))!=0)){
			lines(sigmay2_vals,r_sigmay2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Initial var. of preferences, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=2.7,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of preferences, ',sigma[y]^2)),side=2,at=mean(ylim)-0.18*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
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
	axis(2,at=log(2^(0:9)),labels=2^(0:9))
	mtext(expression(paste('Initial var. of preferences, ',sigma[y]^2~(0))),side=1,at=mean(sigmay2_vals),line=2.7,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
		
}

dev.off()

# source('pearson_kurtosis.R)
# sourece('mutation_test.R')
# source('step_function_example.R')

##### effect of step trait distributions on equilibrium distributions
trait_chunk_num = 281
source('range_setup.R')

col_vec = c(brewer.pal(9,'Set1')[c(1,7)],brewer.pal(9,'Set1')[c(9,5,8)],brewer.pal(9,'Blues')[9],brewer.pal(9,'Set1')[2:4],brewer.pal(8,'Dark2')[1])
lwd = 2
marg = c(0.4,0.4,0.03,0.05)
omarg = c(0.3,0.7,1.3,0.0)

width = 6.8
height = 3.5

pdf('/Users/eleanorbrush/Desktop/effect_of_step_traits_on_peaks.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=TRUE))

load('step_song_equilibrium.Rdata')

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2_vals[1],NA,NA)
m_init = ic$m_init

xlim = c(-6,6)
ylim=c(0,.14)
plot(mrange+1,m_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song, x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('Initial conditions',side=3,line=0,cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		ic = init_conds('step','norm','norm',sigmax2,sigmay2,sigma2_vals[1],k1,k2)
		m_init = ic$m_init		
		
		if(sum(is.na(m_init))==0){			
			lines(mrange+1,m_init,col=col_vec[j+1],lwd=lwd)
		}
}

xlim = c(-6,6)
ylim=c(0,.14)
k=5
plot(mrange+1,equilibrium[[1,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song, x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('Mechanism 3',side=3,line=0,cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[1,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[1,k,j+1]],col=col_vec[j+1],lwd=lwd)
		}
}

xlim = c(-6,6)
ylim=c(0,.14)
k=5
plot(mrange+1,equilibrium[[2,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('Mechanism 7',side=3,line=0,cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[2,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[2,k,j+1]],col=col_vec[j+1],lwd=lwd)
		}
}

load('step_pref_equilibrium.Rdata')

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2_vals[1],NA,NA)
f_init = ic$f_init


xlim = c(-6,6)
ylim=c(0,.14)
plot(mrange+1,f_init,col=col_vec[6],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		ic = init_conds('norm','step','norm',sigmax2,sigmay2,sigma2_vals[1],k1,k2)
		f_init = ic$f_init	
		
		if(sum(is.na(f_init))==0){
			lines(mrange+1,f_init,col=col_vec[j+6],lwd=lwd)
		}
}

xlim = c(-6,6)
ylim=c(0,.9)
k=5
plot(mrange+1,equilibrium[[1,k,1]],col=col_vec[6],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song, x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[1,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[1,k,j+1]],col=col_vec[j+6],lwd=lwd)
		}
}

xlim = c(-6,6)
ylim=c(0,.14)
k=5
plot(mrange+1,equilibrium[[2,k,1]],col=col_vec[6],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
# mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[2,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[2,k,j+1]],col=col_vec[j+6],lwd=lwd)
		}
}


dev.off()

### effect of step distributions on equilibrium variance
trait_chunk_num = 281
source('range_setup.R')

col_vec = c(brewer.pal(9,'Set1')[c(1,7)],brewer.pal(9,'Set1')[c(9,5,8)],brewer.pal(9,'Blues')[9],brewer.pal(9,'Set1')[2:4],brewer.pal(8,'Dark2')[1])
lwd = 2
marg = c(0.4,0.4,0.03,0.05)
omarg = c(0.3,0.7,1.3,0.0)

width = 6.8
height = 3.5

pdf('/Users/eleanorbrush/Desktop/effect_of_step_traits_on_variance.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=TRUE))

load('step_song_equilibrium.Rdata')

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2_vals[1],NA,NA)
m_init = ic$m_init

xlim = c(-6,6)
ylim=c(0,.14)
plot(mrange+1,m_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song, x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext('Initial conditions',side=3,line=0,cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		ic = init_conds('step','norm','norm',sigmax2,sigmay2,sigma2_vals[1],k1,k2)
		m_init = ic$m_init		
		
		if(sum(is.na(m_init))==0){			
			lines(mrange+1,m_init,col=col_vec[j+1],lwd=lwd)
		}
}


for(i in 1:2){
	plot(sigma2_vals,var_mat[i,,1],xlim=round(range(sigma2_vals)),ylim=c(0,2),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='',xaxt='n')
	axis(1,at=seq(0,2,by=0.5))
	mtext(expression(paste('Var. of pref. function, ', sigma^2)),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste('Eq. song var., ',sigma[x]^2, "*")),expression(paste('Eq. pref. var., ',sigma[y]^2, "*")))[[i]],side=2,line=1.3,at=1,cex=largefontsize/smallfontsize)
	mtext(list('Mechanism 3','Mechanism 7')[[i]],side=3,line=0,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+1],t='o',lwd=lwd)}
	}
}

load('step_pref_equilibrium.Rdata')

ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2_vals[1],NA,NA)
f_init = ic$f_init	


xlim = c(-6,6)
ylim=c(0,.14)
plot(mrange+1,f_init,col=col_vec[6],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference, y',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){		
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		ic = init_conds('norm','step','norm',sigmax2,sigmay2,sigma2_vals[1],k1,k2)
		f_init = ic$f_init	
		
		if(sum(is.na(f_init))==0){
			lines(mrange+1,f_init,col=col_vec[j+6],lwd=lwd)
		}
}


for(i in 1:2){
	plot(sigma2_vals,var_mat[i,,1],xlim=round(range(sigma2_vals)),ylim=c(0,2),t='o',lwd=lwd,col=col_vec[6],xlab='',ylab='',xaxt='n')
	axis(1,at=seq(0,2,by=0.5))
	mtext(expression(paste('Var. of pref. function, ', sigma^2)),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste('Eq. song var., ',sigma[x]^2, "*")),expression(paste('Eq. pref. var., ',sigma[y]^2, "*")))[[i]],side=2,line=1.3,at=1,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+6],t='o',lwd=lwd)}
	}
}

dev.off()

##### effect of step pref function on peaks and var
trait_chunk_num = 281
source('range_setup.R')

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.4,0.4,0.03,0.05)
omarg = c(0.3,0.7,1.3,0.0)

width = 6.8
height = 3.5

load('step_pref_fun_equilibrium.Rdata')

pdf('/Users/eleanorbrush/Desktop/effect_of_step_pref_fun.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=FALSE))

for(i in 1:2){
	sigma2 = sigma2_vals[c(3,9)[i]]
	
	ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,NA,NA)
	fixed_weight = ic$fixed_weight	
	
	xlim = c(-6,6)
	ylim=c(0,c(.15,.15)[i])
	plot(mrange+1,fixed_weight,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
	mtext('Preference - song, y-x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
	mtext('Preference',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
	if(i==1){mtext('Preference function',side=3,line=0,cex=largefontsize/smallfontsize)}
	for(j in rev(1:length(k1_vals))){		
			k1 = k1_vals[j]
			k2 = k2_vals[j]
			
			ic = init_conds('norm','norm','step',sigmax2,sigmay2,sigma2,k1,k2)
			fixed_weight = ic$fixed_weight			
			
			if(sum(is.na(fixed_weight))==0){
				lines(mrange+1,fixed_weight,col=col_vec[j+1],lwd=lwd)
			}
	}
}

k=9

for(i in 1:2){
xlim = c(-6,6)
ylim=c(0,0.08)
plot(mrange+1,equilibrium[[i,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext(list('Song, x','Preference, y')[[i]],side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
mtext(list('Mechanism 3','Mechanism 7')[[i]],side=3,line=0,cex=largefontsize/smallfontsize)
for(j in rev(1:length(k1_vals))){					
		if(!is.na(equilibrium[[i,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[i,k,j+1]],col=col_vec[j+1],lwd=lwd)
		}
}

plot(sigma2_vals,var_mat[i,,1],xlim=round(range(sigma2_vals)),ylim=c(0,round(max(var_mat[1:2,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='',xaxt='n')
axis(1,at=seq(0,2,by=0.5))
mtext(expression(paste('Var. of pref. function, ', sigma^2)),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
mtext(list(expression(paste('Eq. song var., ',sigma[x]^2, "*")),expression(paste('Eq. pref. var., ',sigma[y]^2, "*")))[[i]],side=2,line=1.3,cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){
	if(length(which(!is.na(var_mat[i,,j+1])))!=0){
		points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+1],t='o',lwd=lwd)}
}
}


dev.off()



