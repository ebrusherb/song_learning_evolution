setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
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
mtext(expression(sigma[x]^2),side=1,at=mean(sigmax2_lim),line=2.5,cex=largefontsize/smallfontsize)
mtext(expression(sigma[y]^2),side=2,at=mean(sigmay2_lim),line=2,cex=largefontsize/smallfontsize)

dev.off()
