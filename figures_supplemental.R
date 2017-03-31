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
# pdf(file='/Users/eleanorbrush/Desktop/mode7_dynamics.pdf',width=width,height=height,family=fontfamily)

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

# dev.off()

# # # # # mode 9 
# pdf(file='/Users/eleanorbrush/Desktop/mode9_dynamics.pdf',width=width,height=height,family=fontfamily)

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

# dev.off()

#### mode 2 
# pdf(file='/Users/eleanorbrush/Desktop/mode2_dynamics.pdf',width=width,height=height,family=fontfamily)

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

# dev.off()

#### mode 8
# pdf(file='/Users/eleanorbrush/Desktop/mode8_dynamics.pdf',width=width,height=height,family=fontfamily)

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

# dev.off()

#### mode 1
# pdf(file='/Users/eleanorbrush/Desktop/mode1_dynamics.pdf',width=width,height=height,family=fontfamily)

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

# dev.off()

# ####### effect of sigma2 on equilibrium and transient variance

# sigmax2_max = 2
# sigmax2_min = 0.01
# sigmay2_max = 4.5
# sigmay2_min = 0.01
# sigma2_max = 4.5
# sigma2_min = 0.01
# rho = 0.6
# steps_long = 5000
# steps_short = 10
# steps_time = 1000
# disappear_thresh <- 5e-3

# sigmax2 = 0.5
# sigmay2 = 4
# sigma2 = 1

# vec_length = 200

# sigmax2_vals = matrix(seq(sigmax2_min,sigmax2_max,length.out=vec_length),nrow=1)
# sigmay2_vals = matrix(seq(sigmay2_min,sigmay2_max,length.out=vec_length),nrow=1)
# sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=vec_length),nrow=1)
# variables=rbind(sigmax2_vals,sigmay2_vals,sigma2_vals)

# sigmax2_toplot <- as.list(1:9)
# dim(sigmax2_toplot) = c(3,3)

# sigmay2_toplot <- as.list(1:9)
# dim(sigmay2_toplot) = c(3,3)

# sigmax2_pos = c(9,3)
# sigmax2_zero = c(7,9,1,2,8,3)

# sigmay2_pos = c(5)
# sigmay2_zero = c(7,1,2,4,8)

# steps = steps_long

# sigmax2_toplot[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
# sigmax2_toplot[[1]] <- matrix(unlist(lapply(sigmax2_toplot[[1]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[1]][setdiff(1:9,sigmax2_pos),] = NA

# sigmax2_toplot[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
# sigmax2_toplot[[2]] <- matrix(unlist(lapply(sigmax2_toplot[[2]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[2]][setdiff(1:9,sigmax2_pos),] = NA

# sigmax2_toplot[[3]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
# sigmax2_toplot[[3]] <- matrix(unlist(lapply(sigmax2_toplot[[3]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[3]][setdiff(1:9,sigmax2_pos),] = NA

# steps = steps_short

# sigmax2_toplot[[4]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
# sigmax2_toplot[[4]] <- matrix(unlist(lapply(sigmax2_toplot[[4]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[4]][setdiff(1:9,sigmax2_zero),] = NA

# sigmax2_toplot[[5]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
# sigmax2_toplot[[5]] <- matrix(unlist(lapply(sigmax2_toplot[[5]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[5]][setdiff(1:9,sigmax2_zero),] = NA

# sigmax2_toplot[[6]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
# sigmax2_toplot[[6]] <- matrix(unlist(lapply(sigmax2_toplot[[6]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[6]][setdiff(1:9,sigmax2_zero),] = NA

# steps = steps_time 

# sigmax2_toplot[[7]] <- sapply(sigmax2_vals,time_to_zero,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
# sigmax2_toplot[[7]] <- matrix(unlist(lapply(sigmax2_toplot[[7]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[7]][setdiff(1:9,sigmax2_zero),] = NA


# sigmax2_toplot[[8]] <- sapply(sigmay2_vals,time_to_zero,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
# sigmax2_toplot[[8]] <- matrix(unlist(lapply(sigmax2_toplot[[8]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[8]][setdiff(1:9,sigmax2_zero),] = NA


# sigmax2_toplot[[9]] <- sapply(sigma2_vals,time_to_zero,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
# sigmax2_toplot[[9]] <- matrix(unlist(lapply(sigmax2_toplot[[9]],function(y) y[,1])),nrow=9)
# sigmax2_toplot[[9]][setdiff(1:9,sigmax2_zero),] = NA

# steps = steps_long

# sigmay2_toplot[[1]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
# sigmay2_toplot[[1]] <- matrix(unlist(lapply(sigmay2_toplot[[1]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[1]][setdiff(1:9,sigmay2_pos),] = NA

# sigmay2_toplot[[2]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
# sigmay2_toplot[[2]] <- matrix(unlist(lapply(sigmay2_toplot[[2]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[2]][setdiff(1:9,sigmay2_pos),] = NA

# sigmay2_toplot[[3]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
# sigmay2_toplot[[3]] <- matrix(unlist(lapply(sigmay2_toplot[[3]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[3]][setdiff(1:9,sigmay2_pos),] = NA

# steps = steps_short

# sigmay2_toplot[[4]] <- sapply(sigmax2_vals,recursion_all_end,sigmay2=sigmay2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
# sigmay2_toplot[[4]] <- matrix(unlist(lapply(sigmay2_toplot[[4]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[4]][setdiff(1:9,sigmay2_zero),] = NA

# sigmay2_toplot[[5]] <- sapply(sigmay2_vals,recursion_all_end,sigmax2=sigmax2,sigma2=sigma2,equil=FALSE,simplify=FALSE)
# sigmay2_toplot[[5]] <- matrix(unlist(lapply(sigmay2_toplot[[5]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[5]][setdiff(1:9,sigmay2_zero),] = NA

# sigmay2_toplot[[6]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
# sigmay2_toplot[[6]] <- matrix(unlist(lapply(sigmay2_toplot[[6]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[6]][setdiff(1:9,sigmay2_zero),] = NA

# steps = steps_time 

# sigmay2_toplot[[7]] <- sapply(sigmax2_vals,time_to_zero,sigmay2=sigmay2,sigma2=sigma2,simplify=FALSE)
# sigmay2_toplot[[7]] <- matrix(unlist(lapply(sigmay2_toplot[[7]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[7]][setdiff(1:9,sigmay2_zero),] = NA

# sigmay2_toplot[[8]] <- sapply(sigmay2_vals,time_to_zero,sigmax2=sigmax2,sigma2=sigma2,simplify=FALSE)
# sigmay2_toplot[[8]] <- matrix(unlist(lapply(sigmay2_toplot[[8]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[8]][setdiff(1:9,sigmay2_zero),] = NA


# sigmay2_toplot[[9]] <- sapply(sigma2_vals,time_to_zero,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
# sigmay2_toplot[[9]] <- matrix(unlist(lapply(sigmay2_toplot[[9]],function(y) y[,2])),nrow=9)
# sigmay2_toplot[[9]][setdiff(1:9,sigmay2_zero),] = NA

# ltys_vec = c(2,4,6,5,6,NA,1,3,5)
# ltys_vec = c(1,1,1,1,1,NA,1,1,1)
# pal1 = brewer.pal(9,'Set1')
# pal2 = brewer.pal(8,'Accent')
# col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
# old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

# #####
# lwd=2

# marg = c(0.53,0.5,0.04,0.05)
# omarg = c(0.4,0.7,0.3,0.0)

# width = 6.5
# height =3

# par(mfrow=c(3,3),ps=smallfontsize,mai=marg,oma=omarg)

# xlab_list = as.list(rep(c(expression(sigma[x]^2~(0)),expression(sigma[y]^2~(0)),expression(sigma^2)),times=3))
# ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

# fem = matrix(c(1:3,7:9),nrow=3)

# for(k in c(3,1,2)){
	# subset = fem[k,]
# for(i in c(3,6,9)){
	# ylim = range(sigmax2_toplot[rep(3,3)+rep(c(0,3,6),each=3)[i]],na.rm=TRUE)
	# if(i==6){ylim=c(0,1)}
	# plot(variables[rep(1:3,times=3)[i],],variables[rep(1:3,times=3)[i],],ylim=ylim,xlab='',ylab='',t='n',yaxt=c(rep('t',6),rep('n',3))[i])
	# if(is.element(i,7:9)){
		# axis(2,at=log(5*2^(0:9)),labels=5*2^(0:9))
	# }
	# for(j in subset){
		# if(length(which(!is.na(sigmax2_toplot[[i]][j,])))!=0){
			# lines(variables[rep(1:3,times=3)[i],],sigmax2_toplot[[i]][j,],lty=ltys_vec[j],col=col_vec[j],lwd=lwd)
		# }
	# }
	# mtext(xlab_list[[i]],side=1,at=mean(variables[rep(1:3,times=3)[i],]),line=3,cex=largefontsize/smallfontsize)
	# mtext(ylab_list[[i-2]],side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	# if(i==3){
		# legend(2.5,3.7,c(2,3,5,6,8,9),lty=ltys_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],col=col_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],bty='n',lwd=lwd)
	# }
# }

# }

# ####
# lwd=2

# marg = c(0.53,0.5,0.04,0.05)
# omarg = c(0.4,0.7,0.3,0.0)

# width = 6.5
# height =3

# pdf(file='/Users/eleanorbrush/Desktop/sigmax2_sigma2.pdf',width=width,height=height,family=fontfamily)

# par(mfrow=c(1,3),ps=smallfontsize,mai=marg,oma=omarg)

# xlab_list = as.list(rep(c(expression(sigma[x]^2~(0)),expression(sigma[y]^2~(0)),expression(sigma^2)),times=3))
# ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

# for(i in c(3,6,9)){
	# ylim = range(sigmax2_toplot[rep(3,3)+rep(c(0,3,6),each=3)[i]],na.rm=TRUE)
	# if(i==6){ylim=c(0,1)}
	# plot(variables[rep(1:3,times=3)[i],],variables[rep(1:3,times=3)[i],],ylim=ylim,xlab='',ylab='',t='n',yaxt=c(rep('t',6),rep('n',3))[i])
	# if(is.element(i,7:9)){
		# axis(2,at=log(5*2^(0:9)),labels=5*2^(0:9))
	# }
	# for(j in 1:9){
		# if(length(which(!is.na(sigmax2_toplot[[i]][j,])))!=0){
			# lines(variables[rep(1:3,times=3)[i],],sigmax2_toplot[[i]][j,],lty=ltys_vec[j],col=col_vec[j],lwd=lwd)
		# }
	# }
	# mtext(xlab_list[[i]],side=1,at=mean(variables[rep(1:3,times=3)[i],]),line=3,cex=largefontsize/smallfontsize)
	# mtext(ylab_list[[i-2]],side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	# if(i==3){
		# legend(2.5,3.7,c(2,3,5,6,8,9),lty=ltys_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],col=col_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],bty='n',lwd=lwd)
	# }
# }

# dev.off()

# marg = c(0.53,0.3,0.04,0.05)
# omarg = c(0.4,1.5,0.3,0.0)

# width = 6.5
# height =6

# pdf(file='/Users/eleanorbrush/Desktop/sigmax2_full.pdf',width=width,height=height,family=fontfamily)

# par(mfrow=c(3,3),ps=smallfontsize,mai=marg,oma=omarg)

# xlab_list = as.list(rep(c(expression(sigma[x]^2~(0)),expression(sigma[y]^2~(0)),expression(sigma^2)),times=3))
# ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

# for(i in 1:9){
	# ylim = range(sigmax2_toplot[(1:3)+rep(c(0,3,6),each=3)[i]],na.rm=TRUE)
	# plot(variables[rep(1:3,times=3)[i],],variables[rep(1:3,times=3)[i],],ylim=ylim,xlab='',ylab='',t='n',yaxt=c(rep('t',6),rep('n',3))[i])
	# if(is.element(i,7:9)){
		# axis(2,at=log(5*2^(0:9)),labels=5*2^(0:9))
	# }
	# for(j in 1:9){
		# if(length(which(!is.na(sigmax2_toplot[[i]][j,])))!=0){
			# lines(variables[rep(1:3,times=3)[i],],sigmax2_toplot[[i]][j,],lty=ltys_vec[j],col=col_vec[j],lwd=lwd)
		# }
	# }
	# mtext(xlab_list[[i]],side=1,at=mean(variables[rep(1:3,times=3)[i],]),line=3,cex=largefontsize/smallfontsize)
	# mtext(ylab_list[[i]],side=2,at=mean(ylim),line=1.8,cex=largefontsize/smallfontsize)
	# if(i==2){
		# legend(0,3.5,c(2,3,5,6,8,9),lty=ltys_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],col=col_vec[unique(c(sigmax2_pos,sigmax2_zero))][order(old_to_new[unique(c(sigmax2_pos,sigmax2_zero))])],bty='n',lwd=lwd)
	# }
# }


# dev.off()

# pdf(file='/Users/eleanorbrush/Desktop/sigmay2_full.pdf',width=width,height=height,family=fontfamily)

# par(mfrow=c(3,3),ps=smallfontsize,mai=marg,oma=omarg)

# xlab_list = as.list(rep(c(expression(sigma[x]^2~(0)),expression(sigma[y]^2~(0)),expression(sigma^2)),times=3))
# ylab_list = as.list(c(rep(c(expression(sigma[y]^2),'',''),times=2),c('Time','','')))

# for(i in 1:9){
	# ylim = range(sigmay2_toplot[(1:3)+rep(c(0,3,6),each=3)[i]],na.rm=TRUE)
	# plot(variables[rep(1:3,times=3)[i],],variables[rep(1:3,times=3)[i],],ylim=ylim,xlab='',ylab='',t='n',yaxt=c(rep('t',6),rep('n',3))[i])
	# if(is.element(i,7:9)){
		# axis(2,at=log(5*2^(0:9)),labels=5*2^(0:9))
	# }
	# for(j in 1:9){
		# if(length(which(!is.na(sigmay2_toplot[[i]][j,])))!=0){
			# lines(variables[rep(1:3,times=3)[i],],sigmay2_toplot[[i]][j,],lty=ltys_vec[j],col=col_vec[j],lwd=lwd)
		# }
	# }
	# mtext(xlab_list[[i]],side=1,at=mean(variables[rep(1:3,times=3)[i],]),line=3,cex=largefontsize/smallfontsize)
	# mtext(ylab_list[[i]],side=2,at=mean(ylim),line=1.7,cex=largefontsize/smallfontsize)
	# if(i==1){
		# legend(0,1.3,4:9,lty=ltys_vec[unique(c(sigmay2_pos,sigmay2_zero))][order(old_to_new[unique(c(sigmay2_pos,sigmay2_zero))])],col=col_vec[unique(c(sigmay2_pos,sigmay2_zero))][order(old_to_new[unique(c(sigmay2_pos,sigmay2_zero))])],bty='n',lwd=lwd)

	# }
# }

# dev.off()

# source('step_function_example.R')