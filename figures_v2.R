setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
source('dynamics_by_mode.R')
source('init_conds.R')
library(RColorBrewer)
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12
lwd=2
pal1 = brewer.pal(9,'Set1')
pal2 = brewer.pal(8,'Accent')
col_vec = pal1[c(1,2,8,7,4,NA,3,5,9)]
pale_col = c(brewer.pal(9,'RdPu')[4],brewer.pal(9,'Blues')[4])
dark_col = c(brewer.pal(9,'Reds')[7],brewer.pal(9,'Blues')[7])
old_to_new = c(6,9,3,4,7,1,5,8,2) # converts old mode numbering system to new mode numbering system

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

####### effect of sigma2 on equilibrium and transient variance

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

sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=vec_length),nrow=1)

r_sigma2 <- as.list(1:3)

sigmax2_pos = c(2,3)
sigmax2_zero = c(2,3,4,5,6,8,9,11,12)

sigmay2_pos = c(7)
sigmay2_zero = c(4,5,6,8,9)

steps = steps_long

r_sigma2[[1]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
r_sigma2[[1]] <- matrix(unlist(lapply(r_sigma2[[1]],function(y) y[,1])),nrow=12)
r_sigma2[[1]][setdiff(1:12,sigmax2_pos),] = NA

steps = steps_short

r_sigma2[[2]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
r_sigma2[[2]] <- matrix(unlist(lapply(r_sigma2[[2]],function(y) y[,1])),nrow=12)
r_sigma2[[2]][setdiff(1:12,sigmax2_zero),] = NA

steps = steps_time 

r_sigma2[[3]] <- sapply(sigma2_vals,time_to_zero,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
r_sigma2[[3]] <- matrix(unlist(lapply(r_sigma2[[3]],function(y) y[,1])),nrow=12)
r_sigma2[[3]][setdiff(1:12,sigmax2_zero),] = NA

####

marg = c(0.4,0.5,0.1,0.05)
omarg = c(0.5,0.7,0.3,0.0)

width = 6.8
height = 6

pdf('/Users/eleanorbrush/Desktop/sigmax2_by_female_mode.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:8,nrow=4,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(2,3,5,6,8,9,11,12),nrow=4,byrow=TRUE)

for(k in 1:4){
	subset = modes[k,]
	if(is.element(k,c(1,3,4))){
		ylim = range(r_sigma2[1:2],na.rm=TRUE)
	}else if(k==2){
		ylim = c(min(unlist(lapply(r_sigma2[1:2],min,na.rm=TRUE))),0.4)
	}else{
		ylim = c(min(unlist(lapply(r_sigma2[1:2],min,na.rm=TRUE))),1.5)
	}		
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigma2[[1]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigma2[[2]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(paste('Var. of preference function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext(expression(paste('Var. of songs, ',sigma[x]^2)),side=2,at=mean(ylim)-0.19*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(2,4.5,legend=c(expression(paste(sigma[x]^2, "*",', genetic')),expression(paste(sigma[x]^2~(25),', genetic')),expression(paste(sigma[x]^2, "*",', paternally learned')),expression(paste(sigma[x]^2~(25),', paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}

	
	ylim = range(r_sigma2[3],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigma2[[3]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(c(6,25*4^(0:9))),labels=c(6,25*4^(0:9)))
	axis(2,at=log(1600),labels=1600)
	mtext(expression(paste('Var. of preference function, ', sigma^2)),side=1,at=mean(sigma2_vals),line=2.5,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim)-0.06*diff(ylim),line=2,cex=largefontsize/smallfontsize)
	
# # 	if(k==1){
		# legend(-0.25,6.4,legend=c(expression(paste(sigma[x]^2, "*",', song genetic')),expression(paste(sigma[x]^2~(20),', song genetic')),expression(paste(sigma[x]^2, "*",', song paternally learned')),expression(paste(sigma[x]^2~(20),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		# }
}


dev.off()

##### effect of step trait distributions

load('all_modes_numerically_step_song_dist.Rdata')
load('all_modes_numerically_step_song_dist2.Rdata')

D_song = list(d3,d5)
rm(d2,d4,d6,d7,d8,d9,d11,d12)

load('all_modes_numerically_step_pref_dist.Rdata')
load('all_modes_numerically_step_pref_dist2.Rdata')

D_pref = list(d3,d5)
rm(d2,d4,d6,d7,d8,d9,d11,d12)

source('range_setup_long.R')

marg = c(0.35,0.4,0.1,0.05)
omarg = c(0.1,0.7,0.1,0.0)

width = 6.8
height = 3

pdf('/Users/eleanorbrush/Desktop/effect_of_step_trait_distribution.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:6,ncol=3,byrow=TRUE))

song = 'step'
pref = 'norm'
func = 'norm'

rho = 0.9
ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

pmax = 0.2
xlim=c(-6,6)

plot(mrange+1,m_init,t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)
lines(mrange+1,f_init,lwd=lwd,col=col_vec[2])
# mtext('Song',side=1,line=2,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,cex=largefontsize/smallfontsize)
legend(-6,0.2,legend=c('Male song','Female preference'),bty='n',lwd=lwd,col=col_vec,lty=1)

for(i in 1:2){	
	t = min(dim(D_song[[i]]$Pm)[length(dim(D_song[[i]]$Pm))],1001)
	n = length(dim(D_song[[i]]$Pm))
	if(n==2){
		plot(mrange+1,D_song[[i]]$Pm[,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)}else{
		plot(mrange+1,apply(D_song[[i]]$Pm[,,t],1,sum),t='l',col=col_vec[1],lwd=lwd,xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)
	}
	n = length(dim(D_song[[i]]$Pf))
	if(n==2){
		lines(mrange+1,D_song[[i]]$Pf[,t],t='l',lwd=lwd,col=col_vec[2])}else{
		lines(mrange+1,apply(D_song[[i]]$Pf[,,t],2,sum),t='l',col=col_vec[2],lwd=lwd)
	}
	# mtext('Song',side=1,line=2,cex=largefontsize/smallfontsize)
	mtext('Frequency',side=2,line=2,cex=largefontsize/smallfontsize)
}

song = 'norm'
pref = 'step'
func = 'norm'

rho = 0.9
ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

plot(mrange+1,m_init,t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)
lines(mrange+1,f_init,lwd=lwd,col=col_vec[2])
# mtext('Song',side=1,line=2,cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,cex=largefontsize/smallfontsize)

for(i in 1:2){
	t = min(dim(D_pref[[i]]$Pm)[length(dim(D_pref[[i]]$Pm))],1001)
	n = length(dim(D_pref[[i]]$Pm))
	if(n==2){
		plot(mrange+1,D_pref[[i]]$Pm[,t],t='l',lwd=lwd,col=col_vec[1],xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)}else{
		plot(mrange+1,apply(D_pref[[i]]$Pm[,,t],1,sum),t='l',col=col_vec[1],lwd=lwd,xlab='',ylab='',ylim=c(0,pmax),xlim=xlim)
	}
	n = length(dim(D_pref[[i]]$Pf))
	if(n==2){
		lines(mrange+1,D_pref[[i]]$Pf[,t],t='l',lwd=lwd,col=col_vec[2])}else{
		lines(mrange+1,apply(D_pref[[i]]$Pf[,,t],2,sum),t='l',col=col_vec[2],lwd=lwd)
	}	
	# mtext('Song',side=1,line=2,cex=largefontsize/smallfontsize)
	mtext('Frequency',side=2,line=2,cex=largefontsize/smallfontsize)
}

dev.off()

# source('peak_example.R')