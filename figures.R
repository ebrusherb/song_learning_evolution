setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('multi.outer.R')
source('recursion_all.R')
source('dynamics_by_mode.R')
library(RColorBrewer)
library(pracma)
fontfamily = 'Helvetica'
smallfontsize = 10
largefontsize = 12

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

####### effect of sigma2 on equilibrium and transient variance

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

sigma2_vals = matrix(seq(sigma2_min,sigma2_max,length.out=vec_length),nrow=1)

r_sigma2 <- as.list(1:3)

sigmax2_pos = c(9,3)
sigmax2_zero = c(7,9,1,2,8,3)

sigmay2_pos = c(5)
sigmay2_zero = c(7,1,2,4,8)

steps = steps_long

r_sigma2[[1]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
r_sigma2[[1]] <- matrix(unlist(lapply(r_sigma2[[1]],function(y) y[,1])),nrow=9)
r_sigma2[[1]][setdiff(1:9,sigmax2_pos),] = NA

steps = steps_short

r_sigma2[[2]] <- sapply(sigma2_vals,recursion_all_end,sigmax2=sigmax2,sigmay2=sigmay2,equil=FALSE,simplify=FALSE)
r_sigma2[[2]] <- matrix(unlist(lapply(r_sigma2[[2]],function(y) y[,1])),nrow=9)
r_sigma2[[2]][setdiff(1:9,sigmax2_zero),] = NA

steps = steps_time 

r_sigma2[[3]] <- sapply(sigma2_vals,time_to_zero,sigmax2=sigmax2,sigmay2=sigmay2,simplify=FALSE)
r_sigma2[[3]] <- matrix(unlist(lapply(r_sigma2[[3]],function(y) y[,1])),nrow=9)
r_sigma2[[3]][setdiff(1:9,sigmax2_zero),] = NA

ltys_vec = c(2,4,6,5,6,NA,1,3,5)
ltys_vec = c(1,1,1,1,1,NA,1,1,1)
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

pdf('/Users/eleanorbrush/Desktop/sigmax2_by_female_mode.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg)
layout(matrix(1:6,nrow=3,byrow=TRUE))
ylab_list = as.list(c(rep(c(expression(sigma[x]^2),'',''),times=2),c('Time','','')))

modes = matrix(c(1:3,7:9)[order(old_to_new[c(1:3,7:9)])],nrow=3,byrow=TRUE)

for(k in 1:3){
	subset = modes[k,]
	
	ylim = range(r_sigma2[1:2],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigma2[[1]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[1]][subset[j],],lwd=lwd,col=pale_col[j])
		}
		if(length(which(!is.na(r_sigma2[[2]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[2]][subset[j],],lwd=lwd,col=dark_col[j])
		}
	}
	mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext(expression(sigma[x]^2),side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
	
	if(k==1){
		legend(1.4,4.4,legend=c(expression(paste(sigma[x]^2, "*",', song genetic')),expression(paste(sigma[x]^2~(10),', song genetic')),expression(paste(sigma[x]^2, "*",', song paternally learned')),expression(paste(sigma[x]^2~(10),', song paternally learned'))),lty=rep(1,4),col=c(pale_col[1],dark_col[1],pale_col[2],dark_col[2]),bty='n')	
		}
	
	ylim = range(r_sigma2[3],na.rm=TRUE)
	
	plot(sigma2_vals,sigma2_vals,ylim=ylim,xlab='',ylab='',t='n',yaxt='n')
	for(j in 1:2){
		if(length(which(!is.na(r_sigma2[[3]][subset[j],]))!=0)){
			lines(sigma2_vals,r_sigma2[[3]][subset[j],],lwd=lwd,col=dark_col[j])
		}			
	}
	axis(2,at=log(5*4^(0:9)),labels=5*4^(0:9))
	mtext(expression(sigma^2),side=1,at=mean(sigma2_vals),line=3,cex=largefontsize/smallfontsize)
	mtext('Generations',side=2,at=mean(ylim),line=2,cex=largefontsize/smallfontsize)
}


dev.off()

#################
##### effect of step preference function

setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_by_mode.R')
source('saveit.R')

trait_chunk_num = 301
sigmay2 = 2
sigmax2 = 0.8
pf = 1
pm = 1 
rho = 0
minweight = 10^(-320)
mut_prob = 0.01

steps = 5000
store = 1
source('range_setup.R')

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.55,0.43,0.02,0.15)
omarg = c(0.05,1,0.35,0.0)

width = 6.5
height = 4

load('step_pref_fun_equilibrium.Rdata')

# pdf('/Users/eleanorbrush/Desktop/effect_of_step_function.pdf',width=width,height=height,family=fontfamily)

par(mfrow=c(2,2),ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))

for(i in 1:2){
	sigma2 = sigma2_vals[c(2,6)[i]]
		
	continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
	fixed_weight = continuous_weight/sum(continuous_weight)
	
	xlim = c(-4,4)
	ylim=c(0,c(.15,.15)[i])
	plot(mrange+1,fixed_weight,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
	mtext('y-x',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
	mtext('Preference',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){		
			k1 = k1_vals[j]
			k2 = k2_vals[j]
			
			chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))
	
			n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
			n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]
			
			s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
			s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]
			
			m = rbind(n,s,c(1,0,0))
			
			v = c(1,sigma2,0)
			p = solve(m,v)
			p[1] = 0
			
			if(length(which(p<0))==0 && p[3]>=p[2]){
			
				fixed_weight = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
				fixed_weight = c(fixed_weight,rev(fixed_weight[1:(length(fixed_weight)-1)]))
				lines(mrange+1,fixed_weight,col=col_vec[j+1],lwd=lwd)
			}
	}
}

for(i in 1:2){
	plot(sigma2_vals,var_mat[i,,1],ylim=c(0,round(max(var_mat[i,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='')
	mtext(expression(sigma^2),side=1,line=2,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste(sigma[x]^2, "*")),expression(paste(sigma[y]^2, "*")))[[i]],side=2,line=1.7,at=(round(max(var_mat[i,,],na.rm=TRUE),1)+0.1)/2,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+1],t='o',lwd=lwd)}
	}
}

# dev.off()

##### effect of step trait distributions

col_vec = brewer.pal(9,'Set1')[-c(6,7)]
lwd = 2
marg = c(0.45,0.43,0.02,0.15)
omarg = c(0.03,1,0.35,0.0)

width = 6.5
height = 7.5

pdf('/Users/eleanorbrush/Desktop/effect_of_step_distribution.pdf',width=width,height=height,family=fontfamily)

par(ps=smallfontsize,mai=marg,oma=omarg,mgp=c(3,0.7,0))
layout(matrix(1:8,ncol=2))

load('step_song_equilibrium.Rdata')

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

xlim = c(-4,4)
ylim=c(0,.12)
plot(mrange+1,m_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
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


for(i in 1:2){
	plot(sigma2_vals,var_mat[i,,1],ylim=c(0,round(max(var_mat[i,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='')
	mtext(expression(sigma^2),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste(sigma[x]^2, "*")),expression(paste(sigma[y]^2, "*")))[[i]],side=2,line=1.7,at=(round(max(var_mat[i,,],na.rm=TRUE),1)+0.1)/2,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+1],t='o',lwd=lwd)}
	}
}

xlim = c(-4,4)
ylim=c(0,.035)
k=6
plot(mrange+1,equilibrium[[1,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[1,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[1,k,j+1]],col=col_vec[j+1],lwd=lwd)
		}
}

load('step_pref_equilibrium.Rdata')

f_init = dnorm(mrange,mmin,sqrt(sigmay2))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)

xlim = c(-4,4)
ylim=c(0,.12)
plot(mrange+1,f_init,col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Preference',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
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
			
			lines(mrange+1,f_init,col=col_vec[j+4],lwd=lwd)
		}
}


for(i in 1:2){
	plot(sigma2_vals,var_mat[i,,1],ylim=c(0,round(max(var_mat[i,,],na.rm=TRUE),1)+.1),t='o',lwd=lwd,col=col_vec[1],xlab='',ylab='')
	mtext(expression(sigma^2),side=1,line=2.3,at=mean(sigma2_vals),cex=largefontsize/smallfontsize)
	mtext(list(expression(paste(sigma[x]^2, "*")),expression(paste(sigma[y]^2, "*")))[[i]],side=2,line=1.7,at=(round(max(var_mat[i,,],na.rm=TRUE),1)+0.1)/2,cex=largefontsize/smallfontsize)
	for(j in 1:length(k1_vals)){
		if(length(which(!is.na(var_mat[i,,j+1])))!=0){
			points(sigma2_vals,var_mat[i,,j+1],col=col_vec[j+4],t='o',lwd=lwd)}
	}
}

xlim = c(-4,4)
ylim=c(0,.105)
k=6
plot(mrange+1,equilibrium[[1,k,1]],col=col_vec[1],t='l',xlim=xlim,lwd=lwd,ylim=ylim,xlab='',ylab='')
mtext('Song',side=1,line=2,at=mean(xlim),cex=largefontsize/smallfontsize)
mtext('Frequency',side=2,line=2,at=mean(ylim),cex=largefontsize/smallfontsize)
for(j in 1:length(k1_vals)){					
		if(!is.na(equilibrium[[1,k,j+1]][1])){			
			lines(mrange+1,equilibrium[[1,k,j+1]],col=col_vec[j+4],lwd=lwd)
		}
}

dev.off()


# source('step_function_example.R')
# source('peak_example.R')