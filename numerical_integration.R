# setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('ind2sub.R')
source('glue.R')
source('int.R')
# source('range_setup.R')
source('recursion_all.R')
source('image.scale.R')
library(RColorBrewer)

discretize <- function(v,pref_chunk_num = trait_chunk_num){
	center = which.max(v)
	a = floor(trait_chunk_num/pref_chunk_num)
	m = matrix(c(a,a+1,1,1),nrow=2,byrow=TRUE)
	X = round(solve(m,c(trait_chunk_num,pref_chunk_num)))
	if(a%%2!=0){
		odd_chunk_length = a
		even_chunk_length = a+1
		num_odd_chunk = X[1]
		num_even_chunk = X[2]
	} else{
		odd_chunk_length = a+1
		even_chunk_length = a
		num_odd_chunk = X[2]
		num_even_chunk = X[1]
		}
	if(num_even_chunk>0){
	chunk_vec = c(rep(1:(num_even_chunk/2),each=even_chunk_length),rep((1:num_odd_chunk)+(num_even_chunk/2),each=odd_chunk_length),rep(1:(num_even_chunk/2)+(num_odd_chunk+num_even_chunk/2),each=even_chunk_length))
	} else{
		chunk_vec = rep(1:num_odd_chunk,each=odd_chunk_length)
		}
	chunk_vec = chunk_vec - min(chunk_vec)+1
	chunk=list(chunk_vec)
	averaged = aggregate(v,by=chunk,FUN=mean)$x
	discretized = averaged[chunk_vec]
	return(discretized)
}

## ---- dynamics -----------------------

dynamics <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_adults
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_adults
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

#######
steps = 1000
store = 1
pm = 1
pf = 1
mut_prob = 0 
rho_init = 0

sigma2 = 0.01
sigmay2_init = 1.1
sigmax2_init = 0.5

sigma2_vals = c(0.005,0.01,1)
x_sigma2 = length(sigma2_vals)
trait_chunk_num_vals = (round(c(3,3^seq(2,5.3,length.out=20))))
trait_chunk_num_vals[which(trait_chunk_num_vals%%2==0)] = trait_chunk_num_vals[which(trait_chunk_num_vals%%2==0)]+1 
trait_chunk_num_vals = unique(trait_chunk_num_vals)
x_trait = length(trait_chunk_num_vals)
pref_chunk_num_vals = trait_chunk_num_vals
x_pref = length(pref_chunk_num_vals)

# equilibrium = array(NA,dim=c(x_sigma2,x_trait_step,x_step_width,2,Nm))
equilibrium = vector('list',x_sigma2*x_trait*x_pref*2)
dim(equilibrium) = c(x_sigma2,x_trait,x_pref,2)

for(i in 1:x_sigma2){
	for(j in 1:x_trait){
		for(k in 1:j){
			sigma2 = sigma2_vals[i]
			trait_chunk_num = trait_chunk_num_vals[j]
			pref_chunk_num = pref_chunk_num_vals[k]
			source("range_setup.R")
			
			f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
			f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
			f_init = f_init/sum(f_init)
			f_init = pf*f_init+(1-pf)*rev(f_init)
			m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
			m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
			m_init = m_init / sum(m_init)
			m_init = pf*m_init+(1-pf)*rev(m_init)
			
			continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
			minweight = 10^(-320)
			# minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
			# continuous_weight[continuous_weight==0] = minweight
			
			discrete_weight = discretize(continuous_weight,pref_chunk_num)
			
			# fixed_weight = continuous_weight/sum(continuous_weight)
			fixed_weight = discrete_weight/sum(discrete_weight)
			
			t = steps
			analytical = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
			recurs = dnorm(mrange,mean=-1,sd=sqrt(analytical[3,1,t]))
			recurs = recurs / sum(recurs)
			pop_dens = dynamics()
			Pm = pop_dens$Pm[,t]
			
			equilibrium[[i,j,k,1]]=recurs
			equilibrium[[i,j,k,2]]=Pm
		}
	}
}

save(sigma2_vals,trait_chunk_num_vals,pref_chunk_num_vals,sigmay2_init,steps,store,sigmax2_init,equilibrium,file='effect_of_trait_step_and_step_width.Rdata')


ex_mat = array(NA,dim=c(x_sigma2,x_trait,x_pref,2))
var_mat = array(NA,dim=c(x_sigma2,x_trait,x_pref,2))
for(i in 1:x_sigma2){
	for(j in 1:x_trait){
		for(k in 1:j){
			sigma2 = sigma2_vals[i]
			trait_chunk_num = trait_chunk_num_vals[j]
			pref_chunk_num = pref_chunk_num_vals[k]
			source("range_setup.R")
			for(l in 1:2){
			ex = sum(mrange*equilibrium[[i,j,k,l]])
			vx = sum((mrange-(ex))^2*equilibrium[[i,j,k,l]])
			ex_mat[i,j,k,l] = ex
			var_mat[i,j,k,l] = vx
			}
		}
	}
}


layout(matrix(1:(2*x_sigma2),ncol=x_sigma2))
for(i in 1:x_sigma2){
m=min(var_mat[,,,],na.rm=TRUE);
M=max(var_mat[,,,],na.rm=TRUE);
breaks=c(seq(0,0.05,length.out=10),seq(0.1,M,length.out=10));
# breaks= c(m,10^seq(-5,log(M,base=10),length.out=30))
cols=heat.colors(length(breaks)-1);
image(var_mat[i,,,1],breaks=breaks,col=cols,xlab='Trait chunk #',ylab='Pref fun chunk #',xaxt='n',yaxt='n');
axis(1,at=seq(0,1,length.out=x_trait),labels=trait_chunk_num_vals)
axis(2,at=seq(0,1,length.out=x_pref),labels=pref_chunk_num_vals)
abline(0,1)
# image.scale(var_mat[i,,,1],horiz=FALSE,breaks=breaks,col=cols);
# # m=min(var_mat[i,,,2],na.rm=TRUE);
# M=max(var_mat[i,,,2],na.rm=TRUE);
# breaks=seq(m,M,length.out=50);
# cols=heat.colors(length(breaks)-1);
image(var_mat[i,,,2],breaks=breaks,col=cols,xlab='Trait chunk #',ylab='Pref fun chunk #',xaxt='n',yaxt='n');
axis(1,at=seq(0,1,length.out=x_trait),labels=trait_chunk_num_vals)
axis(2,at=seq(0,1,length.out=x_pref),labels=pref_chunk_num_vals)
abline(0,1)
# image.scale(var_mat[i,,,2],horiz=FALSE,breaks=breaks,col=cols)
}
# layout(matrix(1:2,ncol=2))
# i=1
# to_plot = (var_mat[i,,,1]-var_mat[i,,,2])/var_mat[i,,,1]
# m = min(to_plot,na.rm=TRUE)
# M = max(to_plot,na.rm=TRUE)
# M = max(abs(m),M)
# breaks = seq(-M,M,length.out=21)
# cols = colorRampPalette(brewer.pal(9, "RdBu"))(length(breaks)-1)
# image(trait_chunk_num_vals,pref_chunk_num_vals,to_plot,breaks=breaks,col=cols)
# image.scale(to_plot,breaks=breaks,col=cols,horiz=FALSE)

# # # 
# # load("/Users/eleanorbrush/Documents/research/song_learning_evolution/effect_of_sigma2_and_step_width_01_05_2016.Rdata")
# # x_sigma2 = length(sigma2_vals)
# # x_sigmax2 = length(sigmax2_init_vals)
# # x_step_width = length(step_width_vals)
# # ex_mat = array(NA,dim=c(x_sigma2,x_sigmax2,x_step_width,2))
# # var_mat = array(NA,dim=c(x_sigma2,x_sigmax2,x_step_width,2))
# # trait_step = 0.1
# # int_step = trait_step		
# # source("range_setup.R")
# # for(i in 1:x_sigma2){
	# # for(j in 1:x_sigmax2){
		# # for(k in 1:x_step_width){			
			# # for(l in 1:2){
			# # ex = sum(mrange*equilibrium[i,j,k,l,])
			# # vx = sum((mrange-(ex))^2*equilibrium[i,j,k,l,])
			# # ex_mat[i,j,k,l] = ex
			# # var_mat[i,j,k,l] = vx
			# # }
		# # }
	# # }
# # }

# # layout(matrix(1:2,ncol=2))
# # m = min(var_mat[,j,,1]-var_mat[,j,,2])
# # M = max(var_mat[,j,,1]-var_mat[,j,,2])
# # M = max(abs(m),M)
# # breaks = seq(-M,M,length.out=21)
# # cols = colorRampPalette(brewer.pal(9, "RdBu"))(length(breaks)-1)
# # image(sigma2_vals,step_width_vals,var_mat[,j,,1]-var_mat[,j,,2],breaks=breaks,col=cols)
# # image.scale(var_mat[,j,,1]-var_mat[,j,,2],breaks=breaks,col=cols,horiz=FALSE)

# # #####
# # steps = 20000
# # sigma2 = sigma2_vals[1]
# # sigmay2_init = 2
# # sigmax2_init = sigmax2_init_vals[5]
# # step_width = 4

# # f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
# # f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
# # f_init = f_init/sum(f_init)
# # f_init = pf*f_init+(1-pf)*rev(f_init)
# # m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
# # m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
# # m_init = m_init / sum(m_init)
# # m_init = pf*m_init+(1-pf)*rev(m_init)

# # continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
# # minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
# # continuous_weight[continuous_weight==0] = minweight
# # discrete_weight = discretize(continuous_weight,step_width)

# # # fixed_weight = continuous_weight/sum(continuous_weight)
# # fixed_weight = discrete_weight/sum(discrete_weight)

# # recurs = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
# # recurs = dnorm(mrange,mean=-1,sd=sqrt(recurs[3,1,t]))
# # recurs = recurs / sum(recurs)
# # pop_dens = dynamics()

# # t=steps
# # plot(recurs,col='red',t='l')
# # lines(pop_dens$Pm[,t],t='l')

# # #######
# sigmay2_init = 2
# sigmax2_init = 0.5
# sigma2 = 1
# mut_prob = 0
# pf = 1
# pm = 1
# steps = 50
# store = 1

# trait_chunk_num = 109
# source("range_setup.R")

# continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
# minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
# continuous_weight[continuous_weight==0] = minweight

# fixed_weight = continuous_weight/sum(continuous_weight)

# pref_chunk_num = 3
# fixed_weight = discretize(continuous_weight,pref_chunk_num)
# fixed_weight = fixed_weight / sum(fixed_weight)


# f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
# f_init = f_init/sum(f_init)
# f_init = pf*f_init+(1-pf)*rev(f_init)
# j = midpt
# # j = min(which(mrange>0))
# s = sign(midpt-j+0.5)
# toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
# pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
# hold = f_init
# hold[toreplace[1]:toreplace[2]] = f_init[pull[1]:pull[2]]
# if(toreplace[1]>1){
	# hold[1:(toreplace[1]-1)] = f_init[(pull[2]+1):Nf]}
# f_init = hold

# m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
# m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
# m_init = m_init / sum(m_init)
# m_init = pf*m_init+(1-pf)*rev(m_init)

# t=steps
# recurs = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
# recurs = dnorm(mrange,mean=-1,sd=sqrt(recurs[3,1,t]))
# recurs = recurs / sum(recurs)
# pop_dens = dynamics()
# Pm = pop_dens$Pm[,steps]
# # print(round(Pm,5))
# plot(mrange,recurs)
# points(mrange,Pm,col='red')

############## measuring variance of discretized preference function

pref_var_mat = array(NA,dim=c(x_sigma2,x_trait,x_pref))
pref_ex_mat = array(NA,dim=c(x_sigma2,x_trait,x_pref))
y_var_mat = array(NA,dim=c(x_trait,x_pref))
y_ex_mat = array(NA,dim=c(x_trait,x_pref))
x_var_mat = array(NA,dim=c(x_trait,x_pref))
x_ex_mat = array(NA,dim=c(x_trait,x_pref))

for(i in 1:x_sigma2){
	for(j in 1:x_trait){
		for(k in 1:j){
			sigma2 = sigma2_vals[i]
			trait_chunk_num = trait_chunk_num_vals[j]
			pref_chunk_num = pref_chunk_num_vals[k]
			source("range_setup.R")
			
			f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
			f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
			f_init = f_init/sum(f_init)
			f_init = pf*f_init+(1-pf)*rev(f_init)
			m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
			m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
			m_init = m_init / sum(m_init)
			m_init = pf*m_init+(1-pf)*rev(m_init)
			
			continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) #female preference function
			# minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
			# continuous_weight[continuous_weight==0] = minweight
			discrete_weight = discretize(continuous_weight,pref_chunk_num)
			
			# fixed_weight = continuous_weight/sum(continuous_weight)
			fixed_weight = discrete_weight/sum(discrete_weight)
			
			ex = sum(frange*f_init)
			vx = sum((frange-ex)^2*f_init)
			y_ex_mat[j,k] = ex
			y_var_mat[j,k] = vx
			
			ex = sum(mrange*m_init)
			vx = sum((mrange-ex)^2*m_init)
			x_ex_mat[j,k] = ex
			x_var_mat[j,k] = vx
			
			ex = sum(mrange*fixed_weight)
			vx = sum((mrange-(ex))^2*fixed_weight)
			pref_ex_mat[i,j,k] = ex
			pref_var_mat[i,j,k] = vx
			}
	}		
}

layout(matrix(1:4,ncol=4))
i=1
image(pref_var_mat[i,,])
image(y_var_mat)
image(x_var_mat)
image(var_mat[i,,,2])

i=1
plot(pref_var_mat[i,,],var_mat[i,,,2])
