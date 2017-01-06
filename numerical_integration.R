## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)

# # num_cores <- detectCores()-1
# num_cores <-10
# cl <-makeCluster(num_cores)
# registerDoParallel(cl)

setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('ind2sub.R')
source('glue.R')
source('int.R')
source('range_setup.R')
source('recursion_all.R')
source('image.scale.R')
library(RColorBrewer)

discretize <- function(v,step_width=0.5){
	center = which.max(v)
	chunk_length = ceiling(step_width/step)
	hold_vec = (c(rev(seq(center-floor(chunk_length/2),1,by=-chunk_length)),seq(center+ceiling(chunk_length/2)-1,Nf,by=chunk_length)))
	chunk_vec = c(rep(1,hold_vec[1]-1),rep(2:(length(hold_vec)),each=chunk_length),rep(length(hold_vec)+1,Nf-hold_vec[length(hold_vec)]))
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
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
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
sigmay2_init = 2
sigmax2_init = 0.01
step_width = 1

sigma2_vals = seq(0.01,sigmay2_init,length.out=10)
x_sigma2 = length(sigma2_vals)
sigmax2_init_vals = 0.5
x_sigmax2 = length(sigmax2_init_vals)
step_width_vals = seq(step,2,length.out=10)
x_step_width = length(step_width_vals)

equilibrium = array(NA,dim=c(x_sigma2,x_sigmax2,x_step_width,2,Nm))

for(i in 1:x_sigma2){
	for(j in 1:x_sigmax2){
		for(k in 1:x_step_width){
			sigma2 = sigma2_vals[i]
			sigmax2_init = sigmax2_init_vals[j]
			step_width = step_width_vals[k]
			
			f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
			f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
			f_init = f_init/sum(f_init)
			f_init = pf*f_init+(1-pf)*rev(f_init)
			m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
			m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
			m_init = m_init / sum(m_init)
			m_init = pf*m_init+(1-pf)*rev(m_init)
			
			continuous_weight = dnorm(mrange,mean=frange[midpt],sd=sqrt(sigma2)) #female preference function
			minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
			continuous_weight[continuous_weight==0] = minweight
			discrete_weight = discretize(continuous_weight,step_width)
			
			# fixed_weight = continuous_weight/sum(continuous_weight)
			fixed_weight = discrete_weight/sum(discrete_weight)
			
			t = steps
			recurs = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
			recurs = dnorm(mrange,mean=-1,sd=sqrt(recurs[3,1,t]))
			recurs = recurs / sum(recurs)
			pop_dens = dynamics()
			Pm = pop_dens$Pm[,t]
			
			equilibrium[i,j,k,1,]=recurs
			equilibrium[i,j,k,2,]=Pm
		}
	}
}

# save(sigma2_vals,sigmax2_init_vals,step_width_vals,sigmay2_init,steps,equilibrium,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/effect_of_sigma2_and_step_width.Rdata')


ex_mat = array(NA,dim=c(x_sigma2,x_sigmax2,x_step_width,2))
var_mat = array(NA,dim=c(x_sigma2,x_sigmax2,x_step_width,2))
for(i in 1:x_sigma2){
	for(j in 1:x_sigmax2){
		for(k in 1:x_step_width){
			for(l in 1:2){
			ex = sum(mrange*equilibrium[i,j,k,l,])
			vx = sum((mrange-(ex))^2*equilibrium[i,j,k,l,])
			ex_mat[i,j,k,l] = ex
			var_mat[i,j,k,l] = vx
			}
		}
	}
}

layout(matrix(1:2,ncol=2))
m = min(var_mat[,j,,1]-var_mat[,j,,2])
M = max(var_mat[,j,,1]-var_mat[,j,,2])
M = max(abs(m),M)
breaks = seq(-M,M,length.out=21)
cols = colorRampPalette(brewer.pal(9, "RdBu"))(length(breaks)-1)
image(sigma2_vals,step_width_vals,var_mat[,j,,1]-var_mat[,j,,2],breaks=breaks,col=cols)
image.scale(var_mat[,j,,1]-var_mat[,j,,2],breaks=breaks,col=cols,horiz=FALSE)

#####
steps = 20000
sigma2 = sigma2_vals[1]
sigmay2_init = 2
sigmax2_init = sigmax2_init_vals[5]
step_width = 4

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

continuous_weight = dnorm(mrange,mean=frange[midpt],sd=sqrt(sigma2)) #female preference function
minweight = 10^max(floor(log(min(continuous_weight[which(continuous_weight>0)]),base=10)),-320)	
continuous_weight[continuous_weight==0] = minweight
discrete_weight = discretize(continuous_weight,step_width)

# fixed_weight = continuous_weight/sum(continuous_weight)
fixed_weight = discrete_weight/sum(discrete_weight)

recurs = recursion_all(sigmax2_init,sigmay2_init,sigma2,rho_init)
recurs = dnorm(mrange,mean=-1,sd=sqrt(recurs[3,1,t]))
recurs = recurs / sum(recurs)
pop_dens = dynamics()

t=steps
plot(recurs,col='red',t='l')
lines(pop_dens$Pm[,t],t='l')
