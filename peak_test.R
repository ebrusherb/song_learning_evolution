# setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
# source('dynamics_by_mode.R')
# source('saveit.R')

# dynamics_mode3_bothmut <-function(){
# Pm = matrix(0,Nm,steps+1) #probability of male songs over time
# Pm[,1] = m_init

# Pf = matrix(0,Nf,steps+1) #probability of female preferences over time
# Pf[,1] = f_init

# t = 1
# perc = 0
# while(t <= steps){
	# Pm_adults = Pm[,t]
	# Pf_adults = Pf[,t]
	# pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	# ### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	# for(j in 1:Nf){
		# s = sign(midpt-j+0.5)
		# weight = rep(minweight,Nf)
		# toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
		# pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
		# weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		# z = sum(weight*Pm_adults) #normalization factor
		# if(z>z_thresh){
			# pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			# }
	# }
	# Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	# Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	# Pf_adults = (1-mut_prob)*Pf_adults + mut_prob/2*c(Pf_adults[2:Nm],0) + 
		# mut_prob/2*c(0,Pf_adults[1:(Nm-1)]) /sum(Pf_adults)
	# Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		# mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	# nonzero = which(Pm_adults>nonzero_thresh)
	# perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	# if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		# Pm[,t+1] = Pm_aftermut
		# Pf[,t+1] = Pf_adults
		# t = t+1
		# } else{
			# Pm[,(t+1):(steps+1)] = Pm_aftermut
			# Pf[,(t+1):(steps+1)] = Pf_adults
			# t = steps+1
			# }
# }
# pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
# return(pop_dens)
# }

# trait_chunk_num = 301
# sigmay2 = 2
# sigmax2 = 0.8
# pf = 1
# pm = 1 
# rho = 0
# minweight = 10^(-320)
# mut_prob = 0.01

# steps = 5000
# store = 1
# source('range_setup.R')

# f_init = dnorm(frange,fmin,sqrt(sigmay2))
# f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
# f_init = f_init/sum(f_init)

# trait_chunk_num = 3

# sigma2_vals = seq(0.1,1.9,length.out=10)

# xs = length(sigma2_vals)

# k1_vals = c(7,7,15,35,7,15,21,21)
# k2_vals = c(11,35,35,35,21,21,21,35)

# k1_vals = c(7,7,15)
# k2_vals = c(35,25,25)

# x1 = length(k1_vals)

# sigma2 = 0.9


	# continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
	# fixed_weight = continuous_weight/sum(continuous_weight)

 # j=3
 
 # k1 = k1_vals[j]
		# k2 = k2_vals[j]
		
		# chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

		# n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
		# n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]
		
		# s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
		# s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]
		
		# m = rbind(n,s,c(1,0,0))
		
		# v = c(1,sigmax2,0)
		# p = solve(m,v)
		# p[1] = 0
		
		# if(length(which(p<0))==0 && p[3]>=p[2] ){

			# m_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
			# m_init = c(m_init,rev(m_init[1:(length(m_init)-1)]))
			
			# p1 = dynamics_mode3_bothmut()	
			# }
			
			
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

m_init = dnorm(mrange,mmin,sqrt(sigmax2))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init/sum(m_init)

trait_chunk_num = 3

sigma2_vals = seq(0.1,1.9,length.out=10)

xs = length(sigma2_vals)

k1_vals = c(7,7,15,35,7,15,21,21)
k2_vals = c(11,35,35,35,21,21,21,35)

k1_vals = c(35,25,55)
k2_vals = c(35,35,45)

sigma2 = 0.9

j=1

	continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
	fixed_weight = continuous_weight/sum(continuous_weight)
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
			
			p1 = dynamics_mode3_bothmut()	}			