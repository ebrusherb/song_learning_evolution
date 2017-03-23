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
mut_prob = 0.0

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

k1_vals = c(35,25,45)
k2_vals = c(35,35,45)

x1 = length(k1_vals)

equilibrium = as.list(1:((x1+1)*xs*3))
dim(equilibrium) = c(3,xs,x1+1)

for(i in 1:xs){
	sigma2 = sigma2_vals[i]
	
	continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
	fixed_weight = continuous_weight/sum(continuous_weight)
		
	f_init = dnorm(frange,fmin,sqrt(sigmay2))
	f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
	f_init = f_init/sum(f_init)
	
	p1 = dynamics_mode3()	
	p2 = dynamics_mode5()	
	# p3 = dynamics_mode9()
	equilibrium[[1,i,1]] = p1$Pm[,steps]
	equilibrium[[2,i,1]] = p2$Pf[,steps]
	# equilibrium[[3,i,1]] = p3$Pm[,steps]
	
	for(j in 1:x1){
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
			
			p1 = dynamics_mode3()	
			p2 = dynamics_mode5()	
			# p3 = dynamics_mode9()
			equilibrium[[1,i,j+1]] = p1$Pm[,steps]
			equilibrium[[2,i,j+1]] = p2$Pf[,steps]
			# equilibrium[[3,i,j+1]] = p3$Pm[,steps]
		# plot(mrange,fixed_weight,main=c(sigma2,j));points(mrange,continuous_weight/sum(continuous_weight),col='red')
		} else{
			# print(c(sigma2,k1,k2))
			equilibrium[[1,i,j+1]] = NA
			equilibrium[[2,i,j+1]] = NA
			# equilibrium[[3,i,j+1]] = NA
		}
	}
}

var_mat = array(NA,dim(equilibrium))

for(i in 1:xs){
	for(j in 1:(x1+1)){
		for(k in 1:2){
			if(!is.na(equilibrium[[k,i,j]][1])){
				v = sum((mrange+1)^2*(equilibrium[[k,i,j]]))
				var_mat[k,i,j] = v }
		}
	}
}

saveit(sigma2_vals=sigma2_vals,sigmax2=sigmax2,sigmay2=sigmay2,k1_vals=k1_vals,k2_vals=k2_vals,equilibrium=equilibrium,var_mat=var_mat,file='/Users/eleanorbrush/Documents/research/song_learning_evolution/step_pref_equilibrium.Rdata')
