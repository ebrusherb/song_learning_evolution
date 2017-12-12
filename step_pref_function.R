# setwd('/Users/eleanorbrush/Documents/research/song_learning_evolution')
source('dynamics_by_mode_new_numbers.R')
source('saveit.R')
source('init_conds.R')

trait_chunk_num = 281
minweight = 10^(-320)
mut_prob = 0.0
source('range_setup.R')

sigmay2 = 2
sigmax2 = 1
rho = 0

steps = 30000
store = 1

sigma2_vals = seq(0.1,1.9,length.out=10)

xs = length(sigma2_vals)

k1_vals = c(7,15,35,21,7,21)
k2_vals = c(35,35,35,35,21,21)

x1 = length(k1_vals)

equilibrium = as.list(1:((x1+1)*xs*2))
dim(equilibrium) = c(2,xs,x1+1)

for(i in 1:xs){
	sigma2 = sigma2_vals[i]
	
	ic = init_conds('norm','norm','norm',sigmax2,sigmay2,sigma2,NA,NA)
	m_init = ic$m_init
	f_init = ic$f_init
	fixed_weight = ic$fixed_weight
	
	p1 = dynamics_mode3()	
	p2 = dynamics_mode7()	
	equilibrium[[1,i,1]] = p1$Pm[,steps]
	equilibrium[[2,i,1]] = p2$Pf[,steps]
	
	for(j in 1:x1){
		k1 = k1_vals[j]
		k2 = k2_vals[j]
		
		ic = init_conds('norm','norm','step',sigmax2,sigmay2,sigma2,k1,k2)
		fixed_weight = ic$fixed_weight
				
		if(sum(is.na(fixed_weight))==0){		
			
			p1 = dynamics_mode3()	
			p2 = dynamics_mode7()	
			equilibrium[[1,i,j+1]] = p1$Pm[,steps]
			equilibrium[[2,i,j+1]] = p2$Pf[,steps]			
		} else{
			equilibrium[[1,i,j+1]] = NA
			equilibrium[[2,i,j+1]] = NA
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


saveit(sigma2_vals=sigma2_vals,sigmax2=sigmax2,sigmay2=sigmay2,k1_vals=k1_vals,k2_vals=k2_vals,equilibrium=equilibrium,var_mat=var_mat,mut_prob=mut_prob,file='/homes/ebrush/priv/song_learning_evolution/step_pref_fun_equilibrium.Rdata')
