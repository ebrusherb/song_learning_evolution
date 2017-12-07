init_conds  <- function(song,pref,func,sigmax2,sigmay2,sigma2,k1=NA,k2=NA){
	if(song=='norm'){
		m_init = dnorm(mrange,mmin,sqrt(sigmax2))
		m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
		m_init = m_init/sum(m_init)
	}else if(song=='step'){
		chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

		n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
		n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]
		
		s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
		s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]
		
		m = rbind(n,s,c(1,0,0))
		
		v = c(1,sigmax2,0)
		p = solve(m,v)
		p[1] = 0
		
		if(sum(p<0)==0){
		m_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
		m_init = c(m_init,rev(m_init[1:(length(m_init)-1)]))
		} else{m_init = rep(NA,Nm)}
	}else{
		m_init = rep(NA,Nm)
	}
	
	if(pref=='norm'){
		f_init = dnorm(mrange,-1,sqrt(sigmay2))
		f_init = f_init/sum(f_init)
	}else if(pref=='step'){
		chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

		n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
		n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]
		
		s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
		s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]
		
		m = rbind(n,s,c(1,0,0))
		
		v = c(1,sigmay2,0)
		p = solve(m,v)
		p[1] = 0
		if(sum(p<0)==0){
		f_init = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
		f_init = c(f_init,rev(f_init[1:(length(f_init)-1)]))
		} else{f_init = rep(NA,Nm)}
	}else{
		f_init = rep(NA,Nm)
	}
	
	if(func=='norm'){
		continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
		fixed_weight = continuous_weight/sum(continuous_weight)
	}else if(func=='step'){
		chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

		n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
		n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]
		
		s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
		s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]
		
		m = rbind(n,s,c(1,0,0))
		
		v = c(1,sigma2,0)
		p = solve(m,v)
		p[1] = 0
		if(sum(p<0)==0){		
		fixed_weight = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
		fixed_weight = c(fixed_weight,rev(fixed_weight[1:(length(fixed_weight)-1)]))
		} else{fixed_weight = rep(NA,Nm)}
	}else{
		fixed_weight = rep(NA,Nm)
	}
	return(data.frame(m_init = m_init,f_init = f_init,fixed_weight =fixed_weight))
}


P_init <- function(){

	ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
	m_init = ic$m_init
	f_init = ic$f_init	
	
	Pm = array(0,c(Nm,Nm))
	if(song=='step'){
		for(i in 1:Nm){
			Pm[i,] = m_init[i]*dnorm(mrange,-1+rho*sqrt(sigmay2/sigmax2)*(mrange[i]+1),sqrt((1-rho^2)*sigmay2))/sum(dnorm(mrange,-1+rho*sqrt(sigmay2/sigmax2)*(mrange[i]+1),sqrt((1-rho^2)*sigmay2)))
		}
	}else{
		for(i in 1:Nm){
			Pm[,i] = f_init[i]*dnorm(mrange,-1+rho*sqrt(sigmax2/sigmay2)*(mrange[i]+1),sqrt((1-rho^2)*sigmax2))/sum(dnorm(mrange,-1+rho*sqrt(sigmax2/sigmay2)*(mrange[i]+1),sqrt((1-rho^2)*sigmax2)))
		}
	}
	return(Pm)
}

P_init_pearson <- function(){

	ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
	m_init = ic$m_init
	f_init = ic$f_init	
	
	Pm = array(0,c(Nm,Nm))
	if(song=='step'){
		for(i in 1:Nm){
			Pm[i,] = m_init[i]*dnorm(mrange,-1+rho*sqrt(sigmay2/sigmax2)*(mrange[i]+1),sqrt((1-rho^2)*sigmay2))/sum(dnorm(mrange,-1+rho*sqrt(sigmay2/sigmax2)*(mrange[i]+1),sqrt((1-rho^2)*sigmay2)))
		}
	}else{
		for(i in 1:Nm){
			Pm[,i] = f_init[i]*dnorm(mrange,-1+rho*sqrt(sigmax2/sigmay2)*(mrange[i]+1),sqrt((1-rho^2)*sigmax2))/sum(dnorm(mrange,-1+rho*sqrt(sigmax2/sigmay2)*(mrange[i]+1),sqrt((1-rho^2)*sigmax2)))
		}
	}
	return(Pm)
}
