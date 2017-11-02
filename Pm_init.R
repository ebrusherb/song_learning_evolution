Pm_init <- function(){

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
