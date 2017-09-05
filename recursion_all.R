Q_fun <- function(sigmax2,sigmay2,cov=rho*sqrt(sigmax2*sigmay2)){
	Q = sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2+2*cov/(sigma2+sigmax2)+1
	return(Q)
}
Q_fun_pref <- function(sigmax2,sigmay2,cov){
	if(is.nan(cov) || round(cov,10)==0){
		rho = 0
	} else {
	rho = cov/sqrt(sigmax2*sigmay2)}
	Q = (sigma2*(sigma2+sigmax2)+(1-rho^2)*sigmax2*(sigma2+sigmax2)+rho^2*sigmax2*sigmay2)/(sigma2+sigmax2)^2+2*cov/(sigma2+sigmax2)+1
	return(Q)
}


recursion_all <- function(sigmax2_init,sigmay2_init,sigma2,rho_init){
	cov_init = rho_init*sqrt(sigmax2_init*sigmay2_init)
	Q_init = Q_fun(sigmax2_init,sigmay2_init,cov_init)
	variables_mat = array(NA,dim=c(9,5,steps+1)) #sigmax2, sigmay2, cov, rho, Q
	variables_mat[,,1]=matrix(c(sigmax2_init,sigmay2_init,cov_init,rho_init,Q_init),nrow=9,ncol=5,byrow=TRUE)
	
	for(t in 1:steps){
		#mode 1: song learned from father, preference genetic
		m = 1
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 1/2*cov*(sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+1/2*sigmax2*sigmay2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else {
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 2: song learned from father, preference imprinted from father
		m = 2
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		cov_new  = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 3: song learned from father, preference learned from mother
		m = 3
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = sigmay2
		cov_new  = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 4: song learned obliquely, preference genetic
		m = 4
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = sigmax2
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 0 
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 5: song learned obliquely, preference imprinted from father
		m = 5
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2
		sigmay2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		cov_new = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 6: song learned obliquely, preference learned from mother
		m = 6
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2
		sigmay2_new = sigmay2
		cov_new = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 7: song genetic, preference genetic
		m = 7
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 1/4*cov*Q+(1-rho^2)/4*sigmax2*sigmay2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 8: song genetic, preference imprinted from father
		m = 8 
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmax2*((sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2)
		cov_new = 1/2*sigmax2*((sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2)+1/2*cov*sigmax2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 9: song genetic, preference learned from mother
		m = 9
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmay2
		cov_new = 1/2*sigmax2*sigmay2/(sigma2+sigmax2)+1/2*cov
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
	}
	return(variables_mat)
}


recursion_all_new_numbers <- function(sigmax2_init,sigmay2_init,sigma2,rho_init){
	cov_init = rho_init*sqrt(sigmax2_init*sigmay2_init)
	Q_init = Q_fun(sigmax2_init,sigmay2_init,cov_init)
	variables_mat = array(NA,dim=c(12,5,steps+1)) #sigmax2, sigmay2, cov, rho, Q
	variables_mat[,,1]=matrix(c(sigmax2_init,sigmay2_init,cov_init,rho_init,Q_init),nrow=12,ncol=5,byrow=TRUE)
	
	for(t in 1:steps){
		#mode 6: song learned from father, preference genetic
		m = 6
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 1/2*cov*(sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2+1/2*sigmax2*sigmay2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else {
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 9: song learned from father, preference imprinted from father
		m = 9
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		cov_new  = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 3: song learned from father, preference learned from mother
		m = 3
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = sigmay2
		cov_new  = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 4: song learned obliquely, preference genetic
		m = 4
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = sigmax2
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 0 
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 7: song learned obliquely, preference imprinted from father
		m = 7
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2
		sigmay2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		cov_new = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 1: song learned obliquely, preference learned from mother
		m = 1
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2
		sigmay2_new = sigmay2
		cov_new = 0
		rho_new = 0
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 5: song genetic, preference genetic
		m = 5
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		Q_pref = Q_fun_pref(sigmax2,sigmay2,cov)
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = 1/4*sigmay2*Q_pref
		cov_new = 1/4*cov*Q+(1-rho^2)/4*sigmax2*sigmay2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 8: song genetic, preference imprinted from father
		m = 8 
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmax2*((sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2)
		cov_new = 1/2*sigmax2*((sigma2*(sigma2+sigmax2)+sigmax2*sigmay2)/(sigma2+sigmax2)^2)+1/2*cov*sigmax2/(sigma2+sigmax2)
		if(is.nan(cov_new) || round(cov_new,10)<1e-13 ){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 2: song genetic, preference learned from mother
		m = 2
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmay2
		cov_new = 1/2*sigmax2*sigmay2/(sigma2+sigmax2)+1/2*cov
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 10: song from oblique male, preference learned from oblique male
		m = 10
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2
		sigmay2_new = sigmax2
		cov_new = 0
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 11: song genetic, preference learned from oblique male
		m = 11
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2+1)
		sigmay2_new = sigmax2
		cov_new = 0
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		#mode 12: song from father, preference learned from oblique male
		m = 12
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = sigmax2*(sigma2/(sigma2+sigmax2)+sigmax2*sigmay2/(sigma2+sigmax2)^2)
		sigmay2_new = sigmax2
		cov_new = 0
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
	}
	return(variables_mat)
}

recursion_test <- function(sigmax2_init,sigmay2_init,sigma2,rho_init){
	cov_init = rho_init*sqrt(sigmax2_init*sigmay2_init)
	Q_init = Q_fun(sigmax2_init,sigmay2_init,cov_init)
	variables_mat = array(NA,dim=c(2,5,steps+1)) #sigmax2, sigmay2, cov, rho, Q
	variables_mat[,,1]=matrix(c(sigmax2_init,sigmay2_init,cov_init,rho_init,Q_init),nrow=2,ncol=5,byrow=TRUE)
	
	for(t in 1:steps){
		
		#mode 9: song genetic, preference learned from mother
		m = 1
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmay2
		cov_new = 1/2*sigmax2*sigmay2/(sigma2+sigmax2)+1/2*cov
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
		
		
		#mode 10: song genetic, preference learned from mother
		m = 2
		sigmax2 = variables_mat[m,1,t]
		sigmay2 = variables_mat[m,2,t]
		cov = variables_mat[m,3,t]
		rho = variables_mat[m,4,t]
		Q = variables_mat[m,5,t]
		sigmax2_new = 1/4*sigmax2*Q
		sigmay2_new = sigmay2
		cov_new = 0
		if(is.nan(cov_new) || round(cov_new,10)<1e-13){
			rho_new = 0
		} else if(log(cov_new)>100 && log(sigmax2_new*sigmay2_new)>100) {
			rho_new = 1} else { 
		rho_new = cov_new / sqrt(sigmax2_new*sigmay2_new)}
		Q_new = Q_fun(sigmax2_new,sigmay2_new,cov_new)
		variables_mat[m,,t+1]=c(sigmax2_new,sigmay2_new,cov_new,rho_new,Q_new)
	}
	return(variables_mat)
}



