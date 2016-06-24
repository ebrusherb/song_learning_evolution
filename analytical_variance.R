sigma_malefath_femalefath<-function(m,f){
	m = c(m,m[1]*(m[1]*f[1]+(m[1]+sigma^2)*sigma^2)/(m[1]+sigma^2)^2)
	f = c(f,m[1]*(m[1]*f[1]+(m[1]+sigma^2)*sigma^2)/(m[1]+sigma^2)^2)
	for(T in 2:5000){
		m = c(m,m[T]*(m[T]*f[T]+(m[T]+sigma^2)*sigma^2)/(m[T]+sigma^2)^2)
		f = c(f,m[T]*(m[T]*f[T]+(m[T]+sigma^2)*sigma^2)/(m[T]+sigma^2)^2)
	}
	return(cbind(m,f))
}

sigma_malefath_femalemoth<-function(m,f){	
	for(T in 1:5000){
		m = c(m,m[T]*(m[T]*f[T]+(m[T]+sigma^2)*sigma^2)/(m[T]+sigma^2)^2)
		f = c(f,f)
	}
	return(cbind(m,f))
}

var_time <- function(Pm){
	T = dim(Pm)[2]
	var_time = matrix()
	for(t in 1:T){
		var_time[t] = sum((mrange-(mmin))^2*Pm[,t])
	}
	return(var_time)
}