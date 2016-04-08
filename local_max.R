local_max <- function(v){
	l = length(v)
	m = c()
	if(v[1]>v[2]){ m=c(m,1)}
	for(i in 2:(l-1)){
		if(v[i]>v[i-1] && v[i]>v[i+1]){m=c(m,i)}
	}
	if(v[l]>v[l-1]){ m =c(m,l)}
	return(m)
}