delta = 0.00
local_max <- function(v){
	l = length(v)
	m = c()
	if(v[1]>(v[2]+delta)){ m=c(m,1)}
	for(i in 2:(l-1)){
		if(v[i]>(v[i-1]+delta) && v[i]>(v[i+1]+delta)){m=c(m,i)}
	}
	if(v[l]>(v[l-1]+delta)){ m =c(m,l)}
	return(m)
}