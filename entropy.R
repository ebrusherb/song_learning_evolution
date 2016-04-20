#entropy of signals received
ent<-function(v){
	n<-length(v)
	hold<-array(0,c(n,n))
	hold = -v*log(v)
	hold[v==0] = 0
	return(int(hold))
}