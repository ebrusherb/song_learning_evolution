## ---- integral -------------------------
int<-function(v){ #manual integration
	l=length(v)
	int_step/2*(sum(2*v)-v[1]-v[l])
} 