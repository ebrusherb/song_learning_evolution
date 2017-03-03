discretize_constant_variance <- function(v,pref_chunk_num = trait_chunk_num){
	center = which.max(v)
	a = floor(trait_chunk_num/pref_chunk_num)
	m = matrix(c(a,a+1,1,1),nrow=2,byrow=TRUE)
	X = round(solve(m,c(trait_chunk_num,pref_chunk_num)))
	if(a%%2!=0){
		odd_chunk_length = a
		even_chunk_length = a+1
		num_odd_chunk = X[1]
		num_even_chunk = X[2]
	} else{
		odd_chunk_length = a+1
		even_chunk_length = a
		num_odd_chunk = X[2]
		num_even_chunk = X[1]
		}
	if(num_even_chunk>0){
	chunk_vec = c(rep(1:(num_even_chunk/2),each=even_chunk_length),rep((1:num_odd_chunk)+(num_even_chunk/2),each=odd_chunk_length),rep(1:(num_even_chunk/2)+(num_odd_chunk+num_even_chunk/2),each=even_chunk_length))
	} else{
		chunk_vec = rep(1:num_odd_chunk,each=odd_chunk_length)
		}
	chunk_vec = chunk_vec - min(chunk_vec)+1
	chunk=list(chunk_vec)
	averaged = aggregate(v,by=chunk,FUN=mean)$x
	discretized = averaged[chunk_vec]
	return(discretized)
}
