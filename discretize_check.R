center = which.max(v)
a = floor(trait_chunk_num/pref_chunk_num)
m = matrix(c(a,a+1,1,1),nrow=2,byrow=TRUE)
X = solve(m,c(trait_chunk_num,pref_chunk_num))
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


# # hold_vec = (c(rev(seq(center-floor(odd_chunk_length/2),1,by=-side_chunk_length)),seq(center+floor(odd_chunk_length/2),Nf,by=side_chunk_length)))
# chunk_vec = c(rep(1,hold_vec[1]-1),as.vector(unlist(apply(matrix(2:length(hold_vec),nrow=1),2,function(x) rep(x,each=diff(hold_vec)[x-1])))),rep(length(hold_vec)+1,Nf-hold_vec[length(hold_vec)]+1))
# chunk_vec = chunk_vec - min(chunk_vec)+1

print(c(length(chunk_vec),trait_chunk_num))
print(c(max(unique(chunk_vec)),pref_chunk_num))
print(c(chunk_vec[center],ceiling(pref_chunk_num/2)))
print(chunk_vec)