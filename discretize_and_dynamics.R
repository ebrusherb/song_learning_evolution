
discretize <- function(v,pref_chunk_num = trait_chunk_num){
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


## ---- dynamics -----------------------

dynamics <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_adults
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_adults
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}