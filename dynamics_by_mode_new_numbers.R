dynamics_mode3 <-function(){
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


dynamics_mode7 <-function(){
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
	Pm_new = Pm_adults/sum(Pf_adults)
	Pf_new = apply(pxy,1,sum)/sum(Pf_adults)	
	nonzero = which(Pf_adults>nonzero_thresh)
	perc = max(abs(range(Pf_new[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_new
		Pf[,t+1] = Pf_new
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_new
			Pf[,(t+1):(steps+1)] = Pf_new
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}


dynamics_mode2 <-function(){
	Pm = matrix(0,Nm,steps+1) #probability of male songs over time
	Pm[,1] = m_init
	
	Pf = array(0,c(Nm,Nf,steps+1)) #probability of female preferences over time
	Pf[,,1] = sapply(matrix(f_init,nrow=1),function(x) x*m_init,simplify=TRUE)
	
	reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))
	
	t = 1
	perc = 0
	while(t <= steps){
		Pm_adults = Pm[,t]
		Pf_adults = Pf[,,t]
		pxy = array(0,c(Nm,Nm,Nf)) #probability of a x / (x,y) pair
		Pf_new = array(0,c(Nm,Nf))
		for(j in 1:Nf){
			s = sign(midpt-j+0.5)
			weight = rep(minweight,Nf)
			toreplace=sort(c((s)%%(Nf+1),j+s*(midpt-1)))
			pull=sort(c((s)%%(Nf+1)+midpt-j,(-s)%%(Nf+1)))
			weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
			z = sum(weight*Pm_adults) #normalization factor
			if(z>z_thresh){
				pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
			}
			song_pairs = unlist(lapply(split(pxy[,,j],reverse_diagonal),sum))
			song_pairs = song_pairs[seq(1,2*Nm-1,by=2)]+1/2*c(song_pairs[seq(2,2*Nm-1,by=2)],2*song_pairs[2*Nm-1])+1/2*c(2*song_pairs[1],song_pairs[seq(2,2*Nm-1,by=2)])
			Pf_new[,j] = song_pairs
		}		
		Pf_adults[,apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
		Pm_new = apply(Pf_new,1,sum)/sum(Pf_adults)
		Pf_new = Pf_new/sum(Pf_adults) 	
		nonzero = which(Pm_new>nonzero_thresh)
		perc = max(abs(range(Pm_new[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
		if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
			Pm[,t+1] = Pm_new
			Pf[,,t+1] = Pf_new
			t = t+1
			} else{
				Pm[,(t+1):(steps+1)] = Pm_new
				Pf[,,(t+1):(steps+1)] = Pf_new
				print(c(t,'problem'))
				t = steps+1
				}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}
