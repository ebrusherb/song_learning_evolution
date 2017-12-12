dynamics_mode2 <-function(){
	Pm = matrix(0,Nm,steps+1) #probability of male songs over time
	Pm[,1] = m_init
	
	Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
	Pf[,,1] = P_init()
	
	reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))
	
	t = 1
	perc = 0
	while(t <= steps){
		Pm_adults = Pm[,t]
		Pf_adults = Pf[,,t]
		pxy = array(0,c(Nm,Nm,Nm)) #probability of a x / (x,y) pair
		offspring = array(0,c(Nm,Nm))
		for(j in 1:Nm){
			s = sign(midpt-j+0.5)
			weight = rep(minweight,Nm)
			toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
			pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
			weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
			z = sum(weight*Pm_adults) #normalization factor
			if(z>z_thresh){
				pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
			}
			matches = unlist(lapply(split(pxy[,,j],reverse_diagonal),sum))
			song_pairs = matches[seq(1,2*Nm-1,by=2)]+1/2*c(matches[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,matches[seq(2,2*Nm-1,by=2)])			
			# song_pairs = matches[seq(1,2*Nm-1,by=2)]+c(0,matches[seq(2,Nm-3,by=2)],matches[Nm-1]+matches[Nm+1],matches[seq(Nm+3,2*Nm-1,by=2)],0)	
			# song_pairs = matches[seq(1,2*Nm-1,by=2)] # / sum(matches[seq(1,2*Nm-1,by=2)])*sum(matches)
			offspring[,j] = song_pairs
		}		
		successful = sum(Pf_adults[,apply(pxy,3,sum)>0]) #females who dont have anyone to mate with die
		# offspring = offspring / sum(offspring)
		Pm_beforemut = apply(offspring,1,sum)/successful
		Pf_beforemut = offspring/successful 
		Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
			mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
		Pf_aftermut = Pf_beforemut ###FIX FIX FIX
		nonzero = which(Pm_aftermut>nonzero_thresh)
		perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
		# if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
			Pm[,t+1] = Pm_aftermut
			Pf[,,t+1] = Pf_aftermut
			t = t+1
			# } else{
				# Pm[,(t+1):(steps+1)] = Pm_aftermut
				# Pf[,,(t+1):(steps+1)] = Pf_aftermut
				# print(c(t,'problem'))
				# t = steps+1
				# }
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode2_pearson <-function(){
	Pm = matrix(0,Nm,steps+1) #probability of male songs over time
	Pm[,1] = m_init
	
	Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
	Pf[,,1] = P_init_pearson()
	
	reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))
	
	t = 1
	perc = 0
	while(t <= steps){
		Pm_adults = Pm[,t]
		Pf_adults = Pf[,,t]
		pxy = array(0,c(Nm,Nm,Nm)) #probability of a x / (x,y) pair
		offspring = array(0,c(Nm,Nm))
		for(j in 1:Nm){
			s = sign(midpt-j+0.5)
			weight = rep(minweight,Nm)
			toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
			pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
			weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
			z = sum(weight*Pm_adults) #normalization factor
			if(z>z_thresh){
				pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
			}
			matches = unlist(lapply(split(pxy[,,j],reverse_diagonal),sum))
			song_pairs = matches[seq(1,2*Nm-1,by=2)]+1/2*c(matches[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,matches[seq(2,2*Nm-1,by=2)])			
			# song_pairs = matches[seq(1,2*Nm-1,by=2)]+c(0,matches[seq(2,Nm-3,by=2)],matches[Nm-1]+matches[Nm+1],matches[seq(Nm+3,2*Nm-1,by=2)],0)	
			# song_pairs = matches[seq(1,2*Nm-1,by=2)] # / sum(matches[seq(1,2*Nm-1,by=2)])*sum(matches)
			offspring[,j] = song_pairs
		}		
		successful = sum(Pf_adults[,apply(pxy,3,sum)>0]) #females who dont have anyone to mate with die
		# offspring = offspring / sum(offspring)
		Pm_beforemut = apply(offspring,1,sum)/successful
		Pf_beforemut = offspring/successful 
		Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
			mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
		Pf_aftermut = Pf_beforemut ###FIX FIX FIX
		nonzero = which(Pm_aftermut>nonzero_thresh)
		perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
		# if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
			Pm[,t+1] = Pm_aftermut
			Pf[,,t+1] = Pf_aftermut
			t = t+1
			# } else{
				# Pm[,(t+1):(steps+1)] = Pm_aftermut
				# Pf[,,(t+1):(steps+1)] = Pf_aftermut
				# print(c(t,'problem'))
				# t = steps+1
				# }
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}


dynamics_mode3 <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair	
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	successful = sum(Pf_adults[apply(pxy,2,sum)>0]) #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/successful
	Pf_beforemut = Pf_adults/successful
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode4 <-function(){
Pm = array(0,c(Nm,Nm,steps+1)) #probability of male songs over time
Pm[,,1] = P_init()

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

reverse_diagonal = matrix(rep(1:Nm,each=Nm),nrow=Nm,byrow=TRUE)-(Nm+1-matrix(rep(1:Nm,times=Nm),nrow=Nm,byrow=TRUE))

t = 1
perc = 0
while(t <= steps){
	Pf_adults = Pf[,t]
	pxy = (sapply(matrix(Pf_adults,nrow=1),function(x) x*Pf_adults,simplify=TRUE))
	pref_pairs = unlist(lapply(split(pxy,reverse_diagonal),sum))
	pref_pairs = pref_pairs[seq(1,2*Nm,by=2)]	+1/2*c(pref_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,pref_pairs[seq(2,2*Nm-1,by=2)])
	pref_pairs = pref_pairs / sum(pref_pairs)	
	Pm_beforemut =  t(sapply(matrix(m_init,nrow=1),function(x) x*pref_pairs,simplify=TRUE))
	Pf_beforemut = pref_pairs
	Pm_aftermut = Pm_beforemut #and then they change their songs #FIX FIX FIX
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))		
	nonzero = which(Pf_adults>nonzero_thresh)
	perc = max(abs(range(Pf_aftermut[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode5 <-function(){
Pm = array(0,c(Nm,Nm,steps+1)) #probability of male songs over time
Pm[,,1] = P_init()

Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
Pf[,,1] = P_init()

reverse_diagonal = matrix(rep(1:Nm,each=Nm),nrow=Nm,byrow=TRUE)-(Nm+1-matrix(rep(1:Nm,times=Nm),nrow=Nm,byrow=TRUE))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,,t]
	Pf_adults = Pf[,,t]
	pxy = array(0,c(Nm,Nm,Nm)) #probability of an (x,y) / x pair
	offspring1 = array(0,c(Nm,Nm,Nm)) #probability of an (x,y) / y offspring
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*apply(Pm_adults,1,sum)) #normalization factor
		if(z>z_thresh){
			pxy[,,] = array(apply(matrix(Pf_adults[,j],nrow=1),2,function(x) x*matrix(rep(weight,times=Nm),nrow=Nm,byrow=FALSE)*Pm_adults/z),c(Nm,Nm,Nm))
			}
		for(l in 1:Nm){
			song_pairs = unlist(lapply(split(pxy[,l,],reverse_diagonal),sum))
			song_pairs = song_pairs[seq(1,2*Nm-1,by=2)]+1/2*c(song_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,song_pairs[seq(2,2*Nm-1,by=2)])
			offspring1[,l,j] = song_pairs
		}
	}	
	offspring = array(0,c(Nm,Nm)) #probability of an x y offspring
	for(j in 1:Nm){
		pref_pairs = unlist(lapply(split(offspring1[j,,],reverse_diagonal),sum))
		pref_pairs = pref_pairs[seq(1,2*Nm,by=2)]	+1/2*c(pref_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,pref_pairs[seq(2,2*Nm-1,by=2)])
		offspring[j,] = pref_pairs
	}			
	Pm_beforemut =  offspring
	Pf_beforemut = offspring
	Pm_aftermut = Pm_beforemut #and then they change their songs #FIX FIX FIX
	Pf_aftermut = Pf_beforemut #FIX FIX FIX		
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pf_aftermut[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,,t+1] = Pm_aftermut
		Pf[,,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,,(t+1):(steps+1)] = Pm_aftermut
			Pf[,,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode5_pearson <-function(){
Pm = array(0,c(Nm,Nm,steps+1)) #probability of male songs over time
Pm[,,1] = P_init_pearson()

Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
Pf[,,1] = P_init_pearson()

reverse_diagonal = matrix(rep(1:Nm,each=Nm),nrow=Nm,byrow=TRUE)-(Nm+1-matrix(rep(1:Nm,times=Nm),nrow=Nm,byrow=TRUE))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,,t]
	Pf_adults = Pf[,,t]
	pxy = array(0,c(Nm,Nm,Nm)) #probability of an (x,y) / x pair
	offspring1 = array(0,c(Nm,Nm,Nm)) #probability of an (x,y) / y offspring
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*apply(Pm_adults,1,sum)) #normalization factor
		if(z>z_thresh){
			pxy[,,] = array(apply(matrix(Pf_adults[,j],nrow=1),2,function(x) x*matrix(rep(weight,times=Nm),nrow=Nm,byrow=FALSE)*Pm_adults/z),c(Nm,Nm,Nm))
			}
		for(l in 1:Nm){
			song_pairs = unlist(lapply(split(pxy[,l,],reverse_diagonal),sum))
			song_pairs = song_pairs[seq(1,2*Nm-1,by=2)]+1/2*c(song_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,song_pairs[seq(2,2*Nm-1,by=2)])
			offspring1[,l,j] = song_pairs
		}
	}	
	offspring = array(0,c(Nm,Nm)) #probability of an x y offspring
	for(j in 1:Nm){
		pref_pairs = unlist(lapply(split(offspring1[j,,],reverse_diagonal),sum))
		pref_pairs = pref_pairs[seq(1,2*Nm,by=2)]	+1/2*c(pref_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,pref_pairs[seq(2,2*Nm-1,by=2)])
		offspring[j,] = pref_pairs
	}			
	Pm_beforemut =  offspring
	Pf_beforemut = offspring
	Pm_aftermut = Pm_beforemut #and then they change their songs #FIX FIX FIX
	Pf_aftermut = Pf_beforemut #FIX FIX FIX		
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pf_aftermut[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,,t+1] = Pm_aftermut
		Pf[,,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,,(t+1):(steps+1)] = Pm_aftermut
			Pf[,,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode6 <-function(){
Pm = array(0,c(Nm,Nm,steps+1)) #probability of male songs over time
Pm[,,1] = P_init()

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

reverse_diagonal = matrix(rep(1:Nm,each=Nm),nrow=Nm,byrow=TRUE)-(Nm+1-matrix(rep(1:Nm,times=Nm),nrow=Nm,byrow=TRUE))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,,t]
	Pf_adults = Pf[,t]
	pxy = array(0,c(Nm,Nm,Nm)) #probability of an (x,y) / y pair
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*apply(Pm_adults,1,sum)) #normalization factor
		if(z>z_thresh){
			pxy[,,j] = Pf_adults[j]*matrix(rep(weight,times=Nm),nrow=Nm,byrow=FALSE)*Pm_adults/z
			}
	}
	successful = sum(Pf_adults[apply(pxy,3,sum)>0]) #females who dont have anyone to mate with die
	offspring = array(0,c(Nm,Nm))
	for(j in 1:Nm){
		pref_pairs = unlist(lapply(split(pxy[j,,],reverse_diagonal),sum))
		pref_pairs = pref_pairs[seq(1,2*Nm,by=2)]	+1/2*c(pref_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,pref_pairs[seq(2,2*Nm-1,by=2)])
		offspring[j,] = pref_pairs
	}			
	Pm_beforemut =  offspring
	Pf_beforemut = apply(offspring,2,sum)
	Pm_aftermut = Pm_beforemut #and then they change their songs #FIX FIX FIX
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))		
	nonzero = which(Pf_adults>nonzero_thresh)
	perc = max(abs(range(Pf_aftermut[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode7 <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	successful = sum(Pf_adults[apply(pxy,2,sum)>0]) #females who dont have anyone to mate with die
	Pm_beforemut = Pm_adults
	Pf_beforemut = apply(pxy,1,sum)/successful
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))		
	nonzero = which(Pf_adults>nonzero_thresh)
	perc = max(abs(range(Pf_aftermut[nonzero]/Pf_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode8 <-function(){
	Pm = matrix(0,Nm,steps+1) #probability of male songs over time
	Pm[,1] = m_init
	
	Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
	Pf[,,1] = P_init()
	
	reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))
	
	t = 1
	perc = 0
	while(t <= steps){
		Pm_adults = Pm[,t]
		Pf_adults = Pf[,,t]
		pxy = array(0,c(Nm,Nm,Nm)) #probability of a x / (x,y) pair
		offspring = array(0,c(Nm,Nm,Nm)) #probability of an offspring with an i father, who sings j, with a mother who prefers k 
		for(j in 1:Nm){
			s = sign(midpt-j+0.5)
			weight = rep(minweight,Nm)
			toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
			pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
			weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
			z = sum(weight*Pm_adults) #normalization factor
			if(z>z_thresh){
				pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
			}
			song_pairs = t(apply(cbind(c(0,rep(1:(floor(Nm)/2),each=2)),c(0,rep(c(1,0),times=(floor(Nm)/2))),pxy[,,j]),1,function(v) c(rep(0,v[1]),v[seq(3+v[2],Nm+2,by=2)]+1/2*c(2*v[2]*v[3],v[seq(4+v[2],Nm+2-v[2],by=2)])+1/2*c(v[seq(4+v[2],Nm+2-2*v[2],by=2)],2*v[2]*v[Nm+2]),rep(0,v[2]+(floor(Nm)/2)-v[1]))))
			offspring[,,j] = song_pairs
		}		
		successful= sum(Pf_adults[,apply(pxy,3,sum)>0]) #females who dont have anyone to mate with die
		Pm_beforemut = apply(offspring,2,sum)/successful
		Pf_beforemut = t(apply(offspring,c(1,2),sum))/successful 	
		Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
			mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
		Pf_aftermut = Pf_beforemut ###FIX FIX FIX
		nonzero = which(Pm_aftermut>nonzero_thresh)
		perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
		if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
			Pm[,t+1] = Pm_aftermut
			Pf[,,t+1] = Pf_aftermut
			t = t+1
			} else{
				Pm[,(t+1):(steps+1)] = Pm_aftermut
				Pf[,,(t+1):(steps+1)] = Pf_aftermut
				print(c(t,'problem'))
				t = steps+1
				}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode9 <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair	
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	successful = sum(Pf_adults[apply(pxy,2,sum)>0]) #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/successful
	Pf_beforemut = apply(pxy,1,sum)/successful
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode11 <-function(){
	Pm = matrix(0,Nm,steps+1) #probability of male songs over time
	Pm[,1] = m_init
	
	Pf = array(0,c(Nm,Nm,steps+1)) #probability of female preferences over time
	Pf[,,1] = P_init()
	
	reverse_diagonal = row(Pf[,,1])-(Nm+1-col(Pf[,,1]))
	
	t = 1
	perc = 0
	while(t <= steps){
		Pm_adults = Pm[,t]
		Pf_adults = Pf[,,t]
		pxy = array(0,c(Nm,Nm,Nm)) #probability of a x / (x,y) pair
		offspring = array(0,c(Nm,Nm))
		for(j in 1:Nm){
			s = sign(midpt-j+0.5)
			weight = rep(minweight,Nm)
			toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
			pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
			weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
			z = sum(weight*Pm_adults) #normalization factor
			if(z>z_thresh){
				pxy[,,j] = outer(weight*Pm_adults/z,Pf_adults[,j])
			}
			song_pairs = unlist(lapply(split(pxy[,,j],reverse_diagonal),sum))
			song_pairs = song_pairs[seq(1,2*Nm-1,by=2)]+1/2*c(song_pairs[seq(2,2*Nm-1,by=2)],0)+1/2*c(0,song_pairs[seq(2,2*Nm-1,by=2)])
			offspring[,j] = song_pairs
		}		
		successful = sum(Pf_adults[,apply(pxy,3,sum)>0]) #females who dont have anyone to mate with die
		Pm_beforemut = apply(offspring,1,sum)/successful
		Pf_beforemut = sapply(matrix(Pm_adults,nrow=1),function(x) x*apply(offspring,1,sum),simplify=TRUE)/successful 
		Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
			mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
		Pf_aftermut = Pf_beforemut ###FIX FIX FIX
		nonzero = which(Pm_aftermut>nonzero_thresh)
		perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
		if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
			Pm[,t+1] = Pm_aftermut
			Pf[,,t+1] = Pf_aftermut
			t = t+1
			} else{
				Pm[,(t+1):(steps+1)] = Pm_aftermut
				Pf[,,(t+1):(steps+1)] = Pf_aftermut
				print(c(t,'problem'))
				t = steps+1
				}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,,(store):(steps+1)])
return(pop_dens)
}

dynamics_mode12 <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair	
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	successful = sum(Pf_adults[apply(pxy,2,sum)>0]) #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/successful
	Pf_beforemut = Pm_adults/successful
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	Pf_aftermut = ((1-mut_prob)*Pf_beforemut + mut_prob/2*c(Pf_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pf_beforemut[1:(Nm-1)]))
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_aftermut
		t = t+1
		} else{
			Pm[,(t+1):(steps+1)] = Pm_aftermut
			Pf[,(t+1):(steps+1)] = Pf_aftermut
			t = steps+1
			}
}
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)])
return(pop_dens)
}

