dynamics_pxy <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

pxy_store = array(0,c(Nm,Nm,steps+1))
z_store = array(0,c(Nm,steps+1))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		z_store[j,t] = z
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	pxy_store[,,t] = pxy
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
			z_store[,t:steps] = z_store[,t]
			pxy_store[,,t:steps] = pxy_store[,,t]
			t = steps+1
			}
}
t=steps+1
Pm_adults = Pm[,t]
Pf_adults = Pf[,t]
pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
for(j in 1:Nm){
	s = sign(midpt-j+0.5)
	weight = rep(minweight,Nm)
	toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
	pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
	weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
	z = sum(weight*Pm_adults) #normalization factor
	z_store[j,t] = z
	if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
}
pxy_store[,,t] = pxy
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)],pxy=pxy_store[,,(store):(steps+1)],z=z_store[,(store):(steps+1)])
return(pop_dens)
}

dynamics_bothmut <-function(){
Pm = matrix(0,Nm,steps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nm,steps+1) #probability of female preferences over time
Pf[,1] = f_init

pxy_store = array(0,c(Nm,Nm,steps+1))
z_store = array(0,c(Nm,steps+1))

t = 1
perc = 0
while(t <= steps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nm){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nm)
		toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
		pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		z_store[j,t] = z
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	pxy_store[,,t] = pxy
	Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = ((1-mut_prob)*Pf_adults + mut_prob/2*c(Pf_adults[2:Nm],0) + 
		mut_prob/2*c(0,Pf_adults[1:(Nm-1)]))/sum(Pf_adults)
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
			z_store[,t:steps] = z_store[,t]
			pxy[,,t:steps] = pxy[,,t]
			t = steps+1
			}
}
t=steps+1
Pm_adults = Pm[,t]
Pf_adults = Pf[,t]
pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
for(j in 1:Nm){
	s = sign(midpt-j+0.5)
	weight = rep(minweight,Nm)
	toreplace=sort(c((s)%%(Nm+1),j+s*(midpt-1)))
	pull=sort(c((s)%%(Nm+1)+midpt-j,(-s)%%(Nm+1)))
	weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
	z = sum(weight*Pm_adults) #normalization factor
	z_store[j,t] = z
}
pxy_store[,,t] = pxy
pop_dens = list(Pm=Pm[,(store):(steps+1)],Pf=Pf[,(store):(steps+1)],pxy=pxy_store[,,(store):(steps+1)],z=z_store[,(store):(steps+1)])
return(pop_dens)
}

dynamics_memory <-function(store=FALSE){
Pm = matrix(0,Nm) #probability of male songs over time
Pm = m_init

Pf = matrix(0,Nm) #probability of female preferences over time
Pf = f_init

Pm_store = c()
Pf_store=c()

t = 1

if(store && is.element(t,store_vec)){
			Pm_store = cbind(Pm_store,m_init)
			Pf_store = cbind(Pf_store,f_init)
		}

perc = 0
while(t <= steps){
	Pm_adults = Pm
	Pf_adults = Pf
	pxy = matrix(0,Nm,Nm) #probability of a (x,y) pair
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
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
	Pf_adults[apply(pxy,2,sum)==0]=0 #females who dont have anyone to mate with die
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = ((1-mut_prob)*Pf_adults + mut_prob/2*c(Pf_adults[2:Nm],0) + 
		mut_prob/2*c(0,Pf_adults[1:(Nm-1)]))/sum(Pf_adults)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:(Nm-1)]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	#if growth rate is basically 1 everywhere stop simulations to save time
		Pm = Pm_aftermut
		Pf = Pf_adults		
		t = t+1
		if(store && is.element(t,store_vec)){
			Pm_store = cbind(Pm_store,Pm)
			Pf_store = cbind(Pf_store,Pf)
		}
		} else{
			Pm = Pm_aftermut
			Pf = Pf_adults
			t = steps+1
			}
}
if(!store){
	pop_dens = list(Pm=Pm,Pf=Pf)
	} else{ pop_dens =list(Pm=Pm,Pm_store=Pm_store,Pf=Pf,Pf_store=Pf_store)}
return(pop_dens)
}