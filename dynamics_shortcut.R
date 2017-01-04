

dynamics_shortcut <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

midpt = ceiling(Nf/2)

t = 1
perc = 0
while(t <= Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	fixed_weight = fixed_weight/sum(fixed_weight)
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)
	fixed_weight[fixed_weight==0] = minweight
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		# if(j < ceiling(Nf/2)){
			# weight[1:(Nf-ceiling(Nf/2)+j)] = fixed_weight[(ceiling(Nf/2)-j+1):Nf]
		# }
		# if(j == ceiling(Nf/2)){
			# weight = fixed_weight
		# }
		# if(j > ceiling(Nf/2)){
			# weight[(j-ceiling(Nf/2)+1):Nf] = fixed_weight[1:(Nf+ceiling(Nf/2)-j)]
		# }
		z = sum(weight*Pm_adults) #normalization factor
		if(z>1e-20){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pf_adults[apply(pxy,2,sum)==0]=0
	Pm_beforemut = apply(pxy,1,sum)/sum(Pf_adults)
	Pf_adults = Pf_adults/sum(Pf_adults)
	Pm_aftermut = matrix(0,Nm)
	Pm_aftermut = (1-mut_prob)*Pm_beforemut + mut_prob/2*c(Pm_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,Pm_beforemut[1:Nm-1]) #and then they change their songs
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(Pm_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,t+1] = Pm_aftermut
		Pf[,t+1] = Pf_adults
		t = t+1
		} else{
			Pm[,(t+1):(Tsteps+1)] = Pm_aftermut
			Pf[,(t+1):(Tsteps+1)] = Pf_adults
			t = Tsteps+1
			}
}
pop_dens = list(Pm=Pm[,(1):Tsteps],Pf=Pf[,(1):Tsteps])
return(pop_dens)
}