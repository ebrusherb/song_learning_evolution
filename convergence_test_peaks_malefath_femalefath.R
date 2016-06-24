source('ind2sub.R')
source('glue.R')
source('int.R')
source('range_setup.R')


## ---- dynamics -----------------------

dynamics_malefath_femalefath <-function(){
Pm = matrix(0,Nm,Tsteps+1) #probability of male songs over time
Pm[,1] = m_init

Pf = matrix(0,Nf,Tsteps+1) #probability of female preferences over time
Pf[,1] = f_init

t = 1
perc = 0
while(t <= Tsteps){
	Pm_adults = Pm[,t]
	Pf_adults = Pf[,t]
	pxy = matrix(0,Nm,Nf) #probability of a (x,y) pair
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)
	fixed_weight[fixed_weight==0] = minweight
	fixed_weight = fixed_weight/sum(fixed_weight)
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm_adults) #normalization factor
		if(z>z_thresh){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	Pf_adults[apply(pxy,2,sum)==0]=0
	song_beforemut = apply(pxy,1,sum)/sum(Pf_adults)	
	song_aftermut = (1-mut_prob)*song_beforemut + mut_prob/2*c(song_beforemut[2:Nm],0) + 
		mut_prob/2*c(0,song_beforemut[1:Nm-1]) #and then they change their songs
	song_aftermut[which(song_aftermut<=killthresh)] = 0
	song_beforemut[which(song_beforemut<=killthresh)] = 0
	nonzero = which(Pm_adults>nonzero_thresh)
	perc = max(abs(range(song_aftermut[nonzero]/Pm_adults[nonzero],na.rm=TRUE)-c(1,1)))
	if(perc>perc_thresh){	
		Pm[,t+1] = song_aftermut
		Pf[,t+1] = song_beforemut
		t = t+1
		} else{
			Pm[,(t+1):(Tsteps+1)] = song_aftermut
			Pf[,(t+1):(Tsteps+1)] = song_beforemut
			t = Tsteps+1
			}
}
pop_dens = list(Pm=Pm[,(store):Tsteps],Pf=Pf[,(store):Tsteps])
return(pop_dens)
}

# d = dim(Pm_onepop)
# P = prod(d)

# Ns = length(sigma_vals)
# Nfs = length(f_sigma_vals)
# Nms = length(m_sigma_vals)
# Nmp = length(mut_prob_vals)

# subset = 601:1201;
# subset = 1:Nm;
# Tend = dim(Pm_onepop[[1]])[2]

# thresh = 1e-2
# half = floor(Nm/2)

# four_freq = array(NA,dim(Pm_onepop))
# four_mat = list()
# num_peaks = array(NA,dim(Pm_onepop))
# var_mat_m = array(NA,c(dim(Pm_onepop)))
# var_mat_f = array(NA,c(dim(Pf_onepop)))
# ent_mat = array(Inf,dim=dim(Pm_onepop))



# for(i in 1:P){
	# eq_pop = Pm_onepop[[i]][,Tend]
	# fourComps = fft(eq_pop[subset]);
	# fourCoeffs = abs(fourComps);
	# coeffs = fourCoeffs / (length(subset)/2);
	# coeffs[1] = fourCoeffs[1] / length(subset)	
	# mf = mean(coeffs[1:200])
	# ml = intersect(local_max(round(coeffs[1:half],5)),which(coeffs[1:half]>mf/2))
	# if(length(ml)==1 || length(ml)==0){four_freq[i] = 1} else{ four_freq[i] = median(diff(ml))}
	# four_mat[[i]] = coeffs[1:half]
	# m = local_max(eq_pop)
	# m = intersect(m, which(eq_pop>thresh))
	# l = length(m)
	# num_peaks[i] = l
	# sub = ind2sub(dim(Pm_onepop),i)
	# ex = int(mrange*Pm_onepop[[i]][,Tend])
	# vx = int((mrange-(ex))^2*Pm_onepop[[i]][,Tend])
	# var_mat_m[i] = vx
	# ex = int(mrange*Pf_onepop[[i]][,Tend])
	# vx = int((mrange-(ex))^2*Pf_onepop[[i]][,Tend])
	# var_mat_f[i] = vx
	# ent_mat[i] = ent(eq_pop)
# }
# dim(four_mat) = dim(Pm_onepop)

# smat= array(0,dim=dim(Pm_onepop))
# fmat= array(0,dim=dim(Pm_onepop))
# mmat= array(0,dim=dim(Pm_onepop))
# pmat = array(0,dim=dim(Pm_onepop))

# for(i in 1:P){
	# sub = ind2sub(dim(Pm_onepop),i)
	# smat[i] = sigma_vals[sub[1]]
	# fmat[i] = f_sigma_vals[sub[2]]
	# mmat[i] = m_sigma_vals[sub[3]]
	# pmat[i] = mut_prob_vals[sub[4]]
# }


##

# w = intersect(which(num_peaks >1),(intersect(which(fmat>0.001),which(mmat>0.001))))

# subs = c()

# Tsteps = 500

# Pm_test = list()
# Pf_test = list()

# for(i in 1:length(w)){
	# subs = ind2sub(dim(Pm_onepop),w[i])
	# sigma = sigma_vals[subs[1]]
	# f_sigma = f_sigma_vals[subs[2]]
	# m_sigma = m_sigma_vals[subs[3]]
	# mut_prob = mut_prob_vals[subs[4]]
	# f_init = dnorm(frange,fmin,f_sigma)
	# m_init = dnorm(mrange,mmin,m_sigma)
	# pop_dens = dynamics_full()
	# Pm_test[[i]] = pop_dens$Pm
	# Pf_test[[i]] = pop_dens$Pf
# }

# i = 2
# subs = ind2sub(dim(Pm_onepop),w[i])
# sigma = sigma_vals[subs[1]]
# f_sigma = f_sigma_vals[subs[2]]
# m_sigma = m_sigma_vals[subs[3]]
# Pm = Pm_test[[i]]
# Pf = Pf_test[[i]]

killthresh=0
Tsteps = 25
store = 1
sigma = 0.1
f_sigma = 1
m_sigma = 0.1
mut_prob = 0.01
pf = 1
pm = 1
f_init = dnorm(frange,fmin,f_sigma)
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,m_sigma)
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

pop_dens = dynamics_malefath_femalefath()
Pm = pop_dens$Pm
Pf = pop_dens$Pf

Tmin = Tsteps

preferences = array(0,dim=c(Nm,Nf,Tmin))
female_tots = matrix(0,Nf,Tmin)
pxy = array(0,dim=c(Nm,Nf,Tmin))
growth = array(0,dim=c(Nm,Tmin))
growth_desperate = array(0,dim=c(Nm,Tmin))

for(t in 1:Tmin){
	fixed_weight = dnorm(mrange,mean=frange[midpt],sd=sigma) #female preference function
	minweight = 10^max(floor(log(min(fixed_weight[which(fixed_weight>0)]),base=10)),-320)	
	fixed_weight[fixed_weight==0] = minweight
	fixed_weight = fixed_weight/sum(fixed_weight)
	for(j in 1:Nf){
		s = sign(midpt-j+0.5)
		weight = rep(minweight,Nf)
		toreplace=sort(c(((Nf+1)+s)%%(Nf+1),j+s*(midpt-1)))
		pull=sort(c(((Nf+1)+s)%%(Nf+1)+midpt-j,((Nf+1)-s)%%(Nf+1)))
		weight[toreplace[1]:toreplace[2]] = fixed_weight[pull[1]:pull[2]]
		z = sum(weight*Pm[,t])
		preferences[,j,t] = weight #preference given by each female to each male
		
		female_tots[j,t] = z #total preferences given by each female
		if(z>z_thresh){
			pxy[,j,t] = Pf[j,t]*weight*Pm[,t]/z
			}
	}
	zero = which(Pm[,t]==0)
	growth[,t] =  apply(pxy[,,t],1,sum)/Pm[,t]
	growth[zero,t] = 0
	growth_desperate[,t] = apply(preferences[,,t]/matrix(rep(female_tots[,t],times=Nm),nrow=Nm,byrow=TRUE),1,sum)
	growth_desperate[zero,t] = 0
}

# t=1
# plot(female_tots[,t],t='l',xlim=c(450,850),ylim=c(0,2))
# lines(Pm[,t],col='red')
# lines(Pf[,t],col='green')

t=2;plot(Pm[,t]);vars = sigma_malefath_femalefath(m_sigma^2,f_sigma^2);v=dnorm(mrange,mmin,sqrt(vars[t,1]));lines(v/sum(v),col='red')
