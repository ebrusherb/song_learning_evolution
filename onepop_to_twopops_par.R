## --- cluster set up
library(parallel)
library(foreach)
library(doParallel)

# num_cores <- detectCores()-1
num_cores <-20
cl <-makeCluster(num_cores)
registerDoParallel(cl)

source('ind2sub.R')
source('glue.R')
source('int.R')


direc = "/homes/ebrush/priv/song_learning_evolution/"
files = list.files(path=direc)
datafiles = intersect(grep('song_learning_paramsweep_par_2016',files,value=FALSE),grep('Rdata',files,value=FALSE))
totals = matrix(NA,length(datafiles))
for(i in 1:length(datafiles)){
	file = files[datafiles[i]]
	file = substr(file,1,nchar(file)-6)
	day = strtoi(substr(file,nchar(file)-1,nchar(file)))
	month = strtoi(substr(file,nchar(file)-4,nchar(file)-3))
	year = strtoi(substr(file,nchar(file)-7,nchar(file)-6))
	total = day + 31*month +365*year
	totals[i] = total
}
file = paste(direc,files[datafiles[which.max(totals)]],sep="") #find most recent paramsweep file

load(file)


## ---- parameters ----------------
step = 0.01 #step size of trait space
int_step = step #step to use for integration function
# alpha = 0.5 #if preference function is a step fx, strength of preference
# sigma2 = #variance of female preference function
# mut_prob =  #probability a male changes song to one on either side
# mut_delta = #how to implement mutations of different sizes?
# fmix_sigma2 = #variance of female distribution(s)
# mmix_sigma2 = 0.1 #variance of male distribution(s)
# Tsteps = #how many generations

mrange = seq(-7.5,7.5,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-7.5,7.5,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)

nonzero_thresh = 1e-5
perc_thresh = 1e-10

## ---- dynamics -----------------------

dynamics <-function(){
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
	### should I round Pm_adults?!? how?!? is that why I'm getting bumps?!?
	for(j in 1:Nf){
		y = frange[j]
		# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
		weight = dnorm(mrange,mean=y,sd=sqrt(sigma2)) #female preference function
		# weight = matrix (0,Nf,1)
		# weight[c(f0,x1)] = 1
		# weight[j] = 1+alpha
		z = int(weight*Pm_adults) #normalization factor
		if(z!=0){
			pxy[,j] = Pf_adults[j]*weight*Pm_adults/z
			}
	}
	# Pm_beforemut = matrix(0,Nm)
	# for(i in 1:Nm){
		# Pm_beforemut[i] = int(pxy[i,]) #probability of males being born
	# }
	Pm_beforemut = apply(pxy,1,int)
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
pop_dens = list(Pm=Pm[,c(1,(Tsteps-200):Tsteps)],Pf=Pf[,c(1,(Tsteps-200):Tsteps)])
return(pop_dens)
}

## ---- variance_sweep ---------------
Tsteps = 20000
Tend = dim(Pm_onepop[[1]])[2]
pm = 0.6
pf = 0.4

smallf = which(fmix_sigma2_vals<0.01)
fvals = (length(smallf) + 1):(length(fmix_sigma2_vals))
if(length(smallf)>0){
	fmix_sigma2_vals = fmix_sigma2_vals[-smallf]
}

smallp = which(mut_prob_vals == 0)
pvals = (length(smallp)+1):(length(mut_prob_vals))
if(length(smallp)>0){
	mut_prob_vals = mut_prob_vals[-smallp]
}


Ns = length(sigma2_vals)
Nfs = length(fmix_sigma2_vals)
Nms = length(mmix_sigma2_vals)
Nmp = length(mut_prob_vals)

shift_vals = c(0,37)
Nshift = length(shift_vals)
P = Ns*Nfs*Nms*Nmp*Nshift
d = c(Ns,Nfs,Nms,Nmp,Nshift)


Pm_twopop=as.list(1:P)
dim(Pm_twopop)<-d
Pf_twopop=as.list(1:P)
dim(Pf_twopop)<-d

P_twopop<-foreach(ind = 1:P, .combine='glue', .multicombine = TRUE, .init=list(list(),list())) %dopar% {
	v=ind2sub(d,ind)
	s=v[1]
	f=fvals[v[2]]
	m=v[3]
	p=pvals[v[4]]
	shift = shift_vals[v[5]]
	morig = Pm_onepop[[s,f,m,p]][,Tend]
	forig= Pf_onepop[[s,f,m,p]][,Tend]
	
	mpeak_diff = m1-m0+shift
	fpeak_diff = f1-f0+shift

	mshift = matrix(0,Nm)
	mshift[mpeak_diff:Nm] = morig[1:(Nm-mpeak_diff+1)]
	fshift = matrix(0,Nf)
	fshift[fpeak_diff:Nf] = forig[1:(Nf-fpeak_diff+1)]

	m_init = pm*morig+(1-pm)*mshift
	f_init = pf*forig + (1-pf)*fshift

	sigma2 = sigma2_vals[s]	
	mut_prob = mut_prob_vals[v[4]]
	
	pop_dens = dynamics()
	# pop_dens_last = list(pop_dens$Pm[,Tsteps],pop_dens$Pf[,Tsteps])
}

Pm_twopop_hold = P_twopop[[1]]
Pf_twopop_hold = P_twopop[[2]]

for(ind in 1:P){
	Pm_twopop[[ind]] = Pm_twopop_hold[[ind]]
	Pf_twopop[[ind]] = Pf_twopop_hold[[ind]]
}

Date=Sys.Date()
# save(Pm_keep=Pm_keep,Pf_keep=Pf_keep,Pm_onepop=Pm_onepop,Pf_onopop=Pf_onepop,sigma2_vals=sigma2_vals,fmix_sigma2_vals,mmix_sigma2_vals,mut_prob_vals=mut_prob_vals,file='/homes/ebrush/priv/song_learning_evolution/song_learning_paramsweep_par.Rdata')
save(Pm_twopop=Pm_twopop,Pf_twopop=Pf_twopop,sigma2_vals=sigma2_vals,fmix_sigma2_vals,mmix_sigma2_vals,mut_prob_vals=mut_prob_vals,Tsteps=Tsteps,mrange,Nm,frange,Nf,step,int_step,shift_vals,file=paste('/homes/ebrush/priv/song_learning_evolution/onepop_to_twopops_par_',substr(Date,1,4),'_',substr(Date,6,7),'_',substr(Date,9,10),'.Rdata',sep=''))


stopCluster(cl)

quit()




