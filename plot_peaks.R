library(RColorBrewer)
mypal=brewer.pal(5,'Set1')

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

mrange = seq(-10,10,by=step) #range of male songs
Nm = length(mrange) 
mmin = -1
mmax = 1
m0 = which(mrange==mmin)
m1 = which(mrange==mmax)
mrange_orig = seq(mmin,mmax,by=step) 
frange = seq(-10,10,by=step) #range of female preferences
Nf = length(frange)
fmin = -1
fmax = 1
f0 = which(frange==fmin)
f1 = which(frange==fmax)
frange_orig = seq(fmin,fmax,by=step)

P=prod(dim(Pm_keep))
d=dim(Pm_keep[[1]])
Tsteps=d[2]
## --- 

par(mfrow=c(5,5))
for(k in 0:4){
start = 1+25*k
for(ind in start:(start+25-1)){
	x=ind2sub(dim(Pm_keep),ind)
	s=x[1];f=x[2];m=x[3]
	v =Pm_keep[[ind]][,Tsteps]
	plot(v/max(v),type='l',xlab='',ylab='',main=c(sigma2_vals[s],fmix_sigma2_vals[f],mmix_sigma2_vals[m]))
	v =Pm_keep[[ind]][,1]
	lines(v/max(v),col='red')
	v = Pf_keep[[ind]][,1]
	lines(v/max(v),col='green')
}
}

par(mfrow=c(3,4))

Tsteps=d[2]
m=1
for(s in 1:3){
	for(f in 2:5){
	v =Pm_keep[[s,f,m]][,1]
	plot(mrange,v,type='l',xlab='Song',ylab='Density',main=paste("Pref width =",sigma2_vals[s],'Fem var =',fmix_sigma2_vals[f],"Male var =",mmix_sigma2_vals[m]),ylim=c(0,1.5),col='red')
	v =Pm_keep[[s,f,m]][,Tsteps]
	lines(mrange,v,col='black')
	v = Pf_keep[[s,f,m]][,1]
	lines(mrange,v,col='blue')
	}
}


par(mfrow=c(2,5))

Tsteps=d[2]
f = 5
m = 1
for(s in 1:5){
	v =Pm_keep[[s,f,m]][,Tsteps]
	plot(mrange,v,type='l',ylim=c(0,1.5),xlab='Song',ylab='Density',main=paste("Pref width = ",sigma2_vals[s],"Male var = ",mmix_sigma2_vals[m]))
	v =Pm_keep[[s,f,m]][,1]
	lines(mrange,v,col='red')
	v = Pf_keep[[s,f,m]][,1]
	lines(mrange,v,col='blue')
}
m = 2
for(s in 1:5){
	v =Pm_keep[[s,f,m]][,Tsteps]
	plot(mrange,v,type='l',ylim=c(0,1.5),xlab='Song',ylab='Density',main=paste("Pref width = ",sigma2_vals[s],"Male var = ",mmix_sigma2_vals[m]))
	v =Pm_keep[[s,f,m]][,1]
	lines(mrange,v,col='red')
	v = Pf_keep[[s,f,m]][,1]
	lines(mrange,v,col='blue')
}

par(mfrow=c(5,5))
for(k in 0:4){
P=prod(dim(Pm_onepop))
d=dim(Pm_onepop[[1]])
Tsteps=d[2]
start = 1+25*k
for(ind in start:(start+25-1)){
	x=ind2sub(dim(Pm_onepop),ind)
	s=x[1];f=x[2];m=x[3]
	v =Pm_onepop[[ind]][,Tsteps]
	plot(v/max(v),type='l',xlab='',ylab='',main=c(sigma2_vals[s],fmix_sigma2_vals[f],mmix_sigma2_vals[m]))
	v =Pm_onepop[[ind]][,1]
	lines(v/max(v),col='red')
	v = Pf_onepop[[ind]][,1]
	lines(v/max(v),col='green')
}
}

par(mfrow=c(2,3))
m=5
for(f in 1:5){
	s=1
	plot(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_onepop[[s,f,m]][,Tsteps],type='l',col=mypal[s],ylim=c(0,0.5),xlim=c(-1.5,1.5),main=fmix_sigma2_vals[f])
	for(s in 2:5){
		lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_onepop[[s,f,m]][,Tsteps],col=mypal[s])
	}
	lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pf_onepop[[s,f,m]][,Tsteps],col='black')
}
legend(1,.4,legend=sigma2_vals,col=mypal,lty=1,bty='n')
s=5
f=1
plot(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_onepop[[s,f,m]][,Tsteps],type='l',col=mypal[f],ylim=c(0,0.5),xlim=c(-1.5,1.5))
	for(f in 2:5){
		lines(sign(mrange-mrange[m0])*sqrt(abs(mrange-mrange[m0])),Pm_onepop[[s,f,m]][,Tsteps],col=mypal[f])
	}
legend(1,.4,legend=fmix_sigma2_vals,col=mypal,lty=1,bty='n')

par(mfrow=c(3,4))

Tsteps=d[2]
m=1
for(s in 1:3){
	for(f in 2:5){
	v =Pm_onepop[[s,f,m]][,1]
	plot(mrange,v,type='l',xlab='Song',ylab='Density',main=paste("Pref width =",sigma2_vals[s],'Fem var =',fmix_sigma2_vals[f],"Male var =",mmix_sigma2_vals[m]),ylim=c(0,1.5),col='red')
	v =Pm_onepop[[s,f,m]][,Tsteps]
	lines(mrange,v,col='black')
	v = Pf_onepop[[s,f,m]][,1]
	lines(mrange,v,col='blue')
	}
}

par(mfrow=c(2,5))

Tsteps=d[2]
f = 5
m = 1
for(s in 1:5){
	v =Pm_onepop[[s,f,m]][,Tsteps]
	plot(mrange,v,type='l',ylim=c(0,1.5),xlab='Song',ylab='Density',main=paste("Pref width = ",sigma2_vals[s],"Male var = ",mmix_sigma2_vals[m]))
	v =Pm_onepop[[s,f,m]][,1]
	lines(mrange,v,col='red')
	v = Pf_onepop[[s,f,m]][,1]
	lines(mrange,v,col='blue')
}
m = 2
for(s in 1:5){
	v =Pm_onepop[[s,f,m]][,Tsteps]
	plot(mrange,v,type='l',ylim=c(0,1.5),xlab='Song',ylab='Density',main=paste("Pref width = ",sigma2_vals[s],"Male var = ",mmix_sigma2_vals[m]))
	v =Pm_onepop[[s,f,m]][,1]
	lines(mrange,v,col='red')
	v = Pf_onepop[[s,f,m]][,1]
	lines(mrange,v,col='blue')
}

## ----
num_peaks = array(NA,dim(Pm_keep))
wavelength_peaks = array(Inf,dim(Pm_keep))
ent_mat = array(Inf,dim=dim(Pm_keep))

thresh = 1e-4

for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tsteps]
	m = local_max(eq_pop)
	m = intersect(m, which(eq_pop>thresh))
	l = length(m)
	num_peaks[i] = l
	if(l > 1){
		wavelength_peaks[i] = step*median(diff(m))
	}
	ent_mat[i] = ent(eq_pop)
}

subset = 601:1201;
subset = 1:Nm;
half = floor(length(subset)/2)
four_freq = array(NA,dim(Pm_keep))
# four_coef = array(NA,dim(Pm_keep))
four_mat = array(NA,c(dim(Pm_keep),half))
var_mat = array(NA,c(dim(Pm_keep)))
four_var = array(NA,c(dim(Pm_keep)))

for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tsteps]
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourierCoefficients[1] / length(subset)
	# v=order(coeffs[1:floor(length(subset)/2)],decreasing=TRUE)
	# # if(v[1]==2){four_freq[i]=v[2]} else{ four_freq[i]=v[1]}
	# four_coef[i] = coeffs[2]
	sub = ind2sub(dim(Pm_keep),i)
	four_mat[cbind(matrix(rep(sub,half),nrow=half,byrow=TRUE),1:half)] = coeffs[1:half]
	ex = int(mrange*Pm_onepop[[i]][,Tsteps])
	vx = int((mrange-(ex))^2*Pm_onepop[[i]][,Tsteps])
	var_mat[i] = vx
	mf = mean(coeffs[1:200])
	ml = intersect(local_max(coeffs[1:half]),which(coeffs[1:half]>mf/2))
	if(length(ml)==1){four_freq[i] = 1} else{ four_freq[i] = median(diff(ml))}
	l = lm(log(coeffs[2:7]) ~ poly(1:6,2,raw=TRUE));
	four_var[i] = -summary(l)$coefficients[3,1]
}

layout(matrix(1:10,nrow=2,byrow=FALSE))
f=5
m=2
for(s in 1:5){
	plot(Pm_onepop[[s,f,m]][,Tsteps],type='l')
	lines(Pf_onepop[[s,f,m]][,Tsteps],col='red')
	plot(four_mat[s,f,m,2:7],type='l')
}


layout(matrix(1:10,nrow=2,byrow=FALSE))
f=2;m=1;
for(s in 1:5){
	eq_pop = Pm_onepop[[s,f,m]][,Tsteps]
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
		coeffs[1] = fourierCoefficients[1] / length(subset)
 	v =Pm_onepop[[s,f,m]][,Tsteps]
 	plot(v/max(v),type='l',xlab='',ylab='',main=c(sigma2_vals[s],fmix_sigma2_vals[f],mmix_sigma2_vals[m]))
 	v =Pm_onepop[[s,f,m]][,1]
 	lines(v/max(v),col='red')
 	v = Pf_onepop[[ind]][,1]
 	lines(v/max(v),col='green')
 	plot(four_mat[s,f,m,2:200],type='l')
	mf = mean(coeffs[1:200])
 	ml = intersect(local_max(coeffs[1:half]),which(coeffs[1:half]>mf/2))
 	abline(v=ml,col='blue')
 }

smat= array(0,dim=dim(Pm_keep))
fmat= array(0,dim=dim(Pm_keep))
mmat= array(0,dim=dim(Pm_keep))
colmat = array(0,dim=dim(Pm_keep))

for(i in 1:P){
	sub = ind2sub(dim(Pm_keep),i)
	smat[i] = sigma2_vals[sub[1]]
	fmat[i] = fmix_sigma2_vals[sub[2]]
	mmat[i] = mmix_sigma2_vals[sub[3]]
	colmat[i] = mypal[sub[2]]
}


m=2;symbols(rev(as.vector(smat[,,m])),rev(as.vector(num_peaks[,,m])),circles=(rev(as.vector(ent_mat[,,m]-min(ent_mat)+.01)))^2,bg=rev(as.vector(colmat[,,m])),inches=.35)

m=1;radius=sqrt(rev(as.vector(1/wavelength_peaks[,,m]+.05)));symbols(rev(as.vector(smat[,,m])),rev(as.vector(ent_mat[,,m])),circles=radius,bg=rev(as.vector(colmat[,,m])),inches=.1*max(radius),xlab='Selectivity',ylab='Entropy');legend(x=1.5,y=01,legend=(fmix_sigma2_vals),lty=0,pch=16,col=(mypal),bty='n')

smat= array(0,dim=dim(Pm_keep))
fmat= array(0,dim=dim(Pm_keep))
mmat= array(0,dim=dim(Pm_keep))
colmat = array(0,dim=dim(Pm_keep))

for(i in 1:P){
	sub = ind2sub(dim(Pm_keep),i)
	smat[i] = sigma2_vals[sub[1]]
	fmat[i] = fmix_sigma2_vals[sub[2]]
	mmat[i] = mmix_sigma2_vals[sub[3]]
	colmat[i] = mypal[sub[3]]
}

s=1;radius=sqrt((as.vector(1/wavelength_peaks[s,,])));symbols((as.vector(fmat[s,,])),(as.vector(ent_mat[s,,])),circles=radius,bg=(as.vector(colmat[s,,])),inches=.1*max(radius));legend(x=1.5,y=01,legend=(fmix_sigma2_vals),lty=0,pch=16,col=(mypal),bty='n')


t_peaks_mat = array(0,dim(Pm_keep))

for(i in 1:P){
	t = 1
	where_peaks = local_max(Pm_onepop[[i]][,t])
	where_peaks = intersect(where_peaks,which(Pm_onepop[[i]][,t]>thresh))
	peaks = length(where_peaks)
	while(peaks!=num_peaks[i]){
		where_peaks = local_max(Pm_onepop[[i]][,t])
		where_peaks = intersect(where_peaks,which(Pm_onepop[[i]][,t]>thresh))
		peaks = length(where_peaks)
		t = t+1
	}
	t_peaks_mat[i] = t-1
}
## ----

growth_rate = as.list(1:P)
dim(growth_rate) = dim(Pm_keep)

for( i in 1:P){
	growth_rate[[i]] = array(0,c(Nm,Tsteps-1))
	for(t in 1:Tsteps-1){
		growth_rate[[i]][,t] = Pm_keep[[i]][,t+1] / Pm_keep[[i]][,t]
		# growth_rate[[i]][which(Pm_keep[[i]][,t]<1e-10),t] = 0
	}
}

s = 2
f = 3

preferences = array(0,dim=c(Nm,Nf,5))
female_tots = matrix(0,Nf,5)

t=50
for(m in 1:5){
	Pm = Pm_keep[[s,f,m]]
for(j in 1:Nf){
	y = frange[j]
	weight = dnorm(mrange,mean=y,sd=sqrt(sigma2))
	z = int(weight*Pm[,t])
	preferences[,j,m] = weight #preference given by each female to each male
	female_tots[j,m] = z #total preferences given by each female
}
}

m = 1; plot(Pf_keep[[s,f,m]][,1]/female_tots[,m],type='l',ylim=c(0,15),col=mypal[m]);for(m in 2:5){lines(Pf_keep[[s,f,m]][,1]/female_tots[,m],col=mypal[m],type='l')};legend(1700,1,legend=mmix_sigma2_vals,lty=1,col=mypal);lines(Pf_keep[[s,f,m]][,1])

m = 1; plot(1/female_tots[,m],ylim=c(0,15),type='l',col=mypal[m]);for(m in 2:5){lines(1/female_tots[,m],col=mypal[m],type='l')};legend(1700,0.5,legend=mmix_sigma2_vals,lty=1,col=mypal);lines(Pf_keep[[s,f,m]][,1])

m = 1; plot(female_tots[,m],ylim=c(0,1),type='l',col=mypal[m]);for(m in 2:5){lines(female_tots[,m],col=mypal[m],type='l')};legend(1700,0.5,legend=mmix_sigma2_vals,lty=1,col=mypal);lines(Pf_keep[[s,f,m]][,1])

t = t+10;m = 1; plot(growth_rate[[s,f,m]][,t],type='o',ylim=c(0,10),col=mypal[m]);for(m in 2:5){lines(growth_rate[[s,f,m]][,t],col=mypal[m],type='o');};legend(1700,3,legend=mmix_sigma2_vals,lty=1,col=mypal);lines(Pf_keep[[s,f,m]][,1])

t = 1000;m = 1; plot(Pm_keep[[s,f,m]][,t],type='l',ylim=c(0,2),col=mypal[m]);for(m in 2:5){lines(Pm_keep[[s,f,m]][,t],col=mypal[m])};legend(1700,1,legend=mmix_sigma2_vals,lty=1,col=mypal);lines(Pf_keep[[s,f,m]][,1])



## ---- explain difference growth rates -----------------
s = 2
sigma2 = sigma2_vals[s]
f = 3
m = 2
Pm = Pm_keep[[s,f,m]]
Pf = Pf_keep[[s,f,m]]
t = 1
peaks = 2 
while(peaks!=num_peaks[s,f,m]){
	where_peaks = local_max(Pm_keep[[s,f,m]][,t])
	where_peaks = intersect(where_peaks,which(Pm_keep[[s,f,m]][,t]>1e-4))
	peaks = length(where_peaks)
	t = t+1
}
t_peak = t-1
i1 = peaks_now[1]
i2 = peaks_now[1]+5


# break down the mating and preference probabilities
c = 1e-10 #recognition cutoff
recog_d = sqrt(-2*sigma2*log(sqrt(2*pi*sigma2)*c)) 
#^distance / difference at which a female can recognize a male

preferences = array(0,dim=c(Nm,Nf))
female_tots = matrix(0,Nf)

for(j in 1:Nf){
	y = frange[j]
	# x1 = y-recog_d
	# x2 = y+recog_d
	# w1 = which(mrange>=x1)
	# w2 = which(mrange<=x2)
	# recognize[j] = int(Pm[intersect(w1,w2),t_peak])
	# weight = 1/sqrt(2*pi*sigma2)*exp(-(mrange-y)^2/(2*sigma2))
	weight = dnorm(mrange,mean=y,sd=sqrt(sigma2))
	# weight = matrix (0,Nf,1)
	# weight[c(f0,x1)] = 1
	# weight[j] = 1+alpha
	z = int(weight*Pm[,t_peak])
	# if(z!=0){
		# pxy[,j] = Pf[j,t_peak]*weight*Pm[,t_peak]/z
		# }
	preferences[,j] = weight #preference given by each female to each male
	female_tots[j] = z #total preferences given by each female
}

growth_rate = array(0,dim=c(Nm,1))
preferredby = array(0,dim=c(Nm,1))
recognizedby = array(0,dim=c(Nm,1))

for(i in 1:Nm){
	x = mrange[i]
	y1 = x - recog_d
	y2 = x + recog_d
	w1 = which(frange>=y1)
	w2 = which(frange<=y2)
	preferredby[i] = int(Pf[intersect(w1,w2),t_peak]*
		dnorm(x-frange[intersect(w1,w2)],mean=0,sd=sqrt(sigma2))) 
	recognizedby[i] = int(Pf[intersect(w1,w2),t_peak])
	#^how many females recognize each male, weighted by their preferences
	growth_rate[i] = int(preferences[i,]*Pf[,t_peak]/female_tots)
}
growth_rate[Pm[,t_peak]==0]=0

##---

m=2;radius=sqrt(step*rev(as.vector(four_freq[,,m]+.05)));symbols(rev(as.vector(smat[,,m])),rev(as.vector(var_mat[,,m])),circles=radius,bg=rev(as.vector(colmat[,,m])),inches=.5*max(radius),xlab='Selectivity',ylab='Entropy');legend(x=1.5,y=01,legend=(fmix_sigma2_vals),lty=0,pch=16,col=(mypal),bty='n')


m=2;radius=sqrt(rev(as.vector(1/wavelength_peaks[,,m]+.05)));symbols(rev(as.vector(smat[,,m])),rev(as.vector(ent_mat[,,m])),circles=radius,bg=rev(as.vector(colmat[,,m])),inches=.1*max(radius),xlab='Selectivity',ylab='Entropy');legend(x=1.5,y=01,legend=(fmix_sigma2_vals),lty=0,pch=16,col=(mypal),bty='n')





