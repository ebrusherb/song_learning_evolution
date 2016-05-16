library(RColorBrewer)
source('local_max.R')
source('entropy.R')
source('ind2sub.R')
source('int.R')

mypal=brewer.pal(5,'Set1')

d = dim(Pm_onepop)
P = prod(d)
Tend=d[2]


Ns = length(sigma2_vals)
Nfs = length(fmix_sigma2_vals)
Nms = length(mmix_sigma2_vals)
Nmp = length(mut_prob_vals)

subset = 601:1201;
subset = 1:Nm;
Tend = dim(Pm_onepop[[i]])[2]
half = floor(length(subset)/2)
four_freq = array(NA,dim(Pm_onepop))
# four_coef = array(NA,dim(Pm_onepop))
four_mat = array(NA,c(dim(Pm_onepop),half))
var_mat_m = array(NA,c(dim(Pm_onepop)))
var_mat_f = array(NA,c(dim(Pf_onepop)))
four_var = array(NA,c(dim(Pm_onepop)))

for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tend]
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourCoeffs[1] / length(subset)
	# v=order(coeffs[1:floor(length(subset)/2)],decreasing=TRUE)
	# # if(v[1]==2){four_freq[i]=v[2]} else{ four_freq[i]=v[1]}
	# four_coef[i] = coeffs[2]
	sub = ind2sub(dim(Pm_onepop),i)
	four_mat[cbind(matrix(rep(sub,half),nrow=half,byrow=TRUE),1:half)] = coeffs[1:half]
	ex = int(mrange*Pm_onepop[[i]][,Tend])
	vx = int((mrange-(ex))^2*Pm_onepop[[i]][,Tend])
	var_mat_m[i] = vx
	ex = int(mrange*Pf_onepop[[i]][,Tend])
	vx = int((mrange-(ex))^2*Pf_onepop[[i]][,Tend])
	var_mat_f[i] = vx
	mf = mean(coeffs[1:200])
	ml = intersect(local_max(round(coeffs[1:half],5)),which(coeffs[1:half]>mf/2))
	if(length(ml)==1 || length(ml)==0){four_freq[i] = 1} else{ four_freq[i] = median(diff(ml))}
	l = lm(log(coeffs[2:7]) ~ poly(1:6,2,raw=TRUE));
	four_var[i] = -summary(l)$coefficients[3,1]
}

num_peaks = array(NA,dim(Pm_onepop))
wavelength_peaks = array(Inf,dim(Pm_onepop))
ent_mat = array(Inf,dim=dim(Pm_onepop))

thresh = 1e-4

for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tend]
	m = local_max(eq_pop)
	m = intersect(m, which(eq_pop>thresh))
	l = length(m)
	num_peaks[i] = l
	if(l > 1){
		wavelength_peaks[i] = step*median(diff(m))
	}
	ent_mat[i] = ent(eq_pop)
}


smat= array(0,dim=dim(Pm_onepop))
fmat= array(0,dim=dim(Pm_onepop))
mmat= array(0,dim=dim(Pm_onepop))
pmat = array(0,dim=dim(Pm_onepop))
colmat = array(0,dim=dim(Pm_onepop))

for(i in 1:P){
	sub = ind2sub(dim(Pm_onepop),i)
	smat[i] = sigma2_vals[sub[1]]
	fmat[i] = fmix_sigma2_vals[sub[2]]
	mmat[i] = mmix_sigma2_vals[sub[3]]
	pmat[i] = mut_prob_vals[sub[4]]
	colmat[i] = mypal[sub[2]]
}

# p=1
# m=1;radius=sqrt(step*rev(as.vector(four_freq[,2:Nfs,m,p]+.05)));symbols(rev(as.vector(smat[,2:Nfs,m,p])),rev(as.vector(ent_mat[,2:Nfs,m,p])),circles=radius,bg=rev(as.vector(colmat[,1:(Nfs-1),m,p])),inches=.5*max(radius),xlab='Selectivity',ylab='Entropy');legend(x=1.5,y=01,legend=(fmix_sigma2_vals[2:Nfs]),lty=0,pch=16,col=(mypal[1:(Nfs-1)]),bty='n')

p=2;
m=2;radius=sqrt(step*(as.vector(four_freq[,2:Nfs,m,p]+.05)));symbols(log(as.vector(smat[,2:Nfs,m,p])),log(sqrt(as.vector(var_mat_m[,2:Nfs,m,p]))),circles=radius,bg=(as.vector(colmat[,1:(Nfs-1),m,p])),inches=.5*max(radius),xlab='Selectivity',ylab='Standard deviation');ylim=par('yaxp')[1:2];legend(x=1.5,y=mean(ylim),legend=(fmix_sigma2_vals[2:Nfs]),lty=0,pch=16,col=(mypal[1:(Nfs-1)]),bty='n');for(f in 2:Nfs){abline(h=log(fmix_sigma2_vals[f]),col=mypal[f-1])}
#,ylim=range(sqrt(fmix_sigma2_vals))+c(-.1,.1)

m=2;radius=sqrt(step*(as.vector(four_freq[,1:Nfs,m,p]+.05)));symbols((as.vector(smat[,1:Nfs,m,p])),(sqrt(as.vector(var_mat_m[,1:Nfs,m,p]))),circles=radius,bg=(as.vector(colmat[,1:(Nfs),m,p])),inches=.5*max(radius),xlab='Selectivity',ylab='Standard deviation');legend(x=1.5,y=0.1,legend=(fmix_sigma2_vals[1:Nfs]),lty=0,pch=16,col=(mypal[1:(Nfs)]),bty='n');for(f in 1:Nfs){abline(h=fmix_sigma2_vals[f],col=mypal[f])}

#####SOMETHING IS FUNNY WITH HOW I CALCULATED VARIANCE FOR s=1,m=1,f=1,p=2 and perhaps others!!! the expectation isn't -1 as it should be. screwing up bubble plots.

s = 3
m = 2
p = 2
subset = 1:Nm
plot(1:length(subset),1:length(subset),ylim=c(0,1),t='n')
for(f in 2:Nfs){
	lines(Pm_onepop[[s,f,m,p]][subset,Tend],col=mypal[f-1])
}
legend(200,0.8,fmix_sigma2_vals,col=mypal[1:(Nfs-1)],lty=1)