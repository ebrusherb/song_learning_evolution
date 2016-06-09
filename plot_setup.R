library(ggplot2)
library(RColorBrewer)
library(reshape)
library(gridExtra)
source('local_max.R')
source('entropy.R')
source('ind2sub.R')
source('int.R')
source('get_legend.R')
source('range_setup.R')


d = dim(Pm_onepop)
P = prod(d)
Tend=d[2]


Ns = length(sigma_vals)
Nfs = length(f_sigma_vals)
Nms = length(m_sigma_vals)
Nmp = length(mut_prob_vals)

subset = 601:1201;
subset = 1:Nm;
Tend = dim(Pm_onepop[[1]])[2]

thresh = 5e-4
half = floor(Nm/2)

four_freq = array(NA,dim(Pm_onepop))
four_mat = list()
num_peaks = array(NA,dim(Pm_onepop))
var_mat_m = array(NA,c(dim(Pm_onepop)))
ex_mat_m = array(NA,dim(Pm_onepop))
var_mat_f = array(NA,c(dim(Pf_onepop)))
ent_mat = array(Inf,dim=dim(Pm_onepop))



for(i in 1:P){
	eq_pop = Pm_onepop[[i]][,Tend]
	eq_pop = eq_pop/sum(eq_pop)
	eq_pop_f = Pf_onepop[[i]][,Tend]
	eq_pop_f = eq_pop_f/sum(eq_pop_f)
	fourComps = fft(eq_pop[subset]);
	fourCoeffs = abs(fourComps);
	coeffs = fourCoeffs / (length(subset)/2);
	coeffs[1] = fourCoeffs[1] / length(subset)	
	mf = mean(coeffs[1:200])
	ml = intersect(local_max(round(coeffs[1:half],5)),which(coeffs[1:half]>mf/2))
	if(length(ml)==1 || length(ml)==0){four_freq[i] = 1} else{ four_freq[i] = median(diff(ml))}
	four_mat[[i]] = coeffs[1:half]
	m = local_max(eq_pop)
	m = intersect(m, which(eq_pop>thresh))
	l = length(m)
	num_peaks[i] = l
	sub = ind2sub(dim(Pm_onepop),i)
	ex = sum(mrange*eq_pop)
	vx = sum((mrange-(ex))^2*eq_pop)
	ex_mat_m[i] = ex
	var_mat_m[i] = vx
	ex = sum(mrange*eq_pop_f)
	vx = sum((mrange-(ex))^2*eq_pop_f)
	var_mat_f[i] = vx
	ent_mat[i] = ent(eq_pop)
}
dim(four_mat) = dim(Pm_onepop)


smat= array(0,dim=dim(Pm_onepop))
fmat= array(0,dim=dim(Pm_onepop))
mmat= array(0,dim=dim(Pm_onepop))
pmat = array(0,dim=dim(Pm_onepop))

for(i in 1:P){
	sub = ind2sub(dim(Pm_onepop),i)
	smat[i] = sigma_vals[sub[1]]
	fmat[i] = f_sigma_vals[sub[2]]
	mmat[i] = m_sigma_vals[sub[3]]
	pmat[i] = mut_prob_vals[sub[4]]
}