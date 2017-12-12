source('dynamics_by_mode_new_numbers.R')
source('recursion_all.R')
source('init_conds.R')
library(mvtnorm)
library(PearsonDS)

trait_chunk_num = 121
minweight = 10^(-320)
mut_prob = 0.0
source('range_setup_short.R')

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
rho = 0.9
excess_kurt = -0.5

song = 'norm'
pref = 'pearson'
func = 'norm'

k1 = 21
k2 = 15

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
f_init = f_init / sum(f_init)

store = 1
steps = 1000

d2 = dynamics_mode2_pearson()

store = 29900
steps = 30000

d3 = dynamics_mode3()

store=1
steps = 40

d5 = dynamics_mode5_pearson()

save(d2,d3,d5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,excess_kurt,file='/homes/ebrush/priv/song_learning_evolution/peaks_pearson_pref_dist.Rdata')

song = 'pearson'
pref = 'norm'
func = 'norm'

k1 = 21
k2 = 15

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

m_init = dpearson(mrange,moments=c(mean=-1,variance=sigmax2,skewness=0,kurtosis=3+excess_kurt))
			m_init = m_init / sum(m_init)

f_init = dpearson(mrange,moments=c(mean=-1,variance=sigmay2,skewness=0,kurtosis=3+excess_kurt))
f_init = f_init / sum(f_init)

store = 1
steps = 1000

d2 = dynamics_mode2_pearson()

store = 29900
steps = 30000

d3 = dynamics_mode3()

store=1
steps = 40

d5 = dynamics_mode5_pearson()

save(d2,d3,d5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,excess_kurt,file='/homes/ebrush/priv/song_learning_evolution/peaks_pearson_song_dist.Rdata')