source('dynamics_by_mode_new_numbers.R')
source('recursion_all.R')
source('init_conds.R')
library(mvtnorm)

trait_chunk_num = 121
minweight = 10^(-320)
mut_prob = 0.0
source('range_setup_short.R')

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
rho = 0.9

song = 'norm'
pref = 'step'
func = 'norm'

k1 = 21
k2 = 15

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

store = 29900
steps = 30000

d3 = dynamics_mode3()

store=1
steps = 40

d5 = dynamics_mode5()

save(d3,d5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/peaks_step_pref_dist.Rdata')

song = 'step'
pref = 'norm'
func = 'norm'

k1 = 21
k2 = 15

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

store = 29900
steps = 30000

d3 = dynamics_mode3()

store=1
steps = 40

d5 = dynamics_mode5()

save(d3,d5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/peaks_step_song_dist.Rdata')