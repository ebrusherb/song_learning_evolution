source('dynamics_by_mode_new_numbers.R')
source('recursion_all.R')
source('init_conds.R')
library(mvtnorm)

trait_chunk_num = 205
source('range_setup_long.R')
minweight = 10^(-320)
mut_prob = 0.0

store = 1

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
rho = 0.9

# # song = 'norm'
# pref = 'step'
# func = 'norm'

# k1 = 7
# k2 = 7

# ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
# m_init = ic$m_init
# f_init = ic$f_init
# fixed_weight = ic$fixed_weight

# steps = 200

# var_mat2  = array(NA,c(2,steps+1))

# d2 = dynamics_mode2()
# var_mat2[1,] = apply(d2$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat2[2,] = apply(apply(d2$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# steps = 10

# var_mat5  = array(NA,c(2,steps+1))

# d5 = dynamics_mode5()
# var_mat5[1,] = apply(apply(d5$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat5[2,] = apply(apply(d5$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_pref_dist2.Rdata')

# # song = 'step'
# pref = 'norm'
# func = 'norm'

# k1 = 7
# k2 = 7

# ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
# m_init = ic$m_init
# f_init = ic$f_init
# fixed_weight = ic$fixed_weight

# steps = 500

# var_mat2  = array(NA,c(2,steps+1))

# d2 = dynamics_mode2()
# var_mat2[1,] = apply(d2$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat2[2,] = apply(apply(d2$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# steps = 10

# var_mat5  = array(NA,c(2,steps+1))

# d5 = dynamics_mode5()
# var_mat5[1,] = apply(apply(d5$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat5[2,] = apply(apply(d5$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_song_dist2.Rdata')

# song = 'norm'
# pref = 'norm'
# func = 'norm'

# k1 = NA
# k2 = NA

# ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
# m_init = ic$m_init
# f_init = ic$f_init
# fixed_weight = ic$fixed_weight

# steps = 50

# var_mat2  = array(NA,c(2,steps+1))

# d2 = dynamics_mode2()
# var_mat2[1,] = apply(d2$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat2[2,] = apply(apply(d2$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# steps = 10

# var_mat5  = array(NA,c(2,steps+1))

# d5 = NA

# save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_normal2.Rdata')

# d5 = dynamics_mode5()
# var_mat5[1,] = apply(apply(d5$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat5[2,] = apply(apply(d5$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_normal2.Rdata')
# save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/Users/eleanorbrush/Desktop/all_modes_numerically_all_ normal2.Rdata')

song = 'norm'
pref = 'norm'
func = 'step'

k1 = 3
k2 = 5

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

steps = 50

var_mat2  = array(NA,c(2,steps+1))

d2 = dynamics_mode2()
var_mat2[1,] = apply(d2$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat2[2,] = apply(apply(d2$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

steps = 10

var_mat5  = array(NA,c(2,steps+1))

d5 = dynamics_mode5()
var_mat5[1,] = apply(apply(d5$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
var_mat5[2,] = apply(apply(d5$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

save(d2,d5,var_mat2,var_mat5,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_pref_func2.Rdata')
