source('dynamics_by_mode_new_numbers.R')
source('recursion_all.R')
source('init_conds.R')
library(mvtnorm)

trait_chunk_num = 205
source('range_setup_long.R')
minweight = 10^(-320)
mut_prob = 0.0

steps = 1000
store = 1

sigmay2 = 2
sigmax2 = 1
sigma2 = 0.5
rho = 0.9

song = 'norm'
pref = 'step'
func = 'norm'

k1 = 7
k2 = 7

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

# var_mat  = array(NA,c(12,2,steps+1))

# d3 = dynamics_mode3()
# var_mat[3,1,] = apply(d3$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[3,2,] = apply(d3$Pf,2,function(v) sum((mrange+1)^2*v))

# d4 = dynamics_mode4()
# var_mat[4,1,] = apply(apply(d4$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat[4,2,] = apply(d4$Pf,2,function(v) sum((mrange+1)^2*v))

# d6 = dynamics_mode6()
# var_mat[6,1,] = apply(apply(d6$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat[6,2,] = apply(d6$Pf,2,function(v) sum((mrange+1)^2*v))

# d7 = dynamics_mode7()
# var_mat[7,1,] = apply(d7$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[7,2,] = apply(d7$Pf,2,function(v) sum((mrange+1)^2*v))

# d8=dynamics_mode8()
# var_mat[8,1,] = apply(d8$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[8,2,] = apply(apply(d8$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# d9=dynamics_mode9()
# var_mat[9,1,] = apply(d9$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[9,2,] = apply(d9$Pf,2,function(v) sum((mrange+1)^2*v))

# d11=dynamics_mode11()
# var_mat[11,1,] = apply(d11$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[11,2,] = apply(apply(d11$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# d12=dynamics_mode12()
# var_mat[12,1,] = apply(d12$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[12,2,] = apply(d12$Pf,2,function(v) sum((mrange+1)^2*v))

# save(d3,d4,d6,d7,d8,d9,d11,d12,var_mat,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_pref_dist.Rdata')

##
song = 'step'
pref = 'norm'
func = 'norm'

k1 = 7
k2 = 7

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

var_mat  = array(NA,c(12,2,steps+1))

d3 = dynamics_mode3()
# var_mat[3,1,] = apply(d3$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[3,2,] = apply(d3$Pf,2,function(v) sum((mrange+1)^2*v))

# d4 = dynamics_mode4()
# var_mat[4,1,] = apply(apply(d4$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat[4,2,] = apply(d4$Pf,2,function(v) sum((mrange+1)^2*v))

# d6 = dynamics_mode6()
# var_mat[6,1,] = apply(apply(d6$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
# var_mat[6,2,] = apply(d6$Pf,2,function(v) sum((mrange+1)^2*v))

# d7 = dynamics_mode7()
# var_mat[7,1,] = apply(d7$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[7,2,] = apply(d7$Pf,2,function(v) sum((mrange+1)^2*v))

# d8=dynamics_mode8()
# var_mat[8,1,] = apply(d8$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[8,2,] = apply(apply(d8$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# d9=dynamics_mode9()
# var_mat[9,1,] = apply(d9$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[9,2,] = apply(d9$Pf,2,function(v) sum((mrange+1)^2*v))

# d11=dynamics_mode11()
# var_mat[11,1,] = apply(d11$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[11,2,] = apply(apply(d11$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

# d12=dynamics_mode12()
# var_mat[12,1,] = apply(d12$Pm,2,function(v) sum((mrange+1)^2*v))
# var_mat[12,2,] = apply(d12$Pf,2,function(v) sum((mrange+1)^2*v))

# save(d3,d4,d6,d7,d8,d9,d11,d12,var_mat,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_song_dist.Rdata')

song = 'norm'
pref = 'norm'
func = 'norm'

k1 = NA
k2 = NA

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

var_mat  = array(NA,c(12,2,steps+1))

d3 = dynamics_mode3()
var_mat[3,1,] = apply(d3$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[3,2,] = apply(d3$Pf,2,function(v) sum((mrange+1)^2*v))

d4 = dynamics_mode4()
var_mat[4,1,] = apply(apply(d4$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
var_mat[4,2,] = apply(d4$Pf,2,function(v) sum((mrange+1)^2*v))

d6 = dynamics_mode6()
var_mat[6,1,] = apply(apply(d6$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
var_mat[6,2,] = apply(d6$Pf,2,function(v) sum((mrange+1)^2*v))

d7 = dynamics_mode7()
var_mat[7,1,] = apply(d7$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[7,2,] = apply(d7$Pf,2,function(v) sum((mrange+1)^2*v))

d8=dynamics_mode8()
var_mat[8,1,] = apply(d8$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[8,2,] = apply(apply(d8$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

d9=dynamics_mode9()
var_mat[9,1,] = apply(d9$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[9,2,] = apply(d9$Pf,2,function(v) sum((mrange+1)^2*v))

d11=dynamics_mode11()
var_mat[11,1,] = apply(d11$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[11,2,] = apply(apply(d11$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

d12=dynamics_mode12()
var_mat[12,1,] = apply(d12$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[12,2,] = apply(d12$Pf,2,function(v) sum((mrange+1)^2*v))

save(d3,d4,d6,d7,d8,d9,d11,d12,var_mat,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_normal.Rdata')

song = 'norm'
pref = 'norm'
func = 'step'

k1 = 3
k2 = 5

ic = init_conds(song,pref,func,sigmax2,sigmay2,sigma2,k1,k2)
m_init = ic$m_init
f_init = ic$f_init
fixed_weight = ic$fixed_weight

var_mat  = array(NA,c(12,2,steps+1))

d3 = dynamics_mode3()
var_mat[3,1,] = apply(d3$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[3,2,] = apply(d3$Pf,2,function(v) sum((mrange+1)^2*v))

d4 = dynamics_mode4()
var_mat[4,1,] = apply(apply(d4$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
var_mat[4,2,] = apply(d4$Pf,2,function(v) sum((mrange+1)^2*v))

d6 = dynamics_mode6()
var_mat[6,1,] = apply(apply(d6$Pm,c(1,3),sum),2,function(v) sum((mrange+1)^2*v))
var_mat[6,2,] = apply(d6$Pf,2,function(v) sum((mrange+1)^2*v))

d7 = dynamics_mode7()
var_mat[7,1,] = apply(d7$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[7,2,] = apply(d7$Pf,2,function(v) sum((mrange+1)^2*v))

d8=dynamics_mode8()
var_mat[8,1,] = apply(d8$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[8,2,] = apply(apply(d8$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

d9=dynamics_mode9()
var_mat[9,1,] = apply(d9$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[9,2,] = apply(d9$Pf,2,function(v) sum((mrange+1)^2*v))

d11=dynamics_mode11()
var_mat[11,1,] = apply(d11$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[11,2,] = apply(apply(d11$Pf,c(2,3),sum),2,function(v) sum((mrange+1)^2*v))

d12=dynamics_mode12()
var_mat[12,1,] = apply(d12$Pm,2,function(v) sum((mrange+1)^2*v))
var_mat[12,2,] = apply(d12$Pf,2,function(v) sum((mrange+1)^2*v))

save(d3,d4,d6,d7,d8,d9,d11,d12,var_mat,sigmax2,sigmay2,sigma2,k1,k2,trait_chunk_num,steps,file='/homes/ebrush/priv/song_learning_evolution/all_modes_numerically_step_pref_func.Rdata')
