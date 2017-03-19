trait_chunk_num = 301
sigmay2 = 3
sigmay2_init = sigmay2
sigmax2 = 0.8
sigmax2_init = sigmax2
pf = 1
pm = 1 
rho = 0.6
minweight = 10^(-320)
mut_prob = 0
steps = 200
store = 1
source('range_setup.R')


###
trait_chunk_num = 301
sigmay2 = 2
sigmax2 = 0.8
rho = 0.6
steps = 20
source('range_setup.R')

f_init = dnorm(frange,fmin,sqrt(sigmay2_init))
f_init[f_init==0] = 10^max(floor(log(min(f_init[which(f_init>0)]),base=10)),-320)
f_init = f_init/sum(f_init)
f_init = pf*f_init+(1-pf)*rev(f_init)
m_init = dnorm(mrange,mmin,sqrt(sigmax2_init))
m_init[m_init==0] = 10^max(floor(log(min(m_init[which(m_init>0)]),base=10)),-320)
m_init = m_init / sum(m_init)
m_init = pf*m_init+(1-pf)*rev(m_init)

continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
v = continuous_weight

pref_chunk_num = 3

k1 = 35
k2 = 51

k1 = 97
k2 = 51

k1 = 55
k2 = 51

k1 = 35
k2 = 35

chunk_vec = c(rep(1,(Nm-k1-2*k2)/2),rep(2,k2),rep(3,k1),rep(4,k2),rep(5,(Nm-k1-2*k2)/2))

n = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) length(which(chunk_vec==x)))
n[1:(chunk_vec[midpt]-1)] = 2*n[1:(chunk_vec[midpt]-1)]

s = apply(matrix(1:chunk_vec[midpt],nrow=1),2,function(x) sum((mrange[which(chunk_vec==x)]+1)^2))
s[1:(chunk_vec[midpt]-1)] = 2*s[1:(chunk_vec[midpt]-1)]

m = rbind(n,s,c(1,0,0))

# sigma2_vals = seq(1.5,sigmay2,length.out=8)

# equilibrium_step = as.list(1:(2*length(sigma2_vals)))
# dim(equilibrium_step) = c(length(sigma2_vals),2)

# for(i in 1:length(sigma2_vals)){
	sigma2 = 3
	v = c(1,sigma2,0)
	p = solve(m,v)
	p[1] = 0
	print(p)
	# if(length(which(p<0))!=0){
			# p = solve(m,c(1,sigma2,0.001,0))
	# }
	# if(length(which(p<0))!=0){
			# p = solve(m,c(1,sigma2,0,0))
	# }
	fixed_weight1 = apply(matrix(chunk_vec[1:midpt],nrow=1),2,function(x) p[x])
	fixed_weight1 = c(fixed_weight1,rev(fixed_weight1[1:(length(fixed_weight1)-1)]))
	
	continuous_weight = dnorm(mrange,mean=mrange[midpt],sd=sqrt(sigma2)) 
	fixed_weight2 = continuous_weight/sum(continuous_weight)
	
	plot(mrange,fixed_weight1);points(mrange,fixed_weight2,col='red')
	
	fixed_weight = fixed_weight1
	p1 = dynamics()
			
	fixed_weight = fixed_weight2
	p2 = dynamics()	
	
	plot(mrange,p1$Pm[,steps]);points(mrange,p2$Pm[,steps],col='red')
	
	# equilibrium_step[[i,1]] = p1$Pm[,steps]
	# equilibrium_step[[i,2]] = p2$Pm[,steps]
# }