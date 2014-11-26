import numpy as np
from scipy import linalg
from scipy.special import gamma
from scipy.special import gammaln


def mvtpdf(x,mu,sigma,nu):
	
	d = len(x)
	r = 0
	r = r -0.5 * np.log(linalg.det(sigma))
	print r 
	r = r -(d/2)*np.log(nu*3.14)
	print r 
	r = r + gammaln((nu+d)/2) - gammaln(nu/2)
	print r 
	r = r -((nu+d)/2) * np.log(1 +np.dot(np.dot((x-mu) ,linalg.inv(sigma)),(x-mu) ) /nu)
	print r 
	if r < -20:
		return 0.0
	else:
		return np.exp(r)

def dpm_gibbs_sampler(x,alpha,hypers,init_K=3,max_K=200,iter_num=20):
	
	N = len(x)
	dim = len(x[0])
	print dim
	K = init_K

	#init cluster by kmeans
	#wait for imp

	#book keeper data structure
	z =  np.random.choice(np.arange(K),N,p=np.ones(K)/K)

	mu = np.zeros(max_K*dim).reshape((max_K,dim))
	Sigma =  np.zeros(max_K*dim*dim).reshape((max_K,dim,dim))
	#sufficient statistics
	E_x = np.zeros(max_K*dim).reshape((max_K,dim))
	E_xx = np.zeros(max_K*dim*dim).reshape((max_K,dim,dim))
	E_n = np.zeros(max_K)
	

	#sufficient statistics
	for i in range(0,N):
		E_x[z[i]] = E_x[z[i]] + x[i]
		E_xx[z[i]] = E_xx[z[i]] + np.outer(x[i],x[i])
		E_n[z[i]] = E_n[z[i]] + 1

	for j in range(0,K):
		mu[j] = E_x[j] / E_n[j]
		Sigma[j] = E_xx[z[i]] / E_n[j] - np.outer(mu[j],mu[j])
	#main infer process
	for c in range(0,iter_num):
		print "iter" + str(c)
		#gen a permutation R
		R = np.random.permutation(N)
		for r in R:
			p = np.zeros(K+1)
			need_clean = False
			#sample x_r
			k = z[r]
			print "current k = " + str(k)
			E_n[k] = E_n[k] - 1.0
			E_x[k] = E_x[k] -  x[r]
			E_xx[k] = E_xx[k] - np.outer(x[r],x[r])
			
			for j in range(0,K):
				#compute the piror
				p_z = E_n[j]
				# f[j] = Guassain(x[r],mu[k],Sigma[k])
				if j != k:
					f = mvtpdf(x[r],mu[k],Sigma[k],hypers.nu_0 + E_n[k])
				else:	
					if E_n[k]!= 0:
						#compute the likelihood
						kappa_n = hypers.kappa_0 + E_n[k]


						mu_n = (hypers.kappa_0 * hypers.mu_0 +E_x[k])/kappa_n

						nu_n = hypers.nu_0 + E_n[k]

						SS = E_xx[k] - np.outer(E_x[k],E_x[k]) / E_n[k]
						
						print "SS:"
						print SS

						bias =  E_x[k] / E_n[k] - hypers.mu_0
						print "bias:"
						print bias

						Lambda_n =hypers.Lambda_0 + SS + np.outer(bias,bias) * hypers.kappa_0 * E_n[k]/kappa_n
						print "Lambda:"
						print Lambda_n
						mu[k] = mu_n
						Sigma[k] =  Lambda_n*(kappa_n+1)/(kappa_n* (nu_n - dim +1))	

						f = mvtpdf(x[r],mu[k],Sigma[k],nu_n)
					else:
						f = 0.0
						#could be optimized here
						need_clean = True
				p[j] = p_z * f

			p[K] = alpha * mvtpdf(x[r],hypers.mu_0,hypers.Lambda_0,hypers.nu_0)
			#new cluster likelihood 
			p[0:K+1] = p[0:K+1]/np.sum(p[0:K+1]) 
			print p
			z[r] = np.random.choice(np.arange(K+1),1,p=p[0:K+1])
			
			if z[r] == K:
				mu[K] = x[r]
				Sigma[K] = hypers.Lambda_0 
				E_n[K] =0
				E_x[K] = np.zeros(dim)
				E_xx[K] = np.zeros(dim*dim).reshape((dim,dim))

				K = K + 1

			if need_clean:
				print "need clean!!!   r:" +str(r)+"   z[r]:"+ str(z[r])
				print "current clean!!!" + str(K)
				for i in range(0,N):
					if(z[i] == k):
						print "z[i]==k!!!!"

						print "i:" + str(i) + ",z[i]:"+z[i]
					assert(z[i] != k)
					if z[i] > k:
						z[i] = z[i]-1
				for t in range(k,K-1):
					E_n[k] = E_n[k+1]
					E_x[k] = E_x[k+1]
					E_xx[k] = E_xx[k+1]
					mu[k] = mu[k+1]
					Sigma[k] =  Sigma[k+1]

				K = K-1
				E_n[K] =0
				E_x[K] = np.zeros(dim)
				E_xx[K] = np.zeros(dim*dim).reshape((dim,dim))

			E_n[z[r]] = E_n[z[r]] + 1.0
			E_x[z[r]] = E_x[z[r]] +  x[r]
			E_xx[z[r]] = E_xx[z[r]] + np.outer(x[r],x[r])
			print E_n[:K]
	return z;

