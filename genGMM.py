import numpy as np
import matplotlib.pyplot as plt
from dpm import * 

class Hypers:
	def __init__(self,dim,mu_0,Lambda_0,kappa_0,nu_0):
		self.dim = dim
		self.mu_0 = mu_0
		self.Lambda_0 = Lambda_0
		self.kappa_0 = kappa_0
		self.nu_0 = nu_0

def genG(theta,sigma,num):
	r = np.random.multivariate_normal(theta,sigma,num)
	return r

def genGMM(pi,theta,sigma,num):
	n = len(pi)
	r= np.random.choice(np.arange(n),num, p = pi)
	print r
	num = len(r[np.where(r == 0)])
	data =genG(theta[0],sigma[0],num)
	print data
	t = 0
	for c in range(1,n):
		num = len(r[np.where(r == c)])
		data =np.concatenate((data,genG(theta[c],sigma[c],num)))
		print data
	return data

pi = [0.3,0.3,0.4]
theta = np.array([[0.0,10.0],[10.0,0.0],[-10.0,-10.0]])
print theta
sigma = np.array([[[0.2,0.01],[0.01,0.1]],[[0.2,0.05],[0,0.1]],[[0.1,0.05],[0,0.2]]])

data = genGMM(pi,theta,sigma,50)
print data
#x = np.transpose(data)[0]
#y = np.transpose(data)[1]

#plt.scatter(x, y)
#plt.show()

print "111"
mu_0 = np.array([0.0,0.0])
kappa_0 = 1;
dim =2
Lambda_0 = np.array([[1.0,0.0],[0.0,1.0]])
nu_0 = dim;
hypers = Hypers(dim,mu_0,Lambda_0,kappa_0,nu_0)
print mu_0
print hypers.mu_0
z = dpm_gibbs_sampler(data,10,hypers,init_K= 20)
print z
print data
