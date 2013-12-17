import numpy as np
import scipy.optimize as opt
import sys
import matplotlib.pyplot as ppt

"""

A package for solving numerical dynamic programming problems
Runnan Yang 2013

The Bellman equation for an infinite-horizon problem can be written as a system of nonlinear equations:
V_i = max_u [\pi(x_i, u) + \beta \sum_j=1^n q_ij(u) V_j] for i = 1, ..., n
(Judd)

***** FUNCTIONS *****

	Finite State Methods

	valueFunctionIteration() - finds optimal value function via value function iteration
	valueFunctionIterationStep() - steps through single value function iteration

***** EXAMPLES *****

	wealthAccumulation() - deterministic growth model
	wealthAccumulationExample() - prototype for solutions

"""

def valueFunctionIteration():
	# 0: Set up parameters

	# 1: For i = 1 ... n, compute 
	# V_i^{l+1} = max_u \pi(x_i, u) + \beta \sum_j=(1,n) q_ij(u) V_j^l

	# 2: If ||V^l+1 - V^l|| < \epsilon, go to step 3, else go to step 1

	# 3: Compute final solution
	# U* = U V^{l+1}
	# P*_j = \pi(x_i, U*_i)
	# V* = (I-\beta Q^U)^-1 P*
	return 0

def valueFunctionIterationStep(v0, ):

	# 1. Construct function to be maximized

	# 2. Construct maximization problem

	# 3. Solve maximization problem

	return 0

def wealthAccumulation(utilityFn, beta, outputFn):

	"""
	utilityFn(c_t) - compute utility given consumption
	beta - period by period discount factor
	outputFn(k_t) - output given capital stock k_t, net of depreciation
	              - k_{t+1} = k_t + f(k_t) - c_t
	"""

	problem = {}
	#valueFn = 
	return 0

def wealthProblem(valueFnIt = False, policyFnIt = False):
	# 0. Parameters

	# Model Parameters
	alpha = 0.33
	beta = 0.95
	depr = 0.1
	s = 2

	# Optimization Parameters
	tol = 0.01
	maxIts = 300
	dif = tol + 1000
	its = 0
	kgrid = 99

	# utility of consumption
	def utility(c):
		return (1 / (1 - s)) * (c**(1-s) - 1)

	# production and capital net depreciation
	def production(k0):
		return k0**alpha + (1 - depr) * k0

	def steadyState():
		# 1. solve for steady state
		kstar = (alpha/(1/beta - (1 - depr)))**(1/(1-alpha)) # capital at ss
		cstar = kstar**alpha - depr * kstar
		istar = depr * kstar
		ystar = kstar**alpha

		return (kstar, cstar, istar, ystar)

	# We use the steady state value as the guess for our optimization
	# as well as to inform our kmin and max values
	kstar = steadyState()[0]

	# state bounds (within 75% of the steady state)
	kmin = 1.00 #0.25 * kstar
	kmax = 7.00 #1.75 * kstar
	bnds = ((kmin, kmax),)

	# Grid Parameters
	kgrid = 99
	grid = (kmax - kmin) / kgrid

	# Grid of points from kmin to kmax, at intervals of length grid
	kmat = np.transpose(np.arange(kmin, kmax, grid))

	# Number of periods
	n = len(kmat)

	# Temporary Value Function Values
	v0 = np.zeros(n)
	v1 = np.zeros(n)
	kl1 = np.zeros(n)
	p1 = np.zeros(n)

	# 2. define value function
	def valFn(k):

		# Note: uses the following external variables
		# kmat - grid of control points
		# v0 - initial value function estimate - **this changes every round**
		# k0 - current period capital - **this changes every round**
		# utility - utility function
		# production - production function

		# 2.1 DISCRETIZATION

		# Get index of discrete level of 
		# capital below (lo) and above (hi)
		khi = np.argmax(kmat > k)
		klo = khi - 1
		
		# **value function evaluated at capital k**
		# assuming linear form between klo and khi
		gg = v0[klo] + (k - kmat[klo]) * (v0[khi] - v0[klo]) / (kmat[khi] - kmat[klo])

		# 2.2 Compute Value Function
		c = production(k0) - k
		if c <= 0:
			val = -9999999 - 999 * abs(c)
		else:
			val = utility(c) + beta * gg

		return -val

	# 3. Value Function Iteration
	#
	if valueFnIt == True:
		while dif > tol and its < maxIts:

			# For each possible control
			# k_{t+1}
			for i in xrange(n):
				# current state
				k0 = kmat[i]

				# find the control which maximizes the value function
				res = opt.minimize(valFn, x0=(kstar), method='L-BFGS-B', bounds=bnds)
				k1 = res.x

				# set value function
				v1[i] = -valFn(k1)
				kl1[i] = k1

			# calculate norm of new and old value function
			# and update values (careful w/ copying values!)
			vdiff = v1 - v0
			dif = np.sqrt(np.dot(vdiff, vdiff))
			v0 = v0 + vdiff
			its = its + 1

			print "Norm",dif
			print "Finished Iteration",its

	# 4. Policy Function Iteration
	#
	if policyFnIt == True:

		# Run one iteration of value function to get initial value function
		# and policy function
		for i in xrange(n):
			# Current capital k0
			k0 = kmat[i]

			# Find optimal level of k_{t+1}
			# given current iteration of value function
			res = opt.minimize(valFn, x0=(kstar), method='L-BFGS-B', bounds=bnds)
			k1 = res.x

			v1[i] = -valFn(k1)
			kl1[i] = k1

		v0 = np.copy(v1)

		while dif > tol and its < maxIts:

			# 1. Compute next policy function iteration
			for i in xrange(n):
				# Current capital k0
				k0 = kmat[i]

				# Find optimal level of k_{t+1}
				# given current iteration of value function
				res = opt.minimize(valFn, x0=kstar, method='L-BFGS-B', bounds=bnds)
				kl1[i] = res.x
				
			print "kl1",kl1 # current policy function iteration

			# 2. Compute returns for each state
			for i in xrange(n):
				consumption = production(kmat[i]) - kl1[i]
				p1[i] = utility(consumption)
			print "p1",p1

			# 3. Compute next value function iteration

			# Create Q transition matrix given policy u
			# in this case, policy u goes to state u (control is k_{t+1})

			Q = np.zeros((n,n))
			for i in xrange(n):
				# find the entry just below
				klo = np.argmax(kmat > kl1[i])
				khi = klo + 1
				if kl1[i] % grid > grid/2:
					Q[i,khi] = 1
				else:
					Q[i,klo] = 1

			#print "to invert",np.eye(n) - beta * Q
			#print "det",np.linalg.det( np.eye(n) - beta * Q)

			v1 = np.asarray(np.linalg.inv(np.eye(n) - beta * Q) * np.matrix(p1).T)
			v1 = v1.reshape(len(v1))
			#print type(v1), type(v0)
			#print "v1",v1
			print "v0",v0

			# calculate norm of new and old value function
			# and update values (careful w/ copying values!)
			vdiff = v1 - v0
			#print "vdiff", vdiff
			#print np.multiply(vdiff, vdiff)
			dif = np.sqrt(np.dot(vdiff, vdiff))
			v0 = np.copy(v1)
			its = its + 1

			print "Norm",dif
			print "Finished Iteration",its



	return (kmat,v0,kl1)

if __name__ == "__main__":
    print "LOADING DYNAMIC PROGRAMMING PACKAGE"
    x = wealthProblem(valueFnIt = True)
    kmat = x[0]
    v0 = x[1]
    kl1 = x[2]
    ppt.clf()
    ppt.plot(x[0],x[1])
    ppt.savefig('myfig')