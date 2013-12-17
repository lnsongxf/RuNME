## External Libraries

# Numpy: General Numerical and Matrix Library
import numpy as np

# Scipy.optimize: Numerical Optimization
from scipy import optimize

# Scipy.interpolate: (Polynomial) Interpolation
from scipy import interpolate

# Scipy.integrate: Numerical Integration
from scipy import integrate

# numpy.polynomial.polynomial: Handles polynomials
from numpy.polynomial import polynomial

import sys

# From Algorithm 12.5 in Judd 1998, we have the following structure:
# Objective: Solve the Bellman equation
# Step 0. Initialization: Choose function form for Vhat and approximation grid X.
#						  Make initial guess for Vhat and choose stopping criterion
# Step 1. Maximization: Compute v_j using TVhat for each x_j in X
# Step 2. Fitting: Compute a^i+1 such that Vhat approximates the data (v_i, x_i)
# Step 3. Check if stopping criterion met, otherwise go to step 1.

class parametricDP:

	def __init__(self, maximize, fitting, calcError):
		self.tol = 0.01
		self.error = self.tol + 1
		self.maxIt = 300
		self.n = 10
		self.deg = 3

		self.alpha = 0.33
		self.beta = 0.95
		self.depr = 0.1
		self.s = 2
		self.min = 1.0
		self.max = 7.0

		self.a = np.zeros(self.deg + 1)
		self.v = np.zeros(self.n)

		self.grid = np.arange(self.min, self.max, (self.max - self.min) / self.n)

		self.maximize = maximize
		self.fitting = fitting
		self.calcError = calcError

		# compute initial guess
		# guess = no savings
		guessF = lambda k: (1 / (1 - self.s)) * (k**self.alpha - self.depr * k) ** (1 - self.s)

		for i in xrange(self.n):
			# no savings - kprime = k0
			self.v[i] = guessF(self.grid[i])

		self.a = polynomial.polyfit(self.grid, self.v, self.deg)

		print "INITIAL VALUE FUNCTION",self.v
		print "INITIAL APPROX", self.a

	def solve(self):
		# Solve the dynamic programming problem
		# described in parametricDP

		it = 0

		while self.error > self.tol and it < self.maxIt:
			# Save previous a
			a0 = np.copy(self.a)

			# Maximize
			res = self.maximize(self)
			self.v = res[0]
			self.policy = res[1]

			# Fit
			self.a = self.fitting(self)

			# Update error
			self.error = self.calcError(self, self.a, a0)

			print "Iteration", it

			#print "NEW VALUE FUNCTION", self.v
			#print "NEW APPROX", self.a
			
			print "Error", self.error
			print "\n"

			it = it + 1

		return (self.a,self.v)

# Problem Description
# consumption problem as in dypro.py

optParam = {}
optParam['tol'] = 0.01
optParam['maxIt'] = 300
optParam['n'] = 99.0
optParam['deg'] = 5

param = {}
param['alpha'] = 0.33
param['beta'] = 0.95
param['depr'] = 0.1
param['s'] = 2
param['min'] = 1.0
param['max'] = 7.0

initialA = np.zeros(optParam['deg'] + 1)
initialV = np.zeros(optParam['n'])

# Functions



# def maximize0(a, s, alpha, beta, depr, grid, n):
# 	# Find the value of v that maximizes 
# 	# given self.a

# 	# set up objective function

# 	def utility(consumption):
# 		return (1 / (1 - s)) * consumption**(1 - s)

# 	def production(k0):
# 		return k0 ** alpha + (1 - depr) * k0

# 	def obj(kprime):
# 		# objective function for optimization problem

# 		# Compute current period value
# 		consumption = production(k0) - kprime
# 		curr = utility(consumption) # consumption in current period

# 		# Compute continuation value
# 		# use x_i point		
# 		cont = 0
# 		val = polynomial.polyval(kprime, a)

# 		return -(curr + beta * cont)

# 	v = np.zeros(n)

# 	for i in xrange(n):

# 		k0 = grid[i]

# 		x0 = 0
# 		minbds = grid[0]
# 		maxbds = min(grid[-1], production(grid[i])**alpha + (1 - depr) * grid[i])
# 		bounds = ((minbds, maxbds),)

# 		res = optimize.minimize(obj, x0=x0, method='L-BFGS-B', bounds=bounds)
# 		ubest = res.x
# 		v[i] = -obj(ubest)

# 	return v

def maximize(self):
	# Find the maximum value of v
	# given self.a

	# set up objective function

	def utility(consumption):
		return consumption**(1 - self.s) / (1 - self.s)

	def production(k0):
		return k0 ** self.alpha + (1 - self.depr) * k0

	def obj(kprime, k0):
		# objective function for optimization problem

		# Compute current period value
		curr = utility(production(k0) - kprime) # consumption in current period

		# Compute continuation value
		# value of next period being kprime (since control = next pd state  in this model)
		# use x_i point		
		cont = 0
		val = polynomial.polyval(kprime, self.a)
		cont = val

		return -(curr + self.beta * cont)

	v = np.zeros(self.n)

	policy = np.zeros(self.n)

	# for each point x_j, compute new v_j
	for i in xrange(self.n):

		# Compute the value of the maximizatin problem
		# given a continuation function described by
		# polynomials w/ coefficients a, and 
		# capital is k0 = self.grid[i]

		k0 = self.grid[i]
		
		minbds = self.grid[0]
		maxbds = min(self.grid[-1], self.grid[i]**self.alpha + (1 - self.depr) * self.grid[i])
		x0 = (minbds + maxbds) / 2.0 # pick a better starting point
		bounds = ((minbds, maxbds),)

		# Perform optimization to find best control
		res = optimize.minimize(lambda kprime: obj(kprime, k0), x0=x0, method='L-BFGS-B', bounds=bounds)

		# Check optimization results
		if res.success != True:
			print "Optimization Failed"
			print res
			sys.exit()

		u = res.x

		policy[i] = u
		v[i] = -obj(u,k0)
		self.v[i] = -obj(u,k0)

	#print policy

	return (v, policy)

grid = np.arange(param['min'], param['max'], (param['max'] - param['min']) / optParam['n'])

# def fitting0(v, grid, deg):
# 	# Find the value of a that fits 
# 	# given self.v

# 	# Least Squares Approximation with polynomials
# 	coef = polynomial.polyfit(grid, v, deg)

# 	# Approximation with Chebyshev polynomials
# 	# TODO

# 	return coef

def fitting(self):
	coef = polynomial.polyfit(self.grid, self.v, self.deg)
	return coef

# def calcError0(a, aprime, n, grid):
# 	# Evaluates the polynomials described by coefficients
# 	# a and aprime at each point on the grid
# 	# and calculates the norm of the difference between
# 	# the two vectors

# 	diff = np.zeros(n)
# 	for i in xrange(n):
# 		# evaluate polynomials
# 		diff[i] = polynomial.polyval(grid[i], a) - polynomial.polyval(grid[i], aprime)
	
# 	return np.sqrt(np.dot(diff,diff))

def calcError(self, aprime, a):
	# Evaluates the polynomials described by coefficients
	# a and aprime at each point on the grid
	# and calculates the norm of the difference between
	# the two vectors

	diff = np.zeros(self.n)
	for i in xrange(self.n):
		# evaluate polynomials
		diff[i] = polynomial.polyval(self.grid[i], a) - polynomial.polyval(self.grid[i], aprime)
	
	return np.sqrt(np.dot(diff,diff))

dp = parametricDP(maximize, fitting, calcError)
res = dp.solve()
#print res[0]
#print res[1]