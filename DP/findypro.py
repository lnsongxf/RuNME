# findypro

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

class finiteDP:

	# Class for solving finite horizon dynamic programming problems

	def __init__(self):

		# Part One: Discretize and 
		# In progress

		# Part Two: Parameterization
		# To do

		## Optimization Parameters
		self.tol = 0.01
		self.error = self.tol + 1
		self.maxIt = 300
		self.n = 99

		## Model Parameters
		self.alpha = 0.33
		self.beta = 0.95
		self.delta = 0.1
		self.s = 2
		self.min = 1.0
		self.max = 7.0

		# time periods
		self.T = 60

		## Model functions
		self.utility = lambda c: c**(1 - self.s) / (1-self.s)
		self.terminalValue = lambda c: c
		self.production = lambda k: k**(self.alpha)

		# Grid of capital values
		self.grid = np.arange(self.min, self.max, (self.max - self.min) / self.n)

		self.memo = np.zeros((n+1, T))

	def solveDiscrete(self, k, t):

		# value to be maximized
		def obj(kprime, t):

			# lookup value if possible
			pos = np.argx[self.grid, kprime] - 1
			if self.memo[pos, t] != 0:
				return self.memo[pos, t]

			# compute current value
			self.utility = lambda c: c**(1 - self.s) / (1-self.s)
			self.production = lambda k: k**(self.alpha)

			consumption = self.production(k) + (1 - self.delta) * k - kprime
			current = self.utility(consumption)

			# Handle terminal state
			if t == self.T - 1:
				return current + self.beta * self.terminalValue(kprime)

			# compute value from being in state kprime at time t+1
			cont = solveDiscrete(kprime, t+1)

			return current + self.beta * cont
		
		maxSoFar = 1e-6
		kbest = 0
		for kprime in self.grid:
			tmp = self.obj(knext, t+1)

			if tmp > maxSoFar:
				maxSoFar = tmp
				kbest = knext

		return maxSoFar


	def solve(self, k, t):
		if t % 5 == 0:
			print k, ",", t
		def obj(kprime, t):

			# find value just below on grid
			khi = np.argmax(self.grid > kprime)
			klo = khi - 1

			self.utility = lambda c: c**(1 - self.s) / (1-self.s)
			self.production = lambda k: k**(self.alpha)

			# current value
			consumption = self.production(k) + (1 - self.delta) * k - kprime
			current = self.utility(consumption)

			# continuation value
			if t < self.T-1: # if not last period, call next period's optimization
				continuation = self.solve(kprime, t+1)
			else: # if last period, call terminal condition
				continuation = self.terminalValue(kprime)

			return -(current + self.beta * continuation)

		minbnd = self.min
		maxbnd = min(self.max, self.production(k) + (1 - self.delta) * k)

		res = optimize.fminbound(lambda kprime: obj(kprime, t), minbnd, maxbnd)

		return res

# class parametricDP:

# 	def __init__(self, maximize, fitting, calcError):
# 		# Optimization Parameters
# 		self.tol = 0.01
# 		self.error = self.tol + 1
# 		self.maxIt = 300
# 		self.n = 99
# 		self.deg = 3

# 		# Model Parameters
# 		self.alpha = 0.33
# 		self.beta = 0.95
# 		self.depr = 0.1
# 		self.s = 2
# 		self.min = 1.0
# 		self.max = 7.0

# 		self.a = np.zeros(self.deg + 1)
# 		self.v = np.zeros(self.n)

# 		self.grid = np.arange(self.min, self.max, (self.max - self.min) / self.n)

# 		self.maximize = maximize
# 		self.fitting = fitting
# 		self.calcError = calcError

# 		# compute initial guess
# 		# guess = no savings
# 		guessF = lambda k: (1 / (1 - self.s)) * (k**self.alpha - self.depr * k) ** (1 - self.s)

# 		for i in xrange(self.n):
# 			# no savings - kprime = k0
# 			self.v[i] = guessF(self.grid[i])

# 		self.a = polynomial.polyfit(self.grid, self.v, self.deg)

# 		print "INITIAL VALUE FUNCTION",self.v
# 		print "INITIAL APPROX", self.a

# 	def solve(self):
# 		# Solve the dynamic programming problem
# 		# described in parametricDP

# 		it = 0

# 		while self.error > self.tol and it < self.maxIt:
# 			# Save previous a
# 			a0 = np.copy(self.a)

# 			# Maximize
# 			res = self.maximize(self)
# 			self.v = res[0]
# 			self.policy = res[1]

# 			# Fit
# 			self.a = self.fitting(self)

# 			# Update error
# 			self.error = self.calcError(self, self.a, a0)

# 			print "Iteration", it

# 			#print "NEW VALUE FUNCTION", self.v
# 			#print "NEW APPROX", self.a
			
# 			print "Error", self.error
# 			print "\n"

# 			it = it + 1

# 		return (self.a,self.v)