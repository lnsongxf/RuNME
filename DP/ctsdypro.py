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

# From Algorithm 12.5 in Judd 1998, we have the following structure:
# Objective: Solve the Bellman equation
# Step 0. Initialization: Choose function form for Vhat and approximation grid X.
#						  Make initial guess for Vhat and choose stopping criterion
# Step 1. Maximization: Compute v_j using TVhat for each x_j in X
# Step 2. Fitting: Compute a^i+1 such that Vhat approximates the data (v_i, x_i)
# Step 3. Check if stopping criterion met, otherwise go to step 1.

class parametricDP:

	def __init__(self, parameters, optParam, initial, maximize, fitting, calcError):
		# Dictionary of parameters of model
		self.parameters = parameters

		self.min = parameters['min']
		self.max = parameters['max']

		# Dictionary of optimization parameters
		self.optParam = optParam

		# Tolerance specifications
		self.tol = optParam['tol']
		self.error = self.tol + 1

		# Polynomial Degree
		self.deg = optParam.deg

		# Grid length
		self.n = optParam['n']

		# ***Value function realizations*** #
		self.v = np.zeros(n)

		# ***Initial coefficients for approximation*** #
		self.a = np.copy(initial)

		# Function for performing maximization step (Step 1)
		self.maximize = maximize

		# Function for performing fitting step (Step 2)
		self.fitting = fitting

		# Function for calculating the error in a
		# i.e. the error in approximation coefficients
		self.calcError = calcError

	def createGrid(self):
		# computes grid of approximation points

		# Uniform Grid
		self.grid = np.arange(self.min, self.max, (self.max - self.min) / self.n)

		# Chebyshev Grid
		# (to implement)

	def solve(self):
		# Solve the dynamic programming problem
		# described in parametricDP

		it = 0

		while self.error > self.tol and it < maxIt:
			# Save previous a
			a0 = np.copy(self.a)

			# Maximize
			self.v = self.maximize(self)

			# Fit
			self.a = self.fitting(self)

			# Update error
			self.error = self.calcError(a, a0)

			# Update a
			a = np.copy(a0)

			it = it + 1

		return (a,v)

# Problem Description
# consumption problem as in dypro.py

optParam = {}
optParam['tol'] = 0.1
optParam['maxIt'] = 10
optParam['n'] = 99.0
optParam['deg'] = 5

param = {}
param['alpha'] = 0.33
param['beta'] = 0.95
param['depr'] = 0.1
param['s'] = 2
param['min'] = 0.0
param['max'] = 7.0

initial = np.zeros(optParam['n'])

# Functions



def maximize(a, s, alpha, beta, depr, grid):
	# Find the value of v that maximizes 
	# given self.a

	# set up objective function

	def utility(consumption):
		return (1 / (1 - s)) * consumption**(1 - s)

	def production(k0):
		return k0 ** alpha + (1 - depr) * k0

	def obj(kprime, k0):
		# objective function for optimization problem

		# Compute current period value
		consumption = production(k0)
		curr = utility(consumption) # consumption in current period

		# Compute continuation value
		# use x_i point		
		cont = 0
		val = polynomial.polyval(kprime, a)

		# Q * vals
		#p2 = 

		return -(curr + beta * cont)

	x0 = 0
	minbds = grid[0]
	maxbds = min(grid[-1], production(k0)**alpha + (1 - depr) * k0
	bounds = ((minbds, maxbds),)

	res = optimize.minimize(obj, x0=, method='L-BFGS-B', bounds=bounds)

	return 0

grid = np.arange(param['min'], param['max'], (param['max'] - param['min']) / optParam['n'])

def fitting(v, grid, deg):
	# Find the value of a that fits 
	# given self.v

	# Least Squares Approximation with polynomials
	coef = polynomial.polyfit(grid, v, deg)

	# Approximation with Chebyshev polynomials
	# TODO

	return coef


def calcError(a, aprime, n, grid):
	# Evaluates the polynomials described by coefficients
	# a and aprime at each point on the grid
	# and calculates the norm of the difference between
	# the two vectors

	diff = np.zeros(n)
	for i in xrange(n):
		# evaluate polynomials
		diff[i] = polynomial.polyval(grid[i], a) - polynomial.polyval(grid[i], aprime)
	
	return np.sqrt(np.dot(diff,diff))