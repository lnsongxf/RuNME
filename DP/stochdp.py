import numpy as np
from scipy import optimize
import sys

import matplotlib.pyplot as ppt

#from chebfun import *

#import numpy as np

## Runnan Yang 2013
## Solves Life-cycle dynamic programming problems
## by approximating each period's value function
## with a Chebyshev polynomial and performing
## backwards induction from the final to first period


##############################################################
## Convenient functions for working with Chebyshev polynomials
##############################################################

## nodes(init, end, deg)
# Computes Chebyshev nodes over interval [init,end]
#
# Inputs: init - float, beginning of range
#         end  - float, end of range
#         deg  - int, number of points to find
# Output: deg-length array of floats, Chebyshev nodes within the range [init,end]

def nodes(init, end, deg):
    x = np.polynomial.chebyshev.chebgauss(deg)
    return (x[0] + 1.) * (end - init) / 2. + init

## normalize(x, init, end)
# Normalizes x in [init,end] to [init2,end2]
#
# Inputs: x    - float, point to be normalized
#         init - float, beginning of range
#         end  - float, end of range
# Output: float, x normalized to range [init2,end2]

def normalize(x, init, end, init2=-1, end2=1):
	return (float(x) - init) / (end - init) * (end2 - init2) + init2

## chebval(x, c, init, end)
# Evaluates the chebyshev polynomial at point x, after
# normalizing x to [-1,1] from [init,end]
# 
# Inputs: x    - float, point to evaluation the polynomial at
#         c    - array of floats, coefficients of chebyshev polynomial
#         init - float, beginning of range
#         end  - float, end of range
# Output: float, c evaluated at point x

def chebval(x, c, init, end):
	return np.polynomial.chebyshev.chebval(normalize(x,init,end),c)

## plotpoly(grid, poly, init, end)
# Plots a chebyshev polynomial over some set of points grid
# 
# Inputs: grid - array of floats, points to plot the polynomial at
#         c    - array of floats, coefficients of chebyshev polynomial
#         init - float, beginning of range
#         end  - float, end of range
# Output: float, c evaluated at point x 

def plotpoly(grid, poly, init, end):
	vi = map(lambda x: chebval(x, poly, init, end), grid)
	ppt.clf()
	ppt.plot(grid,vi)
	ppt.show()

##############################################################

##############################################################

##############################################################
## Approximation Functions
##############################################################

# Approximates a function whose values at the points xi
# are vi, where xi are in range [init,end]
# Returns a list of coefficients for a 
# Chebyshev polynomial
# Status: DONE
def approxCheb(xi, vi, init, end, deg):
	return np.polynomial.chebyshev.Chebyshev.fit(xi, vi, deg, domain=[init,end]).coef

# Approximates a function whose values at the points xi
# are vi, where xi are in range [init, end], while
# ensuring that the function is monotonic increasing and
# concave
# 
# Return a set of coefficients for a 
# Chebyshev polynomial
# Status: DONE

# Set up the following minimization problem 
# min_{c_j} \sum_{j=0}^n (c_j T_j(x_i) - v_i)**2
# s.t. \sum_{j=0}^n c_j T'_j(y_i) > 0 for i in (1,m') # monotonicity
#      \sum_{j=0}^n c_j T''_j(y_i) < 0 for i in (1,m') # concavity
# where the points y_i are shapePts, and the points x_i
# are the chebyshev interpolation nodes

# Note: I use a number of nested list comprehensions in this
# function, since they are well-optimized in Python, at the risk
# of loss of clarity

def approxChebshape(xi, vi, init, end, deg, shapeGrid=None):
	
	# If grid of shape-preservation points not supplied,
	# generate points uniformly on [init,end]
	if shapeGrid == None:
		shapeGrid = np.arange(init, end, 0.1)

	n = deg + 1
	m = len(shapeGrid)
	nxi = len(xi)
		
	# We will need to first compute T'_j(y_i), T''_j(y_i), and T_j(z_i)
	# given these values, we can set up the constraints for the
	# optimization problem

	# Compute the basis polynomials and the first
	# and second derivatives
	Tj = [np.polynomial.chebyshev.Chebyshev.basis(i, domain=[0,10]) for i in range(0,n+1)]
	Tjd = [p.deriv() for p in Tj]
	Tjdd = [p.deriv(2) for p in Tj]

	# Compute the values of the basis at each point
	Tjx = [[chebval(xi[i], p.coef, init, end) for p in Tj] for i in range(nxi)]
	Tjyd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjd] for i in range(m)]
	Tjydd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjdd] for i in range(m)]

	## Initial Guess
	# Run initial pass with coefficients from
	# fitting without shape preservation

	# Note: in cases where this function is relevant,
	# this polynomial will NOT obey the shape constraints
	# will it be an issue that the shape constraints are
	# violated by the initial guess, or is it more important
	# that the initial guess is good?

	poly0 = approxCheb(xi, vi, init, end, deg)
	
	## Set up objective function for minimization problem

	# Minimize sum of squared errors
	def obj(c):
		return sum([(sum([Tjx[i][j] * c[j] for j in range(n)]) - vi[i])**2. for i in range(nxi)])

	## Define constraints

	# returns an array where each element
	# must be >= 0
	def cons(c):
		constraints = np.zeros(2*m)

		# Monotonocity
		constraints[0:m] = [sum([Tjyd[i][j] * c[j] for j in range(n)]) for i in range(m)]

		# Concavity
		constraints[m:2*m] = [sum([-Tjydd[i][j] * c[j] for j in range(n)]) for i in range(m)]

		return constraints

	# Solve Optimization Problem
	res = optimize.fmin_slsqp(func=obj, 
							  x0=poly0, 
							  f_ieqcons=cons, 
							  iprint=0, 
							  full_output=1)

	# Check optimization status
	if res[3] != 0:
		print "OPTIMIZATION FAILED"
		print res[4] # String describing exit status
		sys.exit()

	return res[0]

##############################################################

##############################################################

##############################################################
## Maximization Functions
##############################################################

# Compute the maximum at each point in xi given that vt
# is a chebyshev polynomial which approximates the
# next period value function

# vt is a list of chebyshev polynomials, s.t.
# vt[0] approximates the high state in next period
# vt[1] approximates the mean state in next period
# vt[2] approximates the low  state in next period

def maxim(xi, vt, state=[0,0], opts=None):

	assert 'utility' in opts, "No Utility Function Specified"
	assert 'production' in opts, "No Production Function Specified"

	wage = state[1]
	state = state[0]

	# Initialize arrays
	n = len(xi)          # num grid points
	vi = np.zeros(n)     # value
	ui = np.zeros((n,2)) # policy

	# Transition matrix between states
	opts['trans'] = opts['trans'] if 'trans' in opts else np.ones(1)

	# Stochastic Wages States
	# computes the wage given the current state
	def statewage(state, wage):
		if state == 0:
			return wage * 0.9
		elif state == 1:
			return wage * 1
		elif state == 2:
			return wage * 1.1

	wage = statewage(state, wage)

	
	def obj(k1,k0):
		# compute current period
		curr = opts['utility'](k1,[k0,wage])

		# call value function approximation for
		# continuation value
		cont = [chebval(k1[0],vt[s],opts['init'],opts['end']) for s in opts['states']]

		# and then take the conditional mean of the value function
		# at each state
		cont = opts['trans'][state,:].dot(cont)

		return -(curr + opts['beta'] * cont)

	for i in xrange(len(xi)):
		k0 = xi[i]
		# no borrowing
		kmin = 0 
		# max possible savings given l=1
		kmax = opts['production'](k0,0) + wage

		lmin = 0
		lmax = 1

		# initial guess: midpoint of max and min
		# Note: can we guarantee within constraints?
		x0 = np.array([(kmin + kmax) / 2.,(lmin + lmax) / 2.])

		options = {'maxiter':1000}

		res = optimize.minimize(fun=lambda x: obj(x,k0), 
								x0=x0,
								method="SLSQP",
								bounds=((kmin,kmax),(lmin, lmax)),
								options=options)

		# Check optimization status
		if res.success == False:
			# If optimization fails, shut down
			print xi[i], k0
			print res
			sys.exit()

		# compute the value of the optimal policy
		vi[i] = -obj(res.x, k0)
		ui[i,:] = res.x

		#print k0, np.round(res.x[0],5), np.round(res.x[1],5), -obj(res.x, k0)

	return [vi,ui]

##############################################################

##############################################################

##############################################################
## Execution Function
##############################################################

## checkOptions(opts)
#
# Checks validity of options and provides default values
# if necessary
#
# Inputs: opts - dictionary
# Output: opts - dictionary

def checkOptions(opts):

	# Assertions about our options dictionary
	# for items where there is no default value

	assert 'utility' in opts, "No Utility Function Specified"
	assert 'bequest' in opts, "No Value for Final Period Specified"
	assert 'production' in opts, "No Production Function Specified"

	assert 'beta' in opts, "No beta value Specified"
	assert 'T' in opts, "No number of periods Specified"

	assert 'init' in opts, "No initial value for Range"
	assert 'end' in opts, "No final value for Range"

	# Variables with default values if not supplied

	opts['states'] = opts['states'] if 'states' in opts else [0]
	opts['wages'] = opts['wages'] if 'wages' in opts else np.zeros(opts['T'])
	opts['trans'] = opts['trans'] if 'trans' in opts else np.ones(1)

	assert np.array(opts['trans']).shape == (len(opts['states']), len(opts['states'])), "Transition matrix malformed."
	assert len(opts['wages']) == opts['T'], "Wages malformed"

	return opts	

# Computes the value of a dynamic programming problem
# over T periods, given range [init,end]
def execute(deg, pts, opts, preserveShape=False):

	##########
	## Processing options
	##########

	opts = checkOptions(opts)

	## Shape-Preservation Option
	# Choose whether to try to preserve shape during the
	# chebyshev approximation step, i.e. which
	# approximation function to call
	if preserveShape:
		approxshape = approxChebshape
	else:
		approxshape = approxCheb

	## Set up grid for approximation
	grid = nodes(opts['init'], opts['end'], pts)

	# List of value and policy function approximations
	value = []
	policy = []

	## Terminal Period
	# Handle last time period separately due to 
	# needing to accomodate bequest

	# Initialize final time period value function
	print "COMPUTING PERIOD", opts['T']-1

	bequest = [opts['bequest'](x) for x in grid]

	poly = [approxshape(grid, bequest, opts['init'], opts['end'], deg) for state in opts['states']]
	value.append(poly)

	# compute second to last period value function points
	vi, ui = zip(*[maxim(grid, poly, state=[state, opts['wages'][-1]], opts=opts) for state in opts['states']])
	policy.append(ui)

	## Iterate back through previous periods
	# Set T=1 for single period
	for t in reversed(range(opts['T'] - 1)):
		print "\nCOMPUTING PERIOD", t

		# approximate next period value function
		poly = [approxshape(grid, vi[state], opts['init'], opts['end'], deg) for state in opts['states']]
		value.append(poly)

		# solve for value function at points
		vi, ui = zip(*[maxim(grid, poly, state=[state, opts['wages'][t]], opts=opts) for state in opts['states']])
		policy.append(ui)

	# Approximate final polynomial (for first period)
	poly = [approxshape(grid, vi[state], opts['init'], opts['end'], deg) for state in opts['states']]
	value.append(poly)

	# Reverse our record of the polynomial approximations
	# so that rec[i] = value function for period i
	# and rec[T] = function for bequeath motive
	value.reverse()
	policy.reverse()

	# Note: value here is a list of coefficients, 
	# policy is a list of values at each grid pt
	return (value, policy)
