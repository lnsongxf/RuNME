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

# Computes Chebyshev nodes over interval [init,end]
def nodes(init, end, deg):
    x = np.polynomial.chebyshev.chebgauss(deg)
    return (x[0] + 1.) * (end - init) / 2. + init

# Normalizes x in [init,end] to [init2,end2]
def normalize(x, init, end):
	return (float(x) - init) / (end - init) * 2. - 1.

def unnormalize(x, init, end):
	return (float(x) + 1.) / (2.) * (end - init) + init

# Evaluates the chebyshev polynomial at point x, after
# normalizing x to [-1,1] from [init,end]
def chebval(x, c, init, end):
	return np.polynomial.chebyshev.chebval(normalize(x,init,end),c)

# Plots a chebyshev polynomial over some set of points grid
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
		cont = mTrans[state,:].dot(cont)

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
## Maximization Functions
##############################################################

# Computes the value of a dynamic programming problem
# over T periods, given range [init,end]
def execute(deg, pts, opts, preserveShape=False):

	# Processing options

	#####
	# the utility function supplied should
	# take two arguments, x and k0
	# where x is a k-length array of the controls
	# and k0 is the single state variable
	#
	# for example, if we were solving a problem
	# where an agent chooses savings and labor,
	# we might let x[0] = k_{t+1} and x[1] = l
	#####

	# Big list of assertions about our
	# options dictionary

	assert 'utility' in opts, "No Utility Function Specified"
	assert 'bequest' in opts, "No Value for Final Period Specified"
	assert 'init' in opts, "No initial value for Range"
	assert 'end' in opts, "No final value for Range"

	init = opts['init']
	end = opts['end']

	opts['states'] = opts['states'] if 'states' in opts else [0]
	opts['wages'] = opts['wages'] if 'wages' in opts else np.zeros(T)

	# initialize grid
	grid = nodes(opts['init'], opts['end'], pts)

	# Choose whether to try to preserve shape during the
	# chebyshev approximation step
	if preserveShape:
		approxshape = approxChebshape
	else:
		approxshape = approxCheb

	# we will save each of the value function approximation
	# polynomials in rec
	value = []
	policy = []

	## Handle last time period separately due to 
	# different final value

	# Initialize final time period value function
	print "PERIOD", opts['T']-1

	bequest = [opts['bequest'](x) for x in grid]

	poly = [approxshape(grid, bequest, opts['init'], opts['end'], deg) for state in opts['states']]
	value.append(poly)

	# compute second to last period value function points
	vi = [maxim(grid, poly, state=[state, wages[-1]], opts=opts) for state in opts['states']]
	ui = [vi[state][1] for state in opts['states']]
	vi = [vi[state][0] for state in opts['states']]
	policy.append(ui)

	## Iterate back through previous periods
	# Set T=1 for single period
	for t in reversed(range(opts['T'] - 1)):
		print "\nPERIOD", t

		# approximate next period value function
		poly = [approxshape(grid, vi[state], opts['init'], opts['end'], deg) for state in opts['states']]
		value.append(poly)

		# solve for value function at points
		vi = [maxim(grid, poly, state=[state, wages[t]], opts=opts) for state in opts['states']]
		ui = [vi[state][1] for state in opts['states']]
		vi = [vi[state][0] for state in opts['states']]
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

# Just some code to test my useful functions

# x - controls
# y - state variables
def utilityDefault(x, y):

	k1 = x[0]
	l = x[1]

	k0 = y[0]
	w = y[1]

	# Make sure we don't try to take
	# the square root of a negative

	if l >= 1:
		return -1
	if k0**0.33 + (0.9 * k0) - k1 + (w * l) < 0:
		return -1

	return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5 + (1. - l)**0.5

def bequestValue(x):
	return x**0.5

def productionDefault(k0, l):
		return k0 ** 0.5 + (1 - 0.1) * k0

T = 10

wages = np.zeros(T)
period = int(np.floor(T/3.))
if period == 0:
	print "Too few periods - setting constant wages"
	wages[:] = 5
else:
	for i in xrange(period+1):
		wages[i] = i * (5. / period)# growing up to 5
	wages[period:2*period] = 5.
	for i in xrange(2*period,T):
		wages[i] = 5. - 2. * (i - 2.*period) / (T - 2.*period)

mTrans = np.zeros((3,3))
mTrans[0,:] = [0.75,0.2,0.05]
mTrans[1,:] = [0.25,.5,.25]
mTrans[2,:] = [0.05,0.2,0.75]

opts = {}

## Model Parameters

opts['init'] = 0.1
opts['end'] = 10.
opts['wages'] = wages
opts['utility'] = utilityDefault
opts['bequest'] = bequestValue
opts['production'] = productionDefault
opts['T'] = T
opts['states'] = [0,1,2]
opts['beta'] = 0.9
opts['trans'] = mTrans

## Optimization Parameters

# Degree of Chebyshev Polynomial
# and number of nodes
opts['deg'] = 5.
opts['pts'] = 40.

# Range of points
opts['init'] = 0.1
opts['end'] = 10.

## Execute the solver
r = execute(deg = 5, pts=40, opts=opts)
#r2 = execute(deg = 5, pts=40, opts=opts, preserveShape=True)

# two periods with labor choice
# analytic soln

grid = np.arange(0.1,10,0.01)

income = lambda k0: k0**0.33 + 0.9 * k0
labor = lambda k0: (25. + 0.81*25. - income(k0)) / (5. + 25. + 25. * 0.81)
c1 = lambda k0: 25. * (1. - labor(k0))
c2 = lambda k0: 0.81 * c1(k0)
u = lambda k0: c1(k0)**0.5 + 0.9 * c2(k0)**0.5 + (1.-labor(k0))**0.5

grid = np.arange(0.1,10,0.01)
grid = nodes(0.1,10,40)
yactual = map(u, grid)
# yapprox = map(lambda x: chebval(x, r[0][0], 0.1, 10), grid)
# yapprox1 = map(lambda x: chebval(x, r[0][1], 0.1, 10), grid)
# yapprox2 = map(lambda x: chebval(x, r[0][2], 0.1, 10), grid)

# yapproxs = map(lambda x: chebval(x, r2[0][0], 0.1, 10), grid)
# yapprox1s = map(lambda x: chebval(x, r2[0][1], 0.1, 10), grid)
# yapprox2s = map(lambda x: chebval(x, r2[0][2], 0.1, 10), grid)

# yapprox = r[0][1]
# yapprox1 = r[30][1]
# yapprox2 = r[59][1]

# yapproxs = r2[0][1]
# yapprox1s = r2[30][1]
# yapprox2s = r2[59][1]

# print len(grid), len(yapprox1)

# ppt.clf()
# ppt.plot(grid,yapprox1)
# ppt.plot(grid,yapprox)
# ppt.plot(grid,yapprox2)

# ppt.plot(grid,yapprox1s)
# ppt.plot(grid,yapproxs)
# ppt.plot(grid,yapprox2s)

# # ppt.plot(grid,yactual,'k')
# ppt.show()