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
# Returns a chebyshev polynomial object
# Status: DONE
def approx(xi, vi, init, end, deg):
	return np.polynomial.chebyshev.Chebyshev.fit(xi, vi, deg, domain=[init,end])

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

def approxshape(xi, vi, init, end, deg, shapeGrid=None):
	n = deg + 1

	# If grid of shape-preservation points not supplied,
	# generate points uniformly on [init,end]
	if shapeGrid == None:
		shapeGrid = np.arange(init, end, 0.1)
		
	# We will need to first compute T'_j(y_i), T''_j(y_i), and T_j(z_i)
	# given these values, we can set up the constraints for the
	# optimization problem

	# Compute the basis polynomials and the first
	# and second derivatives
	Tj = [np.polynomial.chebyshev.Chebyshev.basis(i, domain=[0,10]) for i in range(0,n+1)]
	Tjd = [p.deriv() for p in Tj]
	Tjdd = [p.deriv(2) for p in Tj]

	# Compute the values of the basis at each point
	Tjx = [[chebval(xi[i], p.coef, init, end) for p in Tj] for i in range(len(xi))]
	Tjyd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjd] for i in range(len(shapeGrid))]
	Tjydd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjdd] for i in range(len(shapeGrid))]

	## Initial Guess
	# Run initial pass with coefficients from
	# fitting without shape preservation

	# Observation: in cases where this function is relevant,
	# this polynomial will NOT obey the shape constraints
	# will it be an issue that the shape constraints are
	# violated by the initial guess, or is it more important
	# that the initial guess is good?
	poly0 = approx(xi, vi, init, end, deg)
	
	## Set up objective function for minimization problem

	# Minimize sum of squared errors
	def obj(c):
		return sum([(sum([Tjx[i][j] * c[j] for j in range(n)]) - vi[i])**2. for i in range(len(xi))])

	## Define constraints

	# returns an array where each element
	# must be >= 0
	def cons(c):
		constraints = np.zeros(2*len(shapeGrid))

		# Monotonocity
		constraints[0:len(shapeGrid)] = [sum([Tjyd[i][j] * c[j] for j in range(n)]) for i in range(len(shapeGrid))]

		# Concavity
		constraints[len(shapeGrid):2*len(shapeGrid)] = [sum([-Tjydd[i][j] * c[j] for j in range(n)]) for i in range(len(shapeGrid))]

		return constraints

	# Solve Optimization Problem
	res = optimize.fmin_slsqp(func=obj, x0=poly0.coef, f_ieqcons=cons, iprint=0, full_output=1)

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

def maxim(xi, vt, init, end, utility=None, state=1, wage=0):
	# print "STATE",state

	# Initialize arrays
	n = len(xi)
	vi = np.zeros(n)
	u = np.zeros((n,2))

	beta = 0.9

	mTrans = np.zeros((3,3))
	mTrans[0,:] = [0.75,0.2,0.05]
	mTrans[1,:] = [0.25,.5,.25]
	mTrans[2,:] = [0.05,0.2,0.75]

	# Stochastic Wages States
	if state == 0:
		wage = wage * 0.9
	elif state == 1:
		wage = wage * 1
	elif state == 2:
		wage = wage * 1.1

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



	if utility == None:
		print "ERROR: No Utility Function Supplied"
		print "Exiting function"
		sys.exit()

	def obj(k1,k0):
		curr = utility(k1,[k0,wage])

		# call value function approximation for
		# continuation value
		# ***for each state***
		cont0 = chebval(k1[0],vt[0],init,end)
		cont1 = chebval(k1[0],vt[1],init,end)
		cont2 = chebval(k1[0],vt[2],init,end)

		# and then take the conditional mean of the value function
		# at each state
		cont = mTrans[state,:].dot([cont0, cont1, cont2])

		return -(curr + beta * cont)

	for i in xrange(len(xi)):
		k0 = xi[i]
		# no borrowing
		kmin = 0 
		# max possible savings given l=1
		kmax = k0**0.33 + (1. - 0.1) * k0 + wage

		# initial guess: midpoint of max and min, work 1/2
		x0 = np.array([(k0**0.33 + 0.9*k0) / 2.,0.5])

		# set constraints
		# k1 < f(k0) + w*l

		cons = {}
		cons['type'] = 'ineq'
		cons['fun'] = lambda x: k0**0.5 + 0.9*k0 + wage*x[1] - x[0]

		options = {}
		options['maxiter'] = 1000

		def obj2(k1):
			return obj(k1, k0)

		res = optimize.minimize(fun=obj2, 
								x0=x0,
								method="SLSQP",
								constraints=cons,
								bounds=((kmin,kmax),(0, 1)),
								options=options)

		# Check optimization status
		if res.success == False:
			# If optimization fails, shut down
			print xi[i], k0
			print res
			sys.exit()

		# compute the value of the optimal policy
		vi[i] = -obj(res.x, k0)
		u[i,:] = res.x

		#print k0, np.round(res.x[0],5), np.round(res.x[1],5), -obj(res.x, k0)

	return vi

# Computes the value of a dynamic programming problem
# over T periods, given range [init,end]
def execute(init, end, deg, pts, T, preserveShape=False, utility=None, wages=None):
	# initialize grid
	grid = nodes(init, end, pts)

	states = [0,1,2]

	if wages==None:
		wages = np.zeros(T)

	# we will save each of the value function approximation
	# polynomials in rec
	rec = []

	## Handle last time period separately due to 
	# different final value

	# Initialize final time period value function
	# (no states in final period)
	print "PERIOD", T-1

	poly = [approxshape(grid, grid**0.5, init, end, deg) for state in states]
	rec.append(poly)

	# compute second to last period value function points
	vi = [maxim(grid, poly, init, end, utility=utility, state=state, wage=wages[-1]) for state in states]

	## Iterate back through previous periods
	# Set T=1 for single period
	for t in reversed(range(T - 1)):
		print "\nPERIOD", t

		# approximate next period value function
		poly = [approxshape(grid, vi[state], init, end, deg) for state in states]
		rec.append(poly)

		# solve for value function at points
		vi = [maxim(grid, poly, init, end, utility=utility, state=state, wage=wages[t]) for state in states]

	# Approximate final polynomial (for first period)
	poly = [approxshape(grid, vi[state], init, end, deg) for state in states]
	rec.append(poly)

	# Reverse our record of the polynomial approximations
	# so that rec[i] = value function for period i
	# and rec[T] = function for bequeath motive
	rec.reverse()

	return rec

# Just some code to test my useful functions

def utilityDefault(x, y):

	k1 = x[0]
	l = x[1]

	k0 = y[0]
	w = y[1]

	# Make sure we don't try to take
	# the square root of a negative
	# (the other bounds can be handled
	#  by the optimizer)
	if l >= 1:
		return -1
	if k0**0.33 + (0.9 * k0) - k1 + (w * l) < 0:
		return -1

	return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5 + (1. - l)**0.5

T = 5

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

r = execute(init=0.1, end=10., deg = 5, pts=40, T=5, utility=utilityDefault, wages=wages)

# two periods with labor choice
# analytic soln

grid = np.arange(0.1,10,0.01)

income = lambda k0: k0**0.33 + 0.9 * k0
labor = lambda k0: (25. + 0.81*25. - income(k0)) / (5. + 25. + 25. * 0.81)
c1 = lambda k0: 25. * (1. - labor(k0))
c2 = lambda k0: 0.81 * c1(k0)
u = lambda k0: c1(k0)**0.5 + 0.9 * c2(k0)**0.5 + (1.-labor(k0))**0.5

grid = np.arange(0.1,10,0.01)
yactual = map(u, grid)
yapprox = map(lambda x: chebval(x, r[0][0], 0.1, 10), grid)
yapprox1 = map(lambda x: chebval(x, r[0][1], 0.1, 10), grid)
yapprox2 = map(lambda x: chebval(x, r[0][2], 0.1, 10), grid)
ppt.clf()
ppt.plot(grid,yapprox1)
ppt.plot(grid,yapprox)
ppt.plot(grid,yapprox2)
ppt.plot(grid,yactual,'k')
ppt.show()