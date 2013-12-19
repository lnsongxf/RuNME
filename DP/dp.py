import numpy as np
from scipy import optimize
import sys

## Runnan Yang 2013
## Solves Life-cycle dynamic programming problems
## by approximating each period's value function
## with a Chebyshev polynomial and performing
## backwards induction from the final to first period

## Convenient functions for working with Chebyshev polynomials

# Computes Chebyshev nodes over interval [init,end]
def nodes(init, end, deg):
    x = np.polynomial.chebyshev.chebgauss(deg)
    return (x[0] + 1.) * (end - init) / 2. + init

# Normalizes x in [init,end] to [-1,1]
def normalize(x, init, end):
	return (float(x) - init) / (end - init) * 2. - 1.

def unnormalize(x, init, end):
	return (float(x) + 1.) / (2) * (end - init) + init

def chebval(x, c, init, end):
	return np.polynomial.chebyshev.chebval(normalize(x,init,end),c)

# Approximates a function whose values at the points xi
# are vi, where xi are in range [init,end]
# Returns a chebyshev polynomial object
def approx(xi, vi, init, end, deg):
	# Compress the grid of xi points into [-1,1]
	normxi = np.array(map(lambda x: normalize(x, init, end), xi))

	# Fits a chebyshev polynomial to these normalized points
	poly = np.polynomial.chebyshev.chebfit(normxi, vi, deg)
	return poly

# Compute the maximum at each point in xi given that vt
# is a chebyshev polynomial which approximates the
# next period value function



def maxim(xi,vt, init, end, state=0, wages=None):
	n = len(xi)
	vi = np.zeros(n)
	u = np.zeros((n,2))
	beta = 0.9
	if wages == None:
		w = 0.
	else:
		w = wages

	# mTrans = np.zeros((3,3))
	# mTrans[0,:] = [0.5,0.4,0.1]
	# mTrans[1,:] = [0.25,.5,.25]
	# mTrans[2,:] = [0.1,0.4,0.5]	

	# Stochastic Wages States
	if state == 0:
		w = w
	elif state == 1:
		w = w * 1.1
	elif state == -1:
		w = w * 0.9

	def utility(k1, l, k0):
		if k1 < 0:
			return -1
		if l < 0:
			return -1
		if l >= 1:
			return -1
		if k0**0.33 + (0.9 * k0) - k1 + (w * l) < 0:
			return -1

		return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5 + (1. - l)**0.5

	def obj(k1, k0):
		curr = utility(k1[0],k1[1],k0)

		# call value function approximation for
		# continuation value
		# ***for each state***
		cont = chebval(k1[0],vt,init,end)
		# and then take the conditional mean of these states

		return -(curr + beta * cont)

	for i in xrange(len(xi)):
		k0 = xi[i]
		kmin = 0 # no borrowing
		kmax = k0**0.33 + (1. - 0.1) * k0 + w

		# initial guess: midpoint of max and min, work 1/2
		x0 = np.array([(k0**0.33 + 0.9*k0) / 2.,0.5])

		# set constraints
		# k1 < f(k0) + w*l

		cons = {}
		cons['type'] = 'ineq'
		cons['fun'] = lambda x: k0**0.5 + 0.9*k0 + w*x[1] - x[0]

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
			# If optimization fails, freak out and die
			print xi[i], k0
			print res
			sys.exit()

		# compute the value of the optimal policy
		vi[i] = -obj(res.x, k0)
		u[i,:] = res.x

		print k0, np.round(res.x[0],5), np.round(res.x[1],5), -obj(res.x, k0)

	#print u
	return vi

# Computes the value of a dynamic programming problem
# over T periods, given range [init,end]
def execute(init, end, deg, pts, T):
	# initialize grid
	grid = nodes(init, end, pts)

	# Extension: wages changing over time
	# Case 1: deterministic wages
	# Wages rise in the early period, remain high in the second period
	# and declines in the third period
	wages = np.zeros(T)
	period = int(np.floor(T/3.))
	for i in xrange(period+1):
		wages[i] = i * (5. / period)# growing up to 5
	wages[period:2*period] = 5.
	for i in xrange(2*period,T):
		wages[i] = 5. - 2. * (i - 2.*period) / (T - 2.*period)

	# Extension: wages that move through different states

	# we will save each of the value function approximation
	# polynomials in rec
	rec = []

	## Handle last time period separately due to 
	# different final value

	# Initialize final time period value function
	print "PERIOD", T-1
	poly = approx(grid, grid**0.5, init, end, deg)

	rec.append(poly)
	# compute second to last period value function points	
	vi = maxim(grid, poly, init, end, wages[-1])

	## Iterate back through previous periods
	# Set T=1 for single period
	for t in reversed(range(T - 1)):
		print "\nPERIOD", t
		# approximate next period value function
		poly = approx(grid, vi, init, end, deg)
		rec.append(poly)

		# solve for value function at points
		vi = maxim(grid, poly, init, end, wages[t])

	poly = approx(grid, vi, init, end, deg)
	rec.append(poly)

	# rec[i] = value function for period T-i
	# rec[0] = value function for bequest
	# rec[1] = value function for final period
	# rec[T] = value function for first period

	return rec

# Just some code to test my useful functions

r = execute(init=0.1, end=10., deg = 5, pts=50, T=5)

# Some code to compute things about a wage process which moves
# over time according to a markov chain
# mTrans = np.zeros((3,3))
# mTrans[0,:] = [0.5,0.4,0.1]
# mTrans[1,:] = [0.25,.5,.25]
# mTrans[2,:] = [0.1,0.4,0.5]

# erg = np.linalg.matrix_power(mTrans,30)

# payoff = np.array([[1.1,1,0.9],]).T

# print mTrans
# print erg
# print payoff
# print np.dot(erg[1,:],payoff)
