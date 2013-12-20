import numpy as np
from scipy import optimize
import sys

import pulp
import matplotlib.pyplot as ppt

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

# Normalizes x in [init,end] to [init2,end2]
def normalize(x, init, end):
	return (float(x) - init) / (end - init) * 2. - 1.

def unnormalize(x, init, end):
	return (float(x) + 1.) / (2.) * (end - init) + init

# Evaluates the chebyshev polynomial at point x, after
# normalizing it to range (-1,1)
def chebval(x, c, init, end):
	return np.polynomial.chebyshev.chebval(normalize(x,init,end),c)

# Plots a chebyshev polynomial over some set of points grid
def plotpoly(grid, poly, init, end):
	vi = map(lambda x: chebval(x, poly, init, end), grid)
	ppt.clf()
	ppt.plot(grid,vi)
	ppt.show()

# Approximates a function whose values at the points xi
# are vi, where xi are in range [init,end]
# Returns a chebyshev polynomial object
# Status: DONE
def approx(xi, vi, init, end, deg):
	# Compress the grid of xi points into [-1,1]
	normxi = np.array(map(lambda x: normalize(x, init, end), xi))

	# Fits a chebyshev polynomial to these normalized points
	poly = np.polynomial.chebyshev.chebfit(normxi, vi, deg)
	#poly = np.polynomial.chebyshev.Chebyshev.fit(xi, vi, deg, domain=[init,end])
	return poly

# Approximates a function whose values at the points xi
# are vi, where xi are in range [init, end]
# 
# Return a chebyshev polynomial object
# Status: IN PROGRESS
def approxshape(xi, vi, init, end, deg, shapePts=None):
	n = deg + 1

	if shapePts == None:
		# Generate points at which to enforce shape
		# uniformly on (init,end)
		shapeGrid = np.arange(init, end, 0.1)
	else:
		shapeGrid = shapePts

	# Set up linear programming problem as described in 
	# Judd Cai 2011 "Shape-Preserving Dynamic Programming"
	# min_{c_j} \sum_{j=0}^n (1/(j+1)^2) (c_j^+ + c_j^-)
	# s.t. \sum_{j=0}^n c_j T'_j(y_i) > 0 for i in (1,m') # monotonicity
	#      \sum_{j=0}^n c_j T''_j(y_i) < 0 for i in (1,m') # concavity
	#      \sum_{j=0}^n c_j T_j(z_i) = v_i for i in (1,m) # concavity
	#      c_j = c_j^+ + c_j^-
	#      c_j^+ >= 0, c_j^- >= 0
	# where the points y_i are shapePts, the points z_i are xi, v_i is vi

	# we will need to first compute T'_j(y_i), T''_j(y_i), and T_j(z_i)
	# given these values, we can set up the linear programming problem

	# Create a chebyshev polynomial of degree deg
	# Compute the basis polynomials
	Tj = [np.polynomial.chebyshev.Chebyshev.basis(i, domain=[0,10]) for i in range(0,n+1)]
	Tjd = [p.deriv() for p in Tj]
	Tjdd = [p.deriv(2) for p in Tj]

	# Compute the values of the basis at each point
	Tjx = [[chebval(xi[i], p.coef, init, end) for p in Tj] for i in range(len(xi))]
	Tjyd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjd] for i in range(len(shapeGrid))]
	Tjydd = [[chebval(shapeGrid[i], p.coef, init, end) for p in Tjdd] for i in range(len(shapeGrid))]

	# Run initial pass with coefficients from
	# fitting without shape preservation
	poly0 = approx(xi, vi, init, end, deg)

	# Set up minimization problem

	def obj(c):
		#minimize least squared errors

		# c is a list of weights on the 
		# corresponding degree j of c[j]

		# compute estimate and subtract 
		sumse = 0.
		#print len(xi), n
		for i in range(len(xi)):
			est = 0.
			for j in range(n):
				est += Tjx[i][j] * c[j]
			sumse += (est - vi[i])**2.

		return sumse

	#print obj([2.058001,1.2780141,-0.2078,0.0752])

	# define constraints

	# Monotonocity
	# return an array where each element
	# should be > 0
	def monoCons(c):
		# for each point
		s = np.zeros(len(shapeGrid))
		for i in range(len(shapeGrid)):
			est = 0.
			# for each degree
			for j in range(n):
				#print Tjyd[i]
				est += Tjyd[i][j] * c[j]
			s[i] = est
		return s

	# Concavity
	# returns an array where each element 
	# should be > 0
	def concCons(c):
		s = np.zeros(len(shapeGrid))
		for i in range(len(shapeGrid)):
			est = 0.
			for j in range(n):
				est -= Tjydd[i][j] * c[j]
			s[i] = est
		return s

	def constraint(c):
		#print c
		out = np.zeros(2*len(shapeGrid))
		out[0:len(shapeGrid)] = monoCons(c)
		out[len(shapeGrid):2*len(shapeGrid)] = concCons(c)

		return out
	#print "poly0", poly0
	res = optimize.fmin_slsqp(func=obj, x0=poly0, f_ieqcons=constraint)
	return res

grid = np.arange(0.1,10,0.01)
#print approx(xi=grid, vi=grid**0.5, init=0.1, end=10, deg=3)
#p = approxshape(xi=grid, vi=grid**0.5, init=0.1, end=10, deg=3)

# Compute the maximum at each point in xi given that vt
# is a chebyshev polynomial which approximates the
# next period value function

# vt is a list of chebyshev polynomials, s.t.
# vt[0] approximates the high state in next period
# vt[1] approximates the mean state in next period
# vt[2] approximates the low  state in next period

def maxim(xi, vt, init, end, state=1, wages=None):
	print "STATE",state
	n = len(xi)
	vi = np.zeros(n)
	u = np.zeros((n,2))
	beta = 0.9
	if wages == None:
		w = 0.
	else:
		w = wages

	mTrans = np.zeros((3,3))
	mTrans[0,:] = [0.5,0.4,0.1]
	mTrans[1,:] = [0.25,.5,.25]
	mTrans[2,:] = [0.1,0.4,0.5]	

	# Stochastic Wages States
	if state == 0:
		w = w * 0.9
	elif state == 1:
		w = w * 1
	elif state == 2:
		w = w * 1.2

	def utility(k1, l, k0):
		# Make sure we don't accidentally try to take
		# the square root of a negative
		if l >= 1:
			return -1
		if k0**0.33 + (0.9 * k0) - k1 + (w * l) < 0:
			return -1

		return (k0**0.33 + (0.9 * k0) - k1 + (w * l))**0.5 + (1. - l)**0.5

	def obj(k1,k0):
		curr = utility(k1[0],k1[1],k0)

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

	return vi

# Computes the value of a dynamic programming problem
# over T periods, given range [init,end]
def execute(init, end, deg, pts, T, preserveShape=False):
	# initialize grid
	grid = nodes(init, end, pts)

	# Extension: wages changing over time
	# Case 1: deterministic wages
	# Wages rise in the early period, remain high in the second period
	# and declines in the third period
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
	#wages = np.zeros(T)

	# Extension: wages that move through different states
	# Status: 3 state model done

	# we will save each of the value function approximation
	# polynomials in rec
	rec = []

	## Handle last time period separately due to 
	# different final value

	# Initialize final time period value function
	# (no states in final period)
	print "PERIOD", T-1
	poly = []
	poly.append(approxshape(grid, grid**0.5, init, end, deg))
	poly.append(approxshape(grid, grid**0.5, init, end, deg))
	poly.append(approxshape(grid, grid**0.5, init, end, deg))
	
	rec.append(poly)
	# compute second to last period value function points	
	vi = []
	vi.append(maxim(grid, poly, init, end, state=0, wages=wages[-1]))
	vi.append(maxim(grid, poly, init, end, state=1, wages=wages[-1]))
	vi.append(maxim(grid, poly, init, end, state=2, wages=wages[-1]))

	## Iterate back through previous periods
	# Set T=1 for single period
	for t in reversed(range(T - 1)):
		print "\nPERIOD", t
		# approximate next period value function
		poly[0] = approxshape(grid, vi[0], init, end, deg)
		poly[1] = approxshape(grid, vi[1], init, end, deg)
		poly[2] = approxshape(grid, vi[2], init, end, deg)
		rec.append(poly)

		# solve for value function at points
		vi[0] = maxim(grid, poly, init, end, state=0, wages=wages[t])
		vi[1] = maxim(grid, poly, init, end, state=1, wages=wages[t])
		vi[2] = maxim(grid, poly, init, end, state=2, wages=wages[t])

	# Approximate final polynomial (for first period)
	poly[0] = approxshape(grid, vi[0], init, end, deg)
	poly[1] = approxshape(grid, vi[1], init, end, deg)
	poly[2] = approxshape(grid, vi[2], init, end, deg)
	rec.append(poly)

	rec.reverse()
	# rec[i] = value function for period i
	# from 0 through T (period T = bequeath)

	return rec

# Just some code to test my useful functions

r = execute(init=0.1, end=10., deg = 10, pts=40, T=1)
