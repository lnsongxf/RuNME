# teststochdp.py

from stochdp import *
import numpy as np
import matplotlib.pyplot as ppt

# Just some code to test my useful functions

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

# Stochastic Wages States
# computes the wage given the current state
def statewage(state, wage):
	if state == 0:
		return wage * 0.9
	elif state == 1:
		return wage * 1
	elif state == 2:
		return wage * 1.1

opts = {}

## Model Parameters

# Utility and Production Functions
opts['utility'] = utilityDefault
opts['bequest'] = bequestValue
opts['production'] = productionDefault
opts['beta'] = 0.9

# Mean wages over time
opts['wages'] = wages

# Number of periods
opts['T'] = T

# Number states describing
# change in wages
opts['states'] = [0,1,2]
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
r = execute(deg = 5, pts=100, opts=opts)[0]
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
# grid = nodes(0.1,10,40)
# yactual = map(u, grid)
# yapprox = map(lambda x: chebval(x, r[0][0], 0.1, 10), grid)
yapprox1 = map(lambda x: chebval(x, r[0][1], 0.1, 10), grid)
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

ppt.clf()
ppt.plot(grid,yapprox1)
ppt.axis([0,10,0,20])
ppt.xlabel('k0')
ppt.ylabel('v')
ppt.savefig('myplot')

# ppt.plot(grid,yapprox)
# ppt.plot(grid,yapprox2)

# ppt.plot(grid,yapprox1s)
# ppt.plot(grid,yapproxs)
# ppt.plot(grid,yapprox2s)

# # ppt.plot(grid,yactual,'k')
# ppt.show()